#!/usr/bin/env python3
"""Reference-guided scaffolding for final_finalizer.

Produces chromosome-scale pseudomolecules by:
1. Grouping contigs by assigned reference chromosome
2. Running RagTag scaffold per chromosome (if available)
3. Falling back to built-in reference-position-based ordering otherwise
4. Producing AGP 2.0 output describing scaffold structure
"""
from __future__ import annotations

import shutil
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.models import ContigClassification
from final_finalizer.utils.io_utils import have_exe
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.sequence_utils import (
    read_fasta_sequences,
    reverse_complement,
    write_fasta,
    write_filtered_fasta,
)

logger = get_logger("scaffolding")


# ---------------------------------------------------------------------------
# AGP output
# ---------------------------------------------------------------------------

def write_agp(agp_lines: List[str], output_path: Path) -> None:
    """Write AGP 2.0 file.

    Args:
        agp_lines: List of tab-delimited AGP lines (without trailing newline).
        output_path: Output file path.
    """
    with output_path.open("w") as fh:
        fh.write("##agp-version\t2.0\n")
        for line in agp_lines:
            fh.write(line + "\n")


def _agp_component_line(
    scaffold_name: str,
    scaffold_start: int,
    scaffold_end: int,
    part_number: int,
    component_id: str,
    component_start: int,
    component_end: int,
    orientation: str,
) -> str:
    """Create an AGP W (component) line."""
    return "\t".join(map(str, [
        scaffold_name,
        scaffold_start,
        scaffold_end,
        part_number,
        "W",
        component_id,
        component_start,
        component_end,
        orientation,
    ]))


def _agp_gap_line(
    scaffold_name: str,
    scaffold_start: int,
    scaffold_end: int,
    part_number: int,
    gap_length: int,
) -> str:
    """Create an AGP N (gap) line."""
    return "\t".join(map(str, [
        scaffold_name,
        scaffold_start,
        scaffold_end,
        part_number,
        "N",
        gap_length,
        "scaffold",
        "yes",
        "align_genus",
    ]))


# ---------------------------------------------------------------------------
# Contig grouping
# ---------------------------------------------------------------------------

def _group_contigs_by_ref_chrom(
    classifications: List[ContigClassification],
    best_ref: Dict[str, str],
    contig_refs: Dict[str, Set[str]],
) -> Dict[str, List[str]]:
    """Group contigs by target reference chromosome for scaffolding.

    Chromosome-assigned contigs go to their assigned ref.
    Candidate pool (chrom_debris only) goes to their assigned_ref_id.

    Returns:
        Dict mapping ref_id -> list of contig original names.
    """
    groups: Dict[str, List[str]] = defaultdict(list)

    # First pass: assigned chromosomes
    assigned_contigs: Set[str] = set()
    for clf in classifications:
        if clf.classification == "chrom_assigned" and clf.assigned_ref_id:
            groups[clf.assigned_ref_id].append(clf.original_name)
            assigned_contigs.add(clf.original_name)

    # Second pass: only chrom_debris contigs are scaffolding candidates
    for clf in classifications:
        if clf.original_name in assigned_contigs:
            continue
        if clf.classification != "chrom_debris":
            continue

        ref_id = clf.assigned_ref_id or best_ref.get(clf.original_name, "")
        if not ref_id:
            continue

        # Only include if this ref is already a scaffolding target
        if ref_id in groups:
            groups[ref_id].append(clf.original_name)

    return dict(groups)


def _compute_group_mean_idents(
    groups: Dict[Tuple[str, int], List[str]],
    clf_lookup: Dict[str, ContigClassification],
) -> Dict[Tuple[str, int], float]:
    """Compute mean seq_identity_vs_ref for chrom_assigned contigs in each group.

    Args:
        groups: Dict mapping (ref_id, grp) -> list of contig names.
        clf_lookup: Dict mapping contig name -> ContigClassification.

    Returns:
        Dict mapping (ref_id, grp) -> mean identity.
    """
    result: Dict[Tuple[str, int], float] = {}
    for key, contigs in groups.items():
        idents = []
        for c in contigs:
            clf = clf_lookup.get(c)
            if clf and clf.classification == "chrom_assigned" and clf.seq_identity_vs_ref is not None:
                idents.append(clf.seq_identity_vs_ref)
        if idents:
            result[key] = statistics.mean(idents)
    return result


def _group_contigs_by_haplotype(
    classifications: List[ContigClassification],
    best_ref: Dict[str, str],
    contig_refs: Dict[str, Set[str]],
    qr_best_chain_ident: Dict[Tuple[str, str], float],
) -> Dict[Tuple[str, int], List[str]]:
    """Group contigs by (ref_id, query_subgenome_grp) for haplotype-aware scaffolding.

    Pass 1: Chromosome-assigned contigs grouped by (assigned_ref_id, grp).
    Pass 2: chrom_debris contigs assigned to the nearest haplotype group by identity.

    Args:
        classifications: List of ContigClassification objects.
        best_ref: Dict mapping contig -> best ref_id.
        contig_refs: Dict mapping contig -> set of ref_ids with evidence.
        qr_best_chain_ident: Dict of (contig, ref_id) -> best chain identity.

    Returns:
        Dict mapping (ref_id, grp) -> list of contig original names.
    """
    groups: Dict[Tuple[str, int], List[str]] = defaultdict(list)
    clf_lookup = {clf.original_name: clf for clf in classifications}

    # Pass 1: chromosome-assigned contigs
    assigned_contigs: Set[str] = set()
    for clf in classifications:
        if clf.classification == "chrom_assigned" and clf.assigned_ref_id:
            grp = clf.query_subgenome_grp or 1
            groups[(clf.assigned_ref_id, grp)].append(clf.original_name)
            assigned_contigs.add(clf.original_name)

    # Build index: ref_id -> set of groups present
    ref_groups: Dict[str, Set[int]] = defaultdict(set)
    for (ref_id, grp) in groups:
        ref_groups[ref_id].add(grp)

    # Compute mean identity per group (for debris assignment)
    group_mean_idents = _compute_group_mean_idents(groups, clf_lookup)

    # Pass 2: chrom_debris candidates only
    for clf in classifications:
        if clf.original_name in assigned_contigs:
            continue
        if clf.classification != "chrom_debris":
            continue

        ref_id = clf.assigned_ref_id or best_ref.get(clf.original_name, "")
        if not ref_id:
            continue

        available_grps = ref_groups.get(ref_id)
        if not available_grps:
            continue

        if len(available_grps) == 1:
            # Single haplotype group for this ref - assign directly
            target_grp = next(iter(available_grps))
        else:
            # Multiple groups - assign by closest identity
            debris_ident = qr_best_chain_ident.get((clf.original_name, ref_id), 0.0)
            target_grp = min(
                available_grps,
                key=lambda g: abs(group_mean_idents.get((ref_id, g), 0.0) - debris_ident),
            )

        groups[(ref_id, target_grp)].append(clf.original_name)

    return dict(groups)


# ---------------------------------------------------------------------------
# Built-in fallback scaffolder
# ---------------------------------------------------------------------------

def _builtin_scaffold(
    contig_names: List[str],
    sequences: Dict[str, str],
    orientations: Dict[str, bool],
    ref_ranges: Dict[Tuple[str, str], Tuple[int, int]],
    ref_id: str,
    gap_size: int,
) -> Tuple[str, List[str]]:
    """Built-in scaffolder: sort contigs by reference position midpoint.

    Args:
        contig_names: List of contig names to scaffold.
        sequences: Dict mapping contig name -> DNA sequence.
        orientations: Dict mapping contig name -> True if reversed.
        ref_ranges: Dict mapping (contig, ref_id) -> (ref_min, ref_max).
        ref_id: Reference chromosome ID for this scaffold.
        gap_size: Number of Ns between contigs.

    Returns:
        Tuple of (scaffolded_sequence, agp_lines).
    """
    # Compute midpoint on reference for sorting
    def ref_midpoint(contig: str) -> float:
        rng = ref_ranges.get((contig, ref_id))
        if rng:
            return (rng[0] + rng[1]) / 2.0
        return float("inf")

    sorted_contigs = sorted(contig_names, key=ref_midpoint)

    scaffold_parts: List[str] = []
    agp_lines: List[str] = []
    scaffold_pos = 1  # AGP is 1-based
    part_number = 0

    for i, contig in enumerate(sorted_contigs):
        seq = sequences.get(contig, "")
        if not seq:
            continue

        # Orient contig
        if orientations.get(contig, False):
            seq = reverse_complement(seq)
            orientation = "-"
        else:
            orientation = "+"

        # Add gap before non-first contig
        if scaffold_parts and gap_size > 0:
            part_number += 1
            gap_end = scaffold_pos + gap_size - 1
            agp_lines.append(_agp_gap_line(
                ref_id, scaffold_pos, gap_end, part_number, gap_size,
            ))
            scaffold_parts.append("N" * gap_size)
            scaffold_pos = gap_end + 1

        # Add component
        part_number += 1
        comp_end = scaffold_pos + len(seq) - 1
        agp_lines.append(_agp_component_line(
            ref_id, scaffold_pos, comp_end, part_number,
            contig, 1, len(seq), orientation,
        ))
        scaffold_parts.append(seq)
        scaffold_pos = comp_end + 1

    return "".join(scaffold_parts), agp_lines


# ---------------------------------------------------------------------------
# RagTag integration
# ---------------------------------------------------------------------------

def _run_ragtag_scaffold(
    ref_chr_fasta: Path,
    contigs_fasta: Path,
    work_dir: Path,
    threads: int,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Run ragtag.py scaffold for a single reference chromosome.

    Returns:
        Tuple of (agp_path, fasta_path) or (None, None) on failure.
    """
    # Clean any stale output from previous runs.  RagTag caches intermediate
    # alignments (*.paf) and skips re-alignment when they exist, which causes
    # wrong results when the input contig set has changed.
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "ragtag.py", "scaffold",
        str(ref_chr_fasta),
        str(contigs_fasta),
        "-o", str(work_dir),
        "-t", str(threads),
        "-u",  # Include unplaced contigs
    ]
    logger.info(f"Running RagTag: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            logger.warning(f"RagTag failed (rc={result.returncode}): {result.stderr[:500]}")
            return None, None
    except subprocess.TimeoutExpired:
        logger.warning("RagTag timed out after 3600s")
        return None, None
    except FileNotFoundError:
        logger.warning("ragtag.py not found in PATH")
        return None, None

    agp_path = work_dir / "ragtag.scaffold.agp"
    fasta_path = work_dir / "ragtag.scaffold.fasta"

    if agp_path.exists() and fasta_path.exists():
        return agp_path, fasta_path
    return None, None


def _parse_ragtag_agp(
    agp_path: Path,
    scaffold_name: str,
) -> List[str]:
    """Parse RagTag AGP output and rename scaffold to our naming convention.

    Args:
        agp_path: Path to RagTag's ragtag.scaffold.agp.
        scaffold_name: Desired scaffold name (e.g., 'chr1A').

    Returns:
        List of AGP lines with scaffold renamed.
    """
    agp_lines: List[str] = []
    with agp_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            # Replace RagTag scaffold name with our name
            # RagTag names scaffolds like "refname_RagTag"
            fields[0] = scaffold_name
            agp_lines.append("\t".join(fields))
    return agp_lines


def _parse_ragtag_confidence(
    confidence_path: Path,
) -> Dict[str, Tuple[float, float, float]]:
    """Parse RagTag confidence.txt -> {contig: (grouping, location, orientation)}.

    Args:
        confidence_path: Path to ragtag.scaffold.confidence.txt.

    Returns:
        Dict mapping contig name -> (grouping_conf, location_conf, orientation_conf).
    """
    result: Dict[str, Tuple[float, float, float]] = {}
    if not confidence_path.exists():
        return result
    with confidence_path.open() as fh:
        header = fh.readline()  # skip header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            contig = fields[0]
            try:
                grouping = float(fields[1])
                location = float(fields[2])
                orientation = float(fields[3])
            except (ValueError, IndexError):
                continue
            result[contig] = (grouping, location, orientation)
    return result


def _extract_ragtag_scaffold_seq(
    fasta_path: Path,
) -> Optional[str]:
    """Extract the main scaffold sequence from RagTag output.

    RagTag produces a FASTA with the scaffold and possibly unplaced contigs.
    The scaffold is the first sequence ending in '_RagTag'.

    Returns:
        Scaffold sequence string, or None if not found.
    """
    seqs = read_fasta_sequences(fasta_path)
    # Find the scaffold sequence (ends with _RagTag)
    for name, seq in seqs.items():
        if name.endswith("_RagTag"):
            return seq
    # Fall back to first sequence
    if seqs:
        return next(iter(seqs.values()))
    return None


# ---------------------------------------------------------------------------
# Main scaffolding function
# ---------------------------------------------------------------------------

def scaffold_chromosomes(
    query_fasta: Path,
    classifications: List[ContigClassification],
    contig_orientations: Dict[str, bool],
    qr_ref_ranges: Dict[Tuple[str, str], Tuple[int, int]],
    ref_fasta: Path,
    ref_lengths: Dict[str, int],
    best_ref: Dict[str, str],
    contig_refs: Dict[str, Set[str]],
    work_dir: Path,
    threads: int = 1,
    gap_size: int = 100,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
    qr_best_chain_ident: Optional[Dict[Tuple[str, str], float]] = None,
) -> Tuple[Dict[str, str], List[str], Dict[str, Tuple[float, float, float]]]:
    """Produce scaffolded chromosome pseudomolecules.

    Groups contigs by (reference chromosome, haplotype group), then either:
    - Emits trivial AGP for T2T single-contig groups (even if debris exists)
    - Runs RagTag scaffold per group (if installed)
    - Falls back to built-in ordering by reference position

    Args:
        query_fasta: Path to query assembly FASTA.
        classifications: List of ContigClassification objects.
        contig_orientations: Dict mapping contig -> True if reversed.
        qr_ref_ranges: Dict mapping (contig, ref_id) -> (ref_min, ref_max).
        ref_fasta: Path to reference FASTA.
        ref_lengths: Dict mapping ref_id -> length.
        best_ref: Dict mapping contig -> best ref_id.
        contig_refs: Dict mapping contig -> set of ref_ids with evidence.
        work_dir: Working directory for intermediate files.
        threads: Number of threads for RagTag.
        gap_size: Number of Ns between scaffolded contigs.
        ref_norm_to_orig: Optional mapping from normalized to original ref IDs.
        qr_best_chain_ident: Optional dict of (contig, ref_id) -> identity for
            haplotype-aware debris assignment.

    Returns:
        Tuple of (scaffolded_sequences, agp_lines, scaffold_confidences) where:
        - scaffolded_sequences: Dict mapping scaffold name -> DNA sequence
        - agp_lines: List of AGP lines for all scaffolds
        - scaffold_confidences: Dict mapping contig original_name ->
          (grouping, location, orientation) confidence from RagTag
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Group contigs by (ref_id, haplotype group)
    groups = _group_contigs_by_haplotype(
        classifications, best_ref, contig_refs, qr_best_chain_ident or {},
    )
    all_confidences: Dict[str, Tuple[float, float, float]] = {}

    if not groups:
        logger.warning("No contigs to scaffold")
        return {}, [], {}

    # Build classification lookup
    clf_lookup = {clf.original_name: clf for clf in classifications}

    # Check RagTag availability once
    have_ragtag = have_exe("ragtag.py")
    if not have_ragtag:
        logger.warning(
            "ragtag.py not found in PATH -- using built-in scaffolder. "
            "Install RagTag for better results: conda install -c bioconda ragtag"
        )

    # Load query sequences
    logger.info("Loading query sequences for scaffolding")
    query_seqs = read_fasta_sequences(query_fasta)

    # Load reference sequences (for RagTag)
    ref_seqs: Dict[str, str] = {}
    if have_ragtag:
        ref_seqs = read_fasta_sequences(ref_fasta)

    # Process each (ref_id, grp)
    all_scaffolded: Dict[str, str] = {}
    all_agp: List[str] = []

    # Sort for deterministic output
    for (ref_id, grp) in sorted(groups.keys()):
        contig_names = groups[(ref_id, grp)]
        if not contig_names:
            continue

        # Determine scaffold name
        base_name = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
        if grp <= 1:
            scaffold_name = base_name
        else:
            scaffold_name = f"{base_name}_{chr(ord('A') + grp - 1)}"

        # T2T handling: emit each T2T contig as its own trivial scaffold.
        # When T2T contigs exist, remaining contigs (fragments, haplotigs,
        # debris) are excluded from scaffolded output — the T2T contig(s)
        # already represent the complete chromosome.
        chrom_assigned_in_group = [
            c for c in contig_names
            if clf_lookup.get(c) and clf_lookup[c].classification == "chrom_assigned"
        ]
        t2t_contigs = [
            c for c in chrom_assigned_in_group
            if clf_lookup[c].is_full_length
            and clf_lookup[c].has_5p_telomere
            and clf_lookup[c].has_3p_telomere
        ]

        if t2t_contigs:
            n_other = len(contig_names) - len(t2t_contigs)

            for t2t_contig in t2t_contigs:
                # Use clf.new_name when multiple T2T copies, base scaffold_name when single
                if len(t2t_contigs) == 1:
                    t2t_name = scaffold_name
                else:
                    t2t_name = clf_lookup[t2t_contig].new_name or scaffold_name

                if n_other > 0:
                    logger.info(
                        f"{t2t_name}: T2T contig {t2t_contig}, "
                        f"dropping {n_other} other contig(s) from scaffold"
                    )
                else:
                    logger.info(f"{t2t_name}: T2T contig {t2t_contig}, trivial scaffold")

                seq = query_seqs.get(t2t_contig, "")
                if seq:
                    if contig_orientations.get(t2t_contig, False):
                        seq = reverse_complement(seq)
                        orientation = "-"
                    else:
                        orientation = "+"
                    all_scaffolded[t2t_name] = seq
                    all_agp.append(_agp_component_line(
                        t2t_name, 1, len(seq), 1,
                        t2t_contig, 1, len(seq), orientation,
                    ))

            continue  # Skip scaffolding — T2T contigs cover the chromosome

        # Single-contig group: emit trivial AGP without invoking RagTag
        if len(contig_names) == 1:
            contig = contig_names[0]
            seq = query_seqs.get(contig, "")
            if seq:
                if contig_orientations.get(contig, False):
                    seq = reverse_complement(seq)
                    orientation = "-"
                else:
                    orientation = "+"
                all_scaffolded[scaffold_name] = seq
                all_agp.append(_agp_component_line(
                    scaffold_name, 1, len(seq), 1,
                    contig, 1, len(seq), orientation,
                ))
                logger.info(f"{scaffold_name}: single contig {contig}, trivial scaffold")
            continue

        # Multi-contig group: scaffold
        n_contigs = len(contig_names)
        logger.info(f"{scaffold_name}: scaffolding {n_contigs} contigs")

        if have_ragtag and ref_id in ref_seqs:
            # Extract reference chromosome to temp FASTA
            chr_work = work_dir / scaffold_name
            chr_work.mkdir(parents=True, exist_ok=True)
            ref_chr_fasta = chr_work / "ref_chr.fa"
            # Find original ref ID in ref_seqs
            ref_seq_key = None
            for key in ref_seqs:
                orig = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
                if key == orig or key == ref_id:
                    ref_seq_key = key
                    break
            if ref_seq_key is None:
                # Try to find by iterating ref_norm_to_orig
                if ref_norm_to_orig:
                    orig_id = ref_norm_to_orig.get(ref_id, ref_id)
                    if orig_id in ref_seqs:
                        ref_seq_key = orig_id

            if ref_seq_key and ref_seq_key in ref_seqs:
                write_fasta({ref_seq_key: ref_seqs[ref_seq_key]}, ref_chr_fasta)
                # Remove stale .fai so pysam/samtools re-indexes
                ref_fai = Path(str(ref_chr_fasta) + ".fai")
                if ref_fai.exists():
                    ref_fai.unlink()

                # Extract contigs to temp FASTA
                contigs_fasta = chr_work / "contigs.fa"
                contig_seqs_subset = {}
                for cn in contig_names:
                    if cn in query_seqs:
                        contig_seqs_subset[cn] = query_seqs[cn]
                write_fasta(contig_seqs_subset, contigs_fasta)
                # Remove stale .fai so pysam/samtools re-indexes
                contigs_fai = Path(str(contigs_fasta) + ".fai")
                if contigs_fai.exists():
                    contigs_fai.unlink()

                # Run RagTag
                ragtag_dir = chr_work / "ragtag_out"
                agp_path, fasta_path = _run_ragtag_scaffold(
                    ref_chr_fasta, contigs_fasta, ragtag_dir, threads,
                )

                if agp_path and fasta_path:
                    agp_lines = _parse_ragtag_agp(agp_path, scaffold_name)
                    scaffold_seq = _extract_ragtag_scaffold_seq(fasta_path)
                    if scaffold_seq and agp_lines:
                        all_scaffolded[scaffold_name] = scaffold_seq
                        all_agp.extend(agp_lines)
                        # Parse RagTag confidence scores
                        conf_path = ragtag_dir / "ragtag.scaffold.confidence.txt"
                        conf = _parse_ragtag_confidence(conf_path)
                        all_confidences.update(conf)
                        continue
                    else:
                        logger.warning(f"{scaffold_name}: RagTag produced empty output, falling back")
                else:
                    logger.warning(f"{scaffold_name}: RagTag failed, falling back to built-in scaffolder")
            else:
                logger.warning(f"{scaffold_name}: reference sequence not found, falling back")

        # Fallback: built-in scaffolder
        scaffold_seq, agp_lines = _builtin_scaffold(
            contig_names=contig_names,
            sequences=query_seqs,
            orientations=contig_orientations,
            ref_ranges=qr_ref_ranges,
            ref_id=scaffold_name,
            gap_size=gap_size,
        )
        if scaffold_seq:
            all_scaffolded[scaffold_name] = scaffold_seq
            all_agp.extend(agp_lines)

    logger.done(f"Scaffolded {len(all_scaffolded)} chromosomes")
    return all_scaffolded, all_agp, all_confidences
