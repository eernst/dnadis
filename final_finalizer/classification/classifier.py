#!/usr/bin/env python3
"""
Classification functions for final_finalizer.

Contains functions for:
- Contig orientation determination
- Gene counting and statistics
- Reference-based debris classification
- Contig naming
- Classification pipeline
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.alignment.external_tools import get_minimap2_exe, run_minimap2
from final_finalizer.models import (
    ChainEvidenceResult,
    ContaminantHit,
    ContigClassification,
    DebrisHit,
    OrganelleHit,
    RdnaHit,
)
from final_finalizer.utils.io_utils import merge_intervals, open_maybe_gzip
from final_finalizer.utils.reference_utils import normalize_ref_id, split_chrom_subgenome
from final_finalizer.utils.sequence_utils import read_fasta_sequences, write_fasta


# ----------------------------
# Orientation determination
# ----------------------------
def compute_orientation_votes(
    macro_block_rows: List[Tuple],
    contig: str,
    assigned_ref_id: str,
) -> Tuple[int, int]:
    """Count forward/reverse orientation votes from macro blocks.

    Each macro block provides one vote based on its strand.

    Returns (fwd_count, rev_count)
    """
    fwd_count = 0
    rev_count = 0

    for row in macro_block_rows:
        # Row format: (contig, contig_len, ref_id, chrom_id, subgenome, strand, ...)
        if len(row) < 6:
            continue
        row_contig = row[0]
        row_ref_id = row[2]
        row_strand = row[5]

        if row_contig != contig:
            continue
        if row_ref_id != assigned_ref_id:
            continue

        if row_strand == "+":
            fwd_count += 1
        elif row_strand == "-":
            rev_count += 1

    return fwd_count, rev_count


def determine_contig_orientations(
    macro_block_rows: List[Tuple],
    best_ref: Dict[str, str],
    chromosome_contigs: Set[str],
) -> Dict[str, bool]:
    """Determine which chromosome contigs need to be reverse-complemented.

    Only chromosome-assigned contigs are subject to reorientation based on synteny
    block strand votes. Non-chromosome contigs (debris, contaminants, unclassified)
    are left in their original orientation since they lack reliable synteny evidence.

    Returns:
        Dict[str, bool]: Mapping of contig name to orientation flag.
            True = contig should be reverse-complemented to match reference orientation.
            False = contig is already in correct orientation.
            Only chromosome contigs are included in the returned dict.
    """
    orientations: Dict[str, bool] = {}

    for contig in chromosome_contigs:
        assigned_ref = best_ref.get(contig, "")
        if not assigned_ref:
            orientations[contig] = False
            continue

        fwd, rev = compute_orientation_votes(macro_block_rows, contig, assigned_ref)

        # Reverse if more reverse votes than forward
        should_reverse = rev > fwd
        orientations[contig] = should_reverse

        if should_reverse:
            print(f"[info] Contig {contig} will be reverse-complemented (fwd={fwd}, rev={rev})", file=sys.stderr)

    return orientations


# ----------------------------
# Gene count statistics
# ----------------------------
def count_genes_per_ref_chrom(
    gff3_path: Path,
    ref_id_map: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """Count genes per reference chromosome from GFF3.

    Counts unique gene IDs per chromosome.

    Returns dict: ref_id -> gene_count
    """
    chrom_genes: Dict[str, Set[str]] = defaultdict(set)

    with open_maybe_gzip(gff3_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _src, ftype, _start, _end, _score, _strand, _phase, attrs = parts

            if ftype.lower() != "gene":
                continue

            # Parse gene ID from attributes
            gene_id = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("ID="):
                    gene_id = kv[3:]
                    break

            if gene_id:
                ref_id = ref_id_map.get(seqid, normalize_ref_id(seqid)) if ref_id_map else seqid
                chrom_genes[ref_id].add(gene_id)

    return {chrom: len(genes) for chrom, genes in chrom_genes.items()}


def compute_mean_genes_per_Mbp(
    qr_gene_count: Dict[Tuple[str, str], int],
    query_lengths: Dict[str, int],
    chromosome_contigs: Set[str],
    best_ref: Dict[str, str],
) -> float:
    """Compute mean genes per Mbp across chromosome-assigned contigs."""
    total_genes = 0
    total_bp = 0

    for contig in chromosome_contigs:
        ref_id = best_ref.get(contig, "")
        if not ref_id:
            continue
        gene_count = qr_gene_count.get((contig, ref_id), 0)
        contig_len = query_lengths.get(contig, 0)

        total_genes += gene_count
        total_bp += contig_len

    if total_bp == 0:
        return 0.0

    mean_gpmbp = total_genes / (total_bp / 1_000_000.0)
    print(f"[info] Mean genes per Mbp (query): {mean_gpmbp:.2f} ({total_genes} genes / {total_bp} bp)", file=sys.stderr)
    return mean_gpmbp


# ----------------------------
# Reference-based debris classification
# ----------------------------
def classify_debris_and_unclassified(
    remaining_contigs: Set[str],
    query_fasta: Path,
    query_lengths: Dict[str, int],
    ref_fasta: Path,
    chrs_fasta: Optional[Path],
    qr_gene_count: Dict[Tuple[str, str], int],
    work_dir: Path,
    threads: int,
    min_coverage: float,
    min_protein_hits: int,
) -> Tuple[Set[str], Set[str], Dict[str, DebrisHit]]:
    """Classify remaining contigs as debris or unclassified based on reference homology.

    MOTIVATION: After chromosome debris detection (which catches assembly artifacts that
    are near-identical to assembled chromosomes), some contigs may still have homology to
    the reference genome without being classified. These include:
    - Divergent haplotigs that don't align well at nucleotide level to assembled chromosomes
      but retain protein-coding content (detected via miniprot hits)
    - Sequences from genomic regions not represented in the primary chromosome contigs
      but present in the reference (detected via reference nucleotide alignment)

    This reference-based approach complements chromosome debris detection by catching
    sequences with ancestral/evolutionary homology rather than assembly-specific duplicates.

    Classification criteria:
    - Debris: contigs with 50%+ nucleotide alignment coverage to reference OR 2+ miniprot hits
    - Unclassified: everything else (potential novel sequences, contaminants missed by
      earlier screening, or highly divergent sequences)

    Args:
        remaining_contigs: Set of contigs not yet classified
        query_fasta: Path to query FASTA
        query_lengths: Dict of contig name -> length
        ref_fasta: Path to reference FASTA
        chrs_fasta: Optional path to assembled chromosome FASTA (currently unused)
        qr_gene_count: Dict of (contig, ref_id) -> gene count from miniprot
        work_dir: Working directory
        threads: Number of threads
        min_coverage: Min alignment coverage fraction for debris classification
        min_protein_hits: Min miniprot protein hits for debris classification

    Returns:
        Tuple of (debris_contigs, unclassified_contigs, debris_hits)
        where debris_hits maps contig name to DebrisHit with coverage and identity details
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    debris_contigs: Set[str] = set()
    unclassified_contigs: Set[str] = set()
    debris_hits: Dict[str, DebrisHit] = {}

    # Check protein hits for each remaining contig
    contigs_with_protein_hits: Dict[str, int] = {}  # contig -> protein hit count
    for (contig, _), gene_count in qr_gene_count.items():
        if contig in remaining_contigs and gene_count >= min_protein_hits:
            # Track max hits across all ref chroms
            if contig not in contigs_with_protein_hits or gene_count > contigs_with_protein_hits[contig]:
                contigs_with_protein_hits[contig] = gene_count

    # For contigs without enough protein hits, check nucleotide alignment
    contigs_needing_alignment = remaining_contigs - set(contigs_with_protein_hits.keys())

    # Track coverage/identity from alignments
    contig_coverage: Dict[str, float] = {}
    contig_identity: Dict[str, float] = {}

    if contigs_needing_alignment:
        # Write subset of contigs to align
        subset_fasta = work_dir / "debris_candidates.fa"
        query_seqs = read_fasta_sequences(query_fasta)
        subset_seqs = {k: v for k, v in query_seqs.items() if k in contigs_needing_alignment}
        if subset_seqs:
            write_fasta(subset_seqs, subset_fasta)

            # Align against reference using minimap2/mm2plus if available
            if get_minimap2_exe() and subset_fasta.exists():
                paf_out = work_dir / "debris_vs_ref.paf"

                run_minimap2(
                    ref=ref_fasta,
                    qry=subset_fasta,
                    paf_out=paf_out,
                    threads=threads,
                    preset="asm20",
                    extra_args=["-c"],  # output CIGAR for identity calculation
                )

                # Parse PAF to get coverage and identity
                if paf_out.exists():
                    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
                    query_matches: Dict[str, int] = defaultdict(int)
                    query_alnlen: Dict[str, int] = defaultdict(int)

                    with paf_out.open("r") as fh:
                        for line in fh:
                            fields = line.rstrip("\n").split("\t")
                            if len(fields) < 12:
                                continue
                            qname = fields[0]
                            try:
                                qs = int(fields[2])
                                qe = int(fields[3])
                                matches = int(fields[9])
                                aln_len = int(fields[10])
                            except ValueError:
                                continue
                            if qs > qe:
                                qs, qe = qe, qs
                            query_intervals[qname].append((qs, qe))
                            query_matches[qname] += matches
                            query_alnlen[qname] += aln_len

                    for qname, intervals in query_intervals.items():
                        _, total_bp = merge_intervals(intervals)
                        qlen = query_lengths.get(qname, 0)
                        coverage = (total_bp / qlen) if qlen > 0 else 0.0
                        identity = (query_matches[qname] / query_alnlen[qname]) if query_alnlen[qname] > 0 else 0.0

                        contig_coverage[qname] = coverage
                        contig_identity[qname] = identity

                        if coverage >= min_coverage:
                            debris_contigs.add(qname)
                            debris_hits[qname] = DebrisHit(
                                coverage=coverage,
                                identity=identity,
                                protein_hits=0,
                                source="reference",
                            )

    # Add contigs with protein hits to debris
    for contig, protein_hits in contigs_with_protein_hits.items():
        debris_contigs.add(contig)
        debris_hits[contig] = DebrisHit(
            coverage=contig_coverage.get(contig, 0.0),
            identity=contig_identity.get(contig, 0.0),
            protein_hits=protein_hits,
            source="reference",
        )

    # Everything else is unclassified
    unclassified_contigs = remaining_contigs - debris_contigs

    print(f"[info] Debris contigs: {len(debris_contigs)}", file=sys.stderr)
    print(f"[info] Unclassified contigs: {len(unclassified_contigs)}", file=sys.stderr)

    return debris_contigs, unclassified_contigs, debris_hits


# ----------------------------
# Contig naming
# ----------------------------
def generate_contig_names(
    classifications: List[ContigClassification],
    query_lengths: Dict[str, int],
    add_subgenome_suffix: Optional[str],
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """Generate new contig names based on classification.

    - Chromosome contigs: inherit reference name (e.g., chr5T)
    - Multiple contigs to same ref: add _1, _2 suffix by descending length
    - Non-chromosome: contig_N (N increases with decreasing length)
    """
    name_mapping: Dict[str, str] = {}

    # Group chromosome contigs by assigned reference
    ref_to_contigs: Dict[str, List[str]] = defaultdict(list)
    non_chrom_contigs: List[str] = []

    for clf in classifications:
        if clf.classification == "chrom_assigned" and clf.assigned_ref_id:
            ref_to_contigs[clf.assigned_ref_id].append(clf.original_name)
        else:
            non_chrom_contigs.append(clf.original_name)

    # Name chromosome contigs
    for ref_id, contigs in ref_to_contigs.items():
        # Sort by length descending
        contigs.sort(key=lambda c: query_lengths.get(c, 0), reverse=True)

        base_name = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
        # Add subgenome suffix if reference doesn't have one and user requested it
        if add_subgenome_suffix:
            chrom_id, sub = split_chrom_subgenome(base_name)
            if sub == "NA":
                base_name = f"{chrom_id}{add_subgenome_suffix}"

        if len(contigs) == 1:
            name_mapping[contigs[0]] = base_name
        else:
            for i, contig in enumerate(contigs, start=1):
                name_mapping[contig] = f"{base_name}_{i}"

    # Name non-chromosome contigs
    # Sort by length descending
    non_chrom_contigs.sort(key=lambda c: query_lengths.get(c, 0), reverse=True)
    for i, contig in enumerate(non_chrom_contigs, start=1):
        name_mapping[contig] = f"contig_{i}"

    return name_mapping


# ----------------------------
# Classification pipeline
# ----------------------------
def classify_all_contigs(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    best_ref: Dict[str, str],
    chr_like_minlen: int,
    ev: ChainEvidenceResult,
    ref_gene_counts: Dict[str, int],
    chrC_contig: Optional[str],
    chrM_contig: Optional[str],
    organelle_debris: Set[str],
    rdna_contigs: Set[str],
    contaminants: Dict[str, ContaminantHit],
    chromosome_debris: Set[str],
    other_debris: Set[str],
    add_subgenome_suffix: Optional[str],
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
    query_gc: Optional[Dict[str, float]] = None,
    ref_gc_mean: Optional[float] = None,
    ref_gc_std: Optional[float] = None,
    asm_gc_mean: Optional[float] = None,
    asm_gc_std: Optional[float] = None,
    organelle_hits: Optional[Dict[str, OrganelleHit]] = None,
    rdna_hits: Optional[Dict[str, RdnaHit]] = None,
    chrom_debris_hits: Optional[Dict[str, DebrisHit]] = None,
    other_debris_hits: Optional[Dict[str, DebrisHit]] = None,
) -> List[ContigClassification]:
    """Classify all contigs and assign confidence levels.

    This function classifies each contig into one of nine categories and assigns
    a confidence level (high/medium/low) based on multiple evidence sources.

    Classification categories:
    - chrom_assigned: Chromosome-length contigs with synteny support
    - chrom_unassigned: Chromosome-length contigs without synteny
    - organelle_complete: Complete organelle genomes (chrC, chrM)
    - organelle_debris: Partial organelle sequences
    - rDNA: Ribosomal DNA repeat units
    - contaminant: Sequences from contaminating organisms
    - chrom_debris: High-coverage duplicates of assembled chromosomes
    - debris: Assembly fragments with reference homology
    - unclassified: Sequences with no classification evidence

    Confidence levels are based on:
    - Gene proportion: Fraction of reference genes aligned to contig
    - GC deviation: Standard deviations from GC baseline
      - chrom_assigned: compared to reference nuclear GC
      - other categories: compared to assembly chromosome GC (more appropriate for divergent genomes)
    - Coverage: Alignment coverage fraction (0.0-1.0)
    - Identity: Alignment identity (0.0-1.0)
    - Protein hits: Number of miniprot alignments

    Args:
        query_fasta: Path to query FASTA
        query_lengths: Dict of contig name -> length
        best_ref: Dict of contig -> best reference assignment
        chr_like_minlen: Minimum length for chromosome classification
        ev: Chain evidence results
        ref_gene_counts: Dict of ref_id -> gene count
        chrC_contig: Chloroplast contig name
        chrM_contig: Mitochondrial contig name
        organelle_debris: Set of organelle debris contigs
        rdna_contigs: Set of rDNA contigs
        contaminants: Dict of contig -> ContaminantHit with taxid, name, coverage, score
        chromosome_debris: Set of chromosome debris contigs (from chr-vs-chr alignment)
        other_debris: Set of other debris contigs (from ref/protein alignment)
        add_subgenome_suffix: Optional subgenome suffix to add
        query_gc: Dict of contig name -> GC content (0.0-1.0)
        ref_gc_mean: Mean GC content of reference nuclear chromosomes
        ref_gc_std: Standard deviation of reference GC content
        asm_gc_mean: Mean GC content of assembly chromosome contigs
        asm_gc_std: Standard deviation of assembly chromosome GC content
        organelle_hits: Dict of contig -> OrganelleHit with detection details
        rdna_hits: Dict of contig -> RdnaHit with detection details
        chrom_debris_hits: Dict of contig -> DebrisHit from chromosome debris detection
        other_debris_hits: Dict of contig -> DebrisHit from reference/protein detection

    Returns:
        List of ContigClassification objects with assigned names, categories, and confidence levels
    """

    classifications: List[ContigClassification] = []
    classified_contigs: Set[str] = set()

    # Helper to compute GC deviation vs REFERENCE (for chrom_assigned validation)
    def _gc_deviation_vs_ref(contig: str) -> Optional[float]:
        if query_gc is None or ref_gc_mean is None or ref_gc_std is None:
            return None
        gc = query_gc.get(contig)
        if gc is None:
            return None
        if ref_gc_std == 0:
            return 0.0
        return abs(gc - ref_gc_mean) / ref_gc_std

    # Helper to compute GC deviation vs ASSEMBLY chromosomes (for non-chrom classifications)
    # This is more appropriate for divergent genomes where the assembly may have
    # different GC than the reference
    def _gc_deviation_vs_asm(contig: str) -> Optional[float]:
        if query_gc is None or asm_gc_mean is None or asm_gc_std is None:
            return None
        gc = query_gc.get(contig)
        if gc is None:
            return None
        if asm_gc_std == 0:
            return 0.0
        return abs(gc - asm_gc_mean) / asm_gc_std

    # Helper to get contig GC content
    def _get_gc(contig: str) -> Optional[float]:
        return query_gc.get(contig) if query_gc else None

    # 1. Chromosome-length contigs with reference assignment
    for contig, ref_id in best_ref.items():
        if not ref_id:
            continue
        contig_len = query_lengths.get(contig, 0)
        if contig_len < chr_like_minlen:
            continue

        # Check if it passed synteny gates
        gene_count = ev.qr_gene_count.get((contig, ref_id), 0)
        ref_total_genes = ref_gene_counts.get(ref_id, 0)
        gene_proportion = (gene_count / ref_total_genes) if ref_total_genes > 0 else None
        gc_dev = _gc_deviation_vs_ref(contig)  # Use reference GC for chrom_assigned validation

        # Compute synteny score (0-1) based on gene proportion
        synteny_score = min(1.0, gene_proportion * 2) if gene_proportion else 0.0

        # Compute confidence
        # High: good synteny (>20% genes) AND GC similar to reference (<2 std)
        # Medium: moderate synteny OR borderline GC
        # Low: weak synteny AND/OR very different GC
        confidence = "high"
        if gene_proportion is None or gene_proportion < 0.1:
            confidence = "low"
        elif gene_proportion < 0.2:
            confidence = "medium"
        # GC deviation penalty
        if gc_dev is not None and gc_dev > 3.0:
            confidence = "low"
        elif gc_dev is not None and gc_dev > 2.0 and confidence == "high":
            confidence = "medium"

        classifications.append(ContigClassification(
            original_name=contig,
            new_name="",  # Will be filled later
            classification="chrom_assigned",
            reversed=False,  # Will be filled later
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id=ref_id,
            ref_gene_proportion=gene_proportion,
            contig_len=contig_len,
            gc_content=_get_gc(contig),
            gc_deviation=gc_dev,
            synteny_score=synteny_score,
            classification_confidence=confidence,
        ))
        classified_contigs.add(contig)

    # 1b. Chromosome-length contigs WITHOUT reference assignment
    # (but not contaminants - those are handled separately)
    for contig, contig_len in query_lengths.items():
        if contig in classified_contigs:
            continue
        if contig_len < chr_like_minlen:
            continue
        # Skip known contaminants - they should be classified as contaminants, not chrom_unassigned
        if contig in contaminants:
            continue

        gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline

        # Confidence is inherently low/medium for unassigned contigs
        # Medium: GC similar to assembly chromosomes (might be a real chromosome from a divergent region)
        # Low: GC different from assembly (might be contaminant or unusual sequence)
        confidence = "medium"
        if gc_dev is not None and gc_dev > 2.0:
            confidence = "low"

        # This contig is chromosome-length but has no synteny assignment
        classifications.append(ContigClassification(
            original_name=contig,
            new_name="",  # Will be filled later
            classification="chrom_unassigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id=None,
            ref_gene_proportion=None,
            contig_len=contig_len,
            gc_content=_get_gc(contig),
            gc_deviation=gc_dev,
            synteny_score=0.0,  # No synteny evidence
            classification_confidence=confidence,
        ))
        classified_contigs.add(contig)

    # 2. Organelle contigs
    # organelle_complete: High confidence - passed strict detection criteria
    # (coverage >80%, length within tolerance)
    if chrC_contig and chrC_contig not in classified_contigs:
        hit = organelle_hits.get(chrC_contig) if organelle_hits else None
        # High confidence if coverage >= 90%, medium if >= 80% (our threshold)
        confidence = "high"
        if hit and hit.coverage < 0.90:
            confidence = "medium"
        classifications.append(ContigClassification(
            original_name=chrC_contig,
            new_name="chrC",
            classification="organelle_complete",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chrC",
            ref_gene_proportion=None,
            contig_len=query_lengths.get(chrC_contig, 0),
            gc_content=_get_gc(chrC_contig),
            gc_deviation=None,  # Organelles have different GC, don't compare to nuclear
            classification_confidence=confidence,
        ))
        classified_contigs.add(chrC_contig)

    if chrM_contig and chrM_contig not in classified_contigs:
        hit = organelle_hits.get(chrM_contig) if organelle_hits else None
        # High confidence if coverage >= 90%, medium if >= 80% (our threshold)
        confidence = "high"
        if hit and hit.coverage < 0.90:
            confidence = "medium"
        classifications.append(ContigClassification(
            original_name=chrM_contig,
            new_name="chrM",
            classification="organelle_complete",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chrM",
            ref_gene_proportion=None,
            contig_len=query_lengths.get(chrM_contig, 0),
            gc_content=_get_gc(chrM_contig),
            gc_deviation=None,  # Organelles have different GC, don't compare to nuclear
            classification_confidence=confidence,
        ))
        classified_contigs.add(chrM_contig)

    # Organelle debris: Medium confidence - partial organelle match
    for contig in organelle_debris:
        if contig not in classified_contigs:
            hit = organelle_hits.get(contig) if organelle_hits else None
            # Medium confidence by default (partial match)
            # Low confidence if coverage is very low (<60%)
            confidence = "medium"
            if hit and hit.coverage < 0.60:
                confidence = "low"
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="organelle_debris",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=hit.organelle_type if hit else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=None,  # Organelles have different GC
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 3. rDNA contigs
    # rDNA: Confidence based on coverage and identity
    for contig in rdna_contigs:
        if contig not in classified_contigs:
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            hit = rdna_hits.get(contig) if rdna_hits else None
            # High confidence if coverage >= 80% and identity >= 95%
            # Medium confidence if coverage >= 50% (our threshold)
            # Low confidence if borderline
            confidence = "medium"
            if hit:
                if hit.coverage >= 0.80 and hit.identity >= 0.95:
                    confidence = "high"
                elif hit.coverage < 0.60:
                    confidence = "low"
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="rDNA",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 4. Contaminants
    for contig, hit in contaminants.items():
        if contig not in classified_contigs:
            contig_len = query_lengths.get(contig, 0)
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline

            # Compute confidence based on coverage and GC deviation
            # High confidence: high coverage (>0.8) AND (GC deviation > 2 std from assembly)
            # Medium confidence: moderate coverage (0.5-0.8) OR conflicting evidence
            # Low confidence: low coverage (<0.5) AND some synteny evidence
            confidence = "high"
            if hit.coverage < 0.5:
                confidence = "low"
            elif hit.coverage < 0.8:
                confidence = "medium"
            # GC deviation supports contaminant classification
            if gc_dev is not None and gc_dev > 2.0:
                confidence = "high"

            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="contaminant",
                reversed=False,
                contaminant_taxid=hit.taxid,
                contaminant_sci=hit.sci_name,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=contig_len,
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                contam_score=min(1.0, hit.score / 1e9) if hit.score else None,  # Normalize to 0-1 range
                contam_coverage=hit.coverage,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 5. Chromosome debris (from chr-vs-chr alignment detection)
    # Confidence based on coverage, identity, and GC deviation
    for contig in chromosome_debris:
        if contig not in classified_contigs:
            ref_id = best_ref.get(contig, "")
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            hit = chrom_debris_hits.get(contig) if chrom_debris_hits else None
            # High confidence: ≥90% coverage AND ≥95% identity
            # Medium confidence: passed strict 80% cov, 90% identity threshold
            # Low confidence: borderline or GC very different from assembly chromosomes
            confidence = "high"
            if hit:
                if hit.coverage >= 0.90 and hit.identity >= 0.95:
                    confidence = "high"
                elif hit.coverage < 0.85 or hit.identity < 0.92:
                    confidence = "medium"
            if gc_dev is not None and gc_dev > 3.0:
                confidence = "medium" if confidence == "high" else "low"
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="chrom_debris",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=ref_id if ref_id else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 6. Other debris (from reference/protein alignment detection)
    # Confidence based on coverage, identity, protein hits, and GC deviation
    for contig in other_debris:
        if contig not in classified_contigs:
            # Check if this contig has synteny support
            ref_id = best_ref.get(contig, "")
            classification = "chrom_debris" if ref_id else "debris"
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            hit = other_debris_hits.get(contig) if other_debris_hits else None

            # Confidence based on detection evidence
            # High: ≥80% coverage OR ≥5 protein hits
            # Medium: ≥50% coverage OR ≥2 protein hits (our threshold)
            # Low: borderline OR GC very different from assembly chromosomes with no synteny
            confidence = "medium"
            if hit:
                if hit.coverage >= 0.80 or hit.protein_hits >= 5:
                    confidence = "high"
                elif hit.coverage < 0.55 and hit.protein_hits < 3:
                    confidence = "low"

            # GC deviation penalty
            if gc_dev is not None and gc_dev > 3.0 and not ref_id:
                confidence = "low"

            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification=classification,
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=ref_id if ref_id else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 7. Unclassified (everything else)
    # Low confidence by definition - no evidence for any classification
    for contig in query_lengths.keys():
        if contig not in classified_contigs:
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="unclassified",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence="low",
            ))

    # Generate names
    name_mapping = generate_contig_names(
        classifications,
        query_lengths,
        add_subgenome_suffix,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    # Update classifications with names
    for clf in classifications:
        if not clf.new_name:
            clf.new_name = name_mapping.get(clf.original_name, clf.original_name)

    return classifications
