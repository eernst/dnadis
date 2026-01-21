#!/usr/bin/env python3
"""
Organelle (chrC/chrM) detection for final_finalizer.

Contains functions for identifying chloroplast and mitochondrial contigs
using BLAST-based alignment to reference organelle sequences.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.detection.blast import run_makeblastdb, run_blastn_megablast
from final_finalizer.models import OrganelleHit
from final_finalizer.utils.io_utils import merge_intervals
from final_finalizer.utils.reference_utils import normalize_organelle_id
from final_finalizer.utils.sequence_utils import (
    is_hifiasm_circular,
    read_fasta_lengths,
    read_fasta_sequences,
    write_fasta,
)


def extract_organelle_from_ref(
    ref_fasta: Path,
    organelle_id: str,
    output_path: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> bool:
    """Extract a specific sequence (chrC or chrM) from reference FASTA.

    Returns True if found and extracted, False otherwise.
    """
    sequences = read_fasta_sequences(ref_fasta)
    if ref_norm_to_orig:
        candidate = ref_norm_to_orig.get(organelle_id)
        if candidate and candidate in sequences:
            write_fasta({candidate: sequences[candidate]}, output_path)
            return True

    for name, seq in sequences.items():
        if normalize_organelle_id(name) == organelle_id:
            write_fasta({name: seq}, output_path)
            return True
    return False


def prepare_organelle_references(
    ref_fasta: Path,
    chrC_ref_arg: Optional[str],
    chrM_ref_arg: Optional[str],
    work_dir: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Prepare organelle reference FASTAs.

    Priority:
    1. User-provided path (if it exists as a file)
    2. chrC/chrM from reference FASTA if present
    3. None (skip this organelle)

    Returns (chrC_fasta_path, chrM_fasta_path)
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    chrC_path: Optional[Path] = None
    chrM_path: Optional[Path] = None

    # Handle chloroplast reference
    if chrC_ref_arg:
        chrC_arg_path = Path(chrC_ref_arg)
        if chrC_arg_path.exists():
            chrC_path = chrC_arg_path
            print(f"[info] Using user-provided chrC reference: {chrC_path}", file=sys.stderr)
        else:
            print(f"[warn] chrC reference not found: {chrC_ref_arg}", file=sys.stderr)
    else:
        # Try to extract from reference
        chrC_extracted = work_dir / "chrC_ref.fa"
        if extract_organelle_from_ref(ref_fasta, "chrC", chrC_extracted, ref_norm_to_orig):
            chrC_path = chrC_extracted
            print(f"[info] Extracted chrC from reference: {chrC_path}", file=sys.stderr)
        else:
            print("[info] No chrC found in reference FASTA", file=sys.stderr)

    # Handle mitochondrial reference
    if chrM_ref_arg:
        chrM_arg_path = Path(chrM_ref_arg)
        if chrM_arg_path.exists():
            chrM_path = chrM_arg_path
            print(f"[info] Using user-provided chrM reference: {chrM_path}", file=sys.stderr)
        else:
            print(f"[warn] chrM reference not found: {chrM_ref_arg}", file=sys.stderr)
    else:
        # Try to extract from reference
        chrM_extracted = work_dir / "chrM_ref.fa"
        if extract_organelle_from_ref(ref_fasta, "chrM", chrM_extracted, ref_norm_to_orig):
            chrM_path = chrM_extracted
            print(f"[info] Extracted chrM from reference: {chrM_path}", file=sys.stderr)
        else:
            print("[info] No chrM found in reference FASTA", file=sys.stderr)

    return chrC_path, chrM_path


def detect_organelles(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    chrC_ref: Optional[Path],
    chrM_ref: Optional[Path],
    work_dir: Path,
    threads: int,
    min_coverage: float,
    chrC_len_tol: float,
    chrM_len_tol: float,
) -> Tuple[Optional[str], Optional[str], Set[str], Dict[str, OrganelleHit]]:
    """Identify best organelle candidates.

    Criteria for full organelle:
    - BLAST coverage >= min_coverage
    - Length within tolerance of reference
    - Prefer hifiasm circular contigs (ptg*c pattern)

    Returns (chrC_contig, chrM_contig, debris_contigs, organelle_hits)
    where:
    - debris_contigs are those with 50%+ coverage but not selected as main organelle
    - organelle_hits is a dict mapping contig name to OrganelleHit with detection details
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    chrC_contig: Optional[str] = None
    chrM_contig: Optional[str] = None
    debris_contigs: Set[str] = set()
    organelle_hits: Dict[str, OrganelleHit] = {}

    # Get reference lengths
    chrC_len = 0
    chrM_len = 0
    if chrC_ref:
        chrC_lengths = read_fasta_lengths(chrC_ref)
        chrC_len = list(chrC_lengths.values())[0] if chrC_lengths else 0
    if chrM_ref:
        chrM_lengths = read_fasta_lengths(chrM_ref)
        chrM_len = list(chrM_lengths.values())[0] if chrM_lengths else 0

    # Combine organelle references for BLAST if available
    if not chrC_ref and not chrM_ref:
        print("[info] No organelle references available, skipping organelle detection", file=sys.stderr)
        return None, None, set(), {}

    # Create combined organelle reference for BLAST
    # Important: rename sequences to "chrC" and "chrM" so BLAST subject IDs match expected names
    combined_ref = work_dir / "organelle_refs.fa"
    combined_seqs: Dict[str, str] = {}
    if chrC_ref:
        chrC_seqs = read_fasta_sequences(chrC_ref)
        if chrC_seqs:
            # Take the first (or only) sequence and rename to "chrC"
            first_seq = next(iter(chrC_seqs.values()))
            combined_seqs["chrC"] = first_seq
    if chrM_ref:
        chrM_seqs = read_fasta_sequences(chrM_ref)
        if chrM_seqs:
            # Take the first (or only) sequence and rename to "chrM"
            first_seq = next(iter(chrM_seqs.values()))
            combined_seqs["chrM"] = first_seq
    write_fasta(combined_seqs, combined_ref)

    # Create BLAST database
    organelle_db = work_dir / "organelle_refs"
    run_makeblastdb(combined_ref, organelle_db, err_path=work_dir / "makeblastdb.err")

    # Run BLAST
    # Note: organelle genomes (especially chloroplasts with inverted repeats) can produce
    # many HSPs per query-subject pair. Use higher value to capture all alignments.
    blast_out = work_dir / "organelle_blast.txt"
    run_blastn_megablast(
        query_fasta=query_fasta,
        db_paths=[str(organelle_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=1000,
        err_path=work_dir / "blastn.err",
    )

    # Parse BLAST results and compute coverage per (query, subject) pair
    # Also track alignment identity (matches / alignment_length)
    query_subject_intervals: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    query_subject_matches: Dict[Tuple[str, str], int] = defaultdict(int)
    query_subject_alnlen: Dict[Tuple[str, str], int] = defaultdict(int)

    if blast_out.exists() and blast_out.stat().st_size > 0:
        with blast_out.open("r") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue

                qseqid = fields[0]
                sseqid = fields[1]
                try:
                    pident = float(fields[2])
                    aln_length = int(fields[3])
                    qstart = int(fields[6])
                    qend = int(fields[7])
                except ValueError:
                    continue

                if qstart > qend:
                    qstart, qend = qend, qstart

                key = (qseqid, sseqid)
                query_subject_intervals[key].append((qstart, qend))
                # Compute matches from percent identity and alignment length
                matches = int(pident * aln_length / 100.0)
                query_subject_matches[key] += matches
                query_subject_alnlen[key] += aln_length

    # Compute coverage and identify candidates
    # Extended tuple: (contig, coverage, identity, length, is_circular, length_ratio)
    chrC_candidates: List[Tuple[str, float, float, int, bool, float]] = []
    chrM_candidates: List[Tuple[str, float, float, int, bool, float]] = []

    for (qseqid, sseqid), intervals in query_subject_intervals.items():
        _, total_bp = merge_intervals(intervals)
        qlen = query_lengths.get(qseqid, 0)
        coverage = (total_bp / qlen) if qlen > 0 else 0.0

        # Compute identity
        key = (qseqid, sseqid)
        total_matches = query_subject_matches[key]
        total_alnlen = query_subject_alnlen[key]
        identity = (total_matches / total_alnlen) if total_alnlen > 0 else 0.0

        is_circular = is_hifiasm_circular(qseqid)

        if sseqid == "chrC" and chrC_len > 0:
            len_ratio = qlen / chrC_len if chrC_len > 0 else 0.0
            len_diff = abs(qlen - chrC_len) / chrC_len
            if coverage >= min_coverage and len_diff <= chrC_len_tol:
                chrC_candidates.append((qseqid, coverage, identity, qlen, is_circular, len_ratio))
            elif coverage >= 0.50:
                debris_contigs.add(qseqid)
                # Store hit info for debris
                organelle_hits[qseqid] = OrganelleHit(
                    organelle_type="chrC",
                    coverage=coverage,
                    identity=identity,
                    length_ratio=len_ratio,
                    is_complete=False,
                )

        elif sseqid == "chrM" and chrM_len > 0:
            len_ratio = qlen / chrM_len if chrM_len > 0 else 0.0
            len_diff = abs(qlen - chrM_len) / chrM_len
            if coverage >= min_coverage and len_diff <= chrM_len_tol:
                chrM_candidates.append((qseqid, coverage, identity, qlen, is_circular, len_ratio))
            elif coverage >= 0.50:
                debris_contigs.add(qseqid)
                # Store hit info for debris
                organelle_hits[qseqid] = OrganelleHit(
                    organelle_type="chrM",
                    coverage=coverage,
                    identity=identity,
                    length_ratio=len_ratio,
                    is_complete=False,
                )

    # Select best candidates (prefer circular, then by length similarity)
    # Returns (contig_name, coverage, identity, length_ratio) or None
    def select_best(
        candidates: List[Tuple[str, float, float, int, bool, float]], ref_len: int
    ) -> Optional[Tuple[str, float, float, float]]:
        if not candidates:
            return None
        # Sort by: circular (desc), then length difference (asc)
        # Tuple: (contig, coverage, identity, length, is_circular, length_ratio)
        candidates.sort(key=lambda x: (-int(x[4]), abs(x[3] - ref_len)))
        best = candidates[0]
        return (best[0], best[1], best[2], best[5])  # (name, coverage, identity, length_ratio)

    chrC_result = select_best(chrC_candidates, chrC_len)
    chrM_result = select_best(chrM_candidates, chrM_len)

    # Extract contig names and create OrganelleHit entries for complete organelles
    if chrC_result:
        chrC_contig, chrC_cov, chrC_ident, chrC_len_ratio = chrC_result
        debris_contigs.discard(chrC_contig)
        organelle_hits[chrC_contig] = OrganelleHit(
            organelle_type="chrC",
            coverage=chrC_cov,
            identity=chrC_ident,
            length_ratio=chrC_len_ratio,
            is_complete=True,
        )
        print(f"[info] Selected chrC candidate: {chrC_contig} (cov={chrC_cov:.2f}, ident={chrC_ident:.3f})", file=sys.stderr)

    if chrM_result:
        chrM_contig, chrM_cov, chrM_ident, chrM_len_ratio = chrM_result
        debris_contigs.discard(chrM_contig)
        organelle_hits[chrM_contig] = OrganelleHit(
            organelle_type="chrM",
            coverage=chrM_cov,
            identity=chrM_ident,
            length_ratio=chrM_len_ratio,
            is_complete=True,
        )
        print(f"[info] Selected chrM candidate: {chrM_contig} (cov={chrM_cov:.2f}, ident={chrM_ident:.3f})", file=sys.stderr)

    return chrC_contig, chrM_contig, debris_contigs, organelle_hits
