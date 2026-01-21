#!/usr/bin/env python3
"""
Chromosome debris detection for final_finalizer.

Contains functions for identifying contigs that are assembly artifacts -
near-identical copies of chromosome contigs that should be classified as debris.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.alignment.external_tools import get_minimap2_exe, run_minimap2
from final_finalizer.models import DebrisHit
from final_finalizer.utils.io_utils import merge_intervals
from final_finalizer.utils.sequence_utils import write_filtered_fasta


def detect_chromosome_debris(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    chromosome_contigs: Set[str],
    work_dir: Path,
    threads: int,
    min_coverage: float = 0.80,
    min_identity: float = 0.90,
    exclude_contigs: Optional[Set[str]] = None,
) -> Tuple[Set[str], Dict[str, DebrisHit]]:
    """Identify contigs that are assembly artifacts: near-identical copies of chromosome contigs.

    MOTIVATION: Genome assemblers (especially hifiasm) can produce multiple representations
    of the same genomic region - haplotigs, bubble branches, or other duplicates. These
    contigs are nearly identical to portions of the primary chromosome contigs and should
    be classified as debris rather than left as "unclassified" or mistakenly flagged as
    contaminants. This detection method aligns candidate contigs against the *assembled*
    chromosome contigs (not the reference) to catch these assembly-specific artifacts.

    This complements reference-based debris detection (classify_debris_and_unclassified),
    which catches sequences with homology to the reference genome but may miss:
    - Non-coding duplicates (no protein hits)
    - Assembly artifacts specific to this assembly (not in reference)

    Uses minimap2/mm2plus with asm5 preset (assembly-to-assembly, high identity) to find
    contigs with high coverage (>=80%) and high identity (>=90%) to chromosome contigs.

    Args:
        query_fasta: Path to query FASTA containing all contigs
        query_lengths: Dict of contig name -> length
        chromosome_contigs: Set of contigs already classified as chromosomes
        work_dir: Working directory for intermediate files
        threads: Number of threads for alignment
        min_coverage: Min fraction of query aligned to classify as debris (default 0.80)
        min_identity: Min alignment identity to classify as debris (default 0.90)
        exclude_contigs: Set of contigs to exclude from debris detection

    Returns:
        Tuple of:
        - Set of contig names classified as chromosome_debris
        - Dict mapping contig name to DebrisHit with coverage, identity, and source
    """
    if exclude_contigs is None:
        exclude_contigs = set()

    work_dir.mkdir(parents=True, exist_ok=True)

    if not chromosome_contigs:
        print("[info] No chromosome contigs for debris detection", file=sys.stderr)
        return set(), {}

    # Check for minimap2/mm2plus
    if not get_minimap2_exe():
        print("[warn] Neither minimap2 nor mm2plus found, skipping chromosome debris detection", file=sys.stderr)
        return set(), {}

    # Create FASTA of chromosome contigs as reference (streaming to avoid loading entire genome)
    chrs_fasta = work_dir / "chromosome_contigs.fa"
    if not chrs_fasta.exists():
        write_filtered_fasta(query_fasta, chrs_fasta, chromosome_contigs)
    print(f"[info] Created chromosome reference for debris detection: {len(chromosome_contigs)} contigs", file=sys.stderr)

    # Create FASTA of candidate contigs to test (streaming)
    excluded = chromosome_contigs | exclude_contigs
    candidate_contigs = set(query_lengths.keys()) - excluded

    if not candidate_contigs:
        print("[info] No candidate contigs for debris detection", file=sys.stderr)
        return set(), {}

    candidates_fasta = work_dir / "debris_candidates.fa"
    if not candidates_fasta.exists():
        write_filtered_fasta(query_fasta, candidates_fasta, candidate_contigs)

    # Run minimap2/mm2plus: candidates (query) vs chromosomes (reference)
    paf_out = work_dir / "debris_vs_chrs.paf"
    success = run_minimap2(
        ref=chrs_fasta,
        qry=candidates_fasta,
        paf_out=paf_out,
        threads=threads,
        preset="asm5",  # assembly-to-assembly preset, high identity
        extra_args=["-c"],  # output CIGAR
    )
    if not success:
        return set(), {}

    # Parse PAF and compute coverage + identity per query
    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    query_matches: Dict[str, int] = defaultdict(int)
    query_alnlen: Dict[str, int] = defaultdict(int)

    with paf_out.open("r") as fh:
        for line in fh:
            if not line.strip():
                continue
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

    # Identify debris: high coverage and identity
    debris_contigs: Set[str] = set()
    debris_hits: Dict[str, DebrisHit] = {}

    for qname, intervals in query_intervals.items():
        qlen = query_lengths.get(qname, 0)
        if qlen == 0:
            continue

        _, total_bp = merge_intervals(intervals)
        coverage = total_bp / qlen

        total_matches = query_matches[qname]
        total_alnlen = query_alnlen[qname]
        identity = (total_matches / total_alnlen) if total_alnlen > 0 else 0.0

        if coverage >= min_coverage and identity >= min_identity:
            debris_contigs.add(qname)
            debris_hits[qname] = DebrisHit(
                coverage=coverage,
                identity=identity,
                protein_hits=0,  # No protein hits from assembly-to-assembly alignment
                source="chromosome",
            )
            print(
                f"[info] Chromosome debris: {qname} ({qlen:,} bp, cov={coverage:.2f}, ident={identity:.2f})",
                file=sys.stderr,
            )

    print(f"[info] Chromosome debris contigs: {len(debris_contigs)}", file=sys.stderr)
    return debris_contigs, debris_hits
