#!/usr/bin/env python3
"""
rDNA detection for final_finalizer.

Contains functions for identifying contigs with significant ribosomal DNA content.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.detection.blast import (
    run_blastn_megablast,
    run_makeblastdb,
)
from final_finalizer.models import RdnaHit
from final_finalizer.utils.io_utils import merge_intervals
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("rdna")


def prepare_rdna_reference(
    rdna_ref_arg: Optional[str],
    script_dir: Path,
) -> Optional[Path]:
    """Prepare rDNA reference FASTA.

    If rdna_ref_arg is None or 'default', search for data/athal-45s-ref.fa in:
      1. script_dir/data/ (for package installations)
      2. script_dir/../data/ (for running from repo root via shim)
    Otherwise, use the provided path.
    """
    if rdna_ref_arg is None or rdna_ref_arg.lower() == "default":
        # Search multiple locations for the default reference
        search_paths = [
            script_dir / "data" / "athal-45s-ref.fa",
            script_dir.parent / "data" / "athal-45s-ref.fa",
        ]
        for default_path in search_paths:
            if default_path.exists():
                logger.info(f"Using default rDNA reference: {default_path}")
                return default_path

        logger.warning(f"Default rDNA reference not found in: {[str(p) for p in search_paths]}")
        return None

    rdna_path = Path(rdna_ref_arg)
    if rdna_path.exists():
        logger.info(f"Using user-provided rDNA reference: {rdna_path}")
        return rdna_path

    logger.warning(f"rDNA reference not found: {rdna_ref_arg}")
    return None


def detect_rdna_contigs(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    rdna_ref: Path,
    work_dir: Path,
    threads: int,
    min_coverage: float,
    exclude_contigs: Set[str],
) -> Tuple[Set[str], Dict[str, RdnaHit]]:
    """Identify contigs with significant rDNA content.

    Returns:
        Tuple of:
        - set of contig names with coverage >= min_coverage
        - dict mapping contig name to RdnaHit with coverage and identity details
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Create BLAST database for rDNA
    rdna_db = work_dir / "rdna_ref"
    run_makeblastdb(rdna_ref, rdna_db, err_path=work_dir / "makeblastdb_rdna.err")

    # Run BLAST
    # rDNA arrays are highly repetitive; use higher max_hsps to capture full coverage
    blast_out = work_dir / "rdna_blast.txt"
    run_blastn_megablast(
        query_fasta=query_fasta,
        db_paths=[str(rdna_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=100,
        err_path=work_dir / "blastn_rdna.err",
    )

    # Parse results with coverage and identity tracking
    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    query_matches: Dict[str, int] = defaultdict(int)
    query_alnlen: Dict[str, int] = defaultdict(int)

    if blast_out.exists() and blast_out.stat().st_size > 0:
        with blast_out.open("r") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue

                qseqid = fields[0]
                try:
                    pident = float(fields[2])
                    aln_length = int(fields[3])
                    qstart = int(fields[6])
                    qend = int(fields[7])
                except ValueError:
                    continue

                if qstart > qend:
                    qstart, qend = qend, qstart

                query_intervals[qseqid].append((qstart, qend))
                # Compute matches from percent identity and alignment length
                matches = int(pident * aln_length / 100.0)
                query_matches[qseqid] += matches
                query_alnlen[qseqid] += aln_length

    rdna_contigs: Set[str] = set()
    rdna_hits: Dict[str, RdnaHit] = {}

    for qseqid, intervals in query_intervals.items():
        if qseqid in exclude_contigs:
            continue

        _, total_bp = merge_intervals(intervals)
        qlen = query_lengths.get(qseqid, 0)
        coverage = (total_bp / qlen) if qlen > 0 else 0.0

        # Compute identity
        total_matches = query_matches[qseqid]
        total_alnlen = query_alnlen[qseqid]
        identity = (total_matches / total_alnlen) if total_alnlen > 0 else 0.0

        if coverage >= min_coverage:
            rdna_contigs.add(qseqid)
            rdna_hits[qseqid] = RdnaHit(coverage=coverage, identity=identity)
            logger.info(f"rDNA contig: {qseqid} ({qlen:,} bp, cov={coverage:.2f}, ident={identity:.3f})")

    return rdna_contigs, rdna_hits
