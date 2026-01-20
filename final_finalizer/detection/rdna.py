#!/usr/bin/env python3
"""
rDNA detection for final_finalizer.

Contains functions for identifying contigs with significant ribosomal DNA content.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional, Set

from final_finalizer.detection.blast import (
    parse_blast_coverage,
    run_blastn_megablast,
    run_makeblastdb,
)


def prepare_rdna_reference(
    rdna_ref_arg: Optional[str],
    script_dir: Path,
) -> Optional[Path]:
    """Prepare rDNA reference FASTA.

    If rdna_ref_arg is None or 'default', use data/athal-45s-ref.fa.
    Otherwise, use the provided path.
    """
    if rdna_ref_arg is None or rdna_ref_arg.lower() == "default":
        default_path = script_dir / "data" / "athal-45s-ref.fa"
        if default_path.exists():
            print(f"[info] Using default rDNA reference: {default_path}", file=sys.stderr)
            return default_path
        else:
            print(f"[warn] Default rDNA reference not found: {default_path}", file=sys.stderr)
            return None

    rdna_path = Path(rdna_ref_arg)
    if rdna_path.exists():
        print(f"[info] Using user-provided rDNA reference: {rdna_path}", file=sys.stderr)
        return rdna_path

    print(f"[warn] rDNA reference not found: {rdna_ref_arg}", file=sys.stderr)
    return None


def detect_rdna_contigs(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    rdna_ref: Path,
    work_dir: Path,
    threads: int,
    min_coverage: float,
    exclude_contigs: Set[str],
) -> Set[str]:
    """Identify contigs with significant rDNA content.

    Returns set of contig names with coverage >= min_coverage.
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

    # Parse results
    blast_results = parse_blast_coverage(blast_out, query_lengths)

    rdna_contigs: Set[str] = set()
    for contig, summary in blast_results.items():
        if contig in exclude_contigs:
            continue
        if summary.total_coverage >= min_coverage:
            rdna_contigs.add(contig)
            print(f"[info] rDNA contig: {contig} (coverage={summary.total_coverage:.2f})", file=sys.stderr)

    return rdna_contigs
