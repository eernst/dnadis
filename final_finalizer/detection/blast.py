#!/usr/bin/env python3
"""
BLAST infrastructure for final_finalizer.

Contains functions for running BLAST searches and parsing results.
"""
from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from final_finalizer.models import BlastHitSummary
from final_finalizer.utils.io_utils import have_exe, merge_intervals
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("blast")


def run_makeblastdb(
    fasta_path: Path,
    db_path: Path,
    dbtype: str = "nucl",
    err_path: Optional[Path] = None,
) -> None:
    """Create BLAST database from FASTA file."""
    if not have_exe("makeblastdb"):
        raise RuntimeError("makeblastdb executable not found in PATH")

    # Check if database already exists and is newer than input FASTA
    db_files = list(db_path.parent.glob(f"{db_path.name}.n*"))
    if db_files:
        input_mtime = fasta_path.stat().st_mtime
        oldest_db = min(f.stat().st_mtime for f in db_files)
        if oldest_db >= input_mtime:
            logger.info(f"BLAST database exists, reusing: {db_path}")
            return
        # Input is newer than DB; remove stale DB files and rebuild
        logger.info(f"Input FASTA newer than BLAST database; rebuilding: {db_path}")
        for f in db_files:
            f.unlink()

    db_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-out", str(db_path),
        "-dbtype", dbtype,
        "-hash_index",
    ]

    logger.info(f"Creating BLAST database: {db_path}")

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)
        with err_path.open("wb") as err_fh:
            ret = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=err_fh)
    else:
        ret = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if ret != 0:
        raise RuntimeError(f"makeblastdb failed with return code {ret}")


def run_blastn_megablast(
    query_fasta: Path,
    db_paths: List[str],
    output_path: Path,
    threads: int,
    max_hsps: int = 500,
    max_target_seqs: int = 5,
    taxidlist: Optional[Path] = None,
    err_path: Optional[Path] = None,
) -> None:
    """Run blastn with megablast task.

    Output format: "6 std staxids" (standard 12 columns + taxonomy IDs)
    Supports gzipped query FASTA files by piping decompressed data via stdin.

    Args:
        query_fasta: Query FASTA file (may be gzipped)
        db_paths: List of BLAST database paths
        output_path: Output file path
        threads: Number of threads
        max_hsps: Maximum HSPs per subject
        max_target_seqs: Maximum target sequences per query
        taxidlist: Optional path to file containing taxids to restrict search
        err_path: Optional path for stderr output
    """
    if output_path.exists():
        # Check if query or any DB is newer than output
        out_mtime = output_path.stat().st_mtime
        input_mtime = query_fasta.stat().st_mtime if query_fasta.exists() else 0
        for dbp in db_paths:
            dbp_path = Path(dbp)
            for f in dbp_path.parent.glob(f"{dbp_path.name}.n*"):
                input_mtime = max(input_mtime, f.stat().st_mtime)
        if out_mtime >= input_mtime:
            logger.info(f"BLAST output exists, reusing: {output_path}")
            return
        logger.info(f"Inputs newer than BLAST output; re-running: {output_path}")

    if not have_exe("blastn"):
        raise RuntimeError("blastn executable not found in PATH")

    # Join database paths with spaces
    db_str = " ".join(db_paths)

    # Check if query is gzipped
    is_gzipped = query_fasta.suffix in (".gz", ".bgz")

    # Build common command options
    base_cmd = [
        "blastn",
        "-num_threads", str(threads),
        "-task", "megablast",
        "-db", db_str,
        "-max_hsps", str(max_hsps),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", "6 std staxids",
        "-out", str(output_path),
    ]

    # Add taxid restriction if provided
    if taxidlist and taxidlist.exists():
        base_cmd += ["-taxidlist", str(taxidlist)]

    if is_gzipped:
        # Use stdin for gzipped files
        cmd = base_cmd + ["-query", "-"]
    else:
        cmd = base_cmd + ["-query", str(query_fasta)]

    logger.info(f"Running blastn megablast -> {output_path}")

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)

    if is_gzipped:
        # Use gunzip -c to decompress and pipe to BLAST
        gunzip_cmd = ["gunzip", "-c", str(query_fasta)]
        gunzip_proc = subprocess.Popen(gunzip_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        if err_path:
            with err_path.open("wb") as err_fh:
                ret = subprocess.call(cmd, stdin=gunzip_proc.stdout, stderr=err_fh)
        else:
            ret = subprocess.call(cmd, stdin=gunzip_proc.stdout, stderr=subprocess.DEVNULL)
        gunzip_proc.stdout.close()
        gunzip_proc.wait()
    else:
        if err_path:
            with err_path.open("wb") as err_fh:
                ret = subprocess.call(cmd, stderr=err_fh)
        else:
            ret = subprocess.call(cmd, stderr=subprocess.DEVNULL)

    if ret != 0:
        raise RuntimeError(f"blastn failed with return code {ret}")


def parse_blast_coverage(
    blast_path: Path,
    query_lengths: Dict[str, int],
) -> Dict[str, BlastHitSummary]:
    """Parse BLAST tabular output and compute per-query coverage.

    Expected format: "6 std staxids" (qseqid, sseqid, pident, length, mismatch,
    gapopen, qstart, qend, sstart, send, evalue, bitscore, staxids)

    Returns dict: contig_name -> BlastHitSummary
    """
    # Collect alignment intervals per query
    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    query_best_hit: Dict[str, Tuple[str, Optional[int], float]] = {}  # (subject, taxid, evalue)

    if not blast_path.exists() or blast_path.stat().st_size == 0:
        return {}

    with blast_path.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            qseqid = fields[0]
            sseqid = fields[1]
            try:
                qstart = int(fields[6])
                qend = int(fields[7])
                evalue = float(fields[10])
            except ValueError:
                continue

            # Parse taxid if present (column 13, may contain multiple semicolon-separated)
            taxid: Optional[int] = None
            if len(fields) >= 13 and fields[12]:
                taxid_str = fields[12].split(";")[0]
                try:
                    taxid = int(taxid_str)
                except ValueError:
                    taxid = None

            # Ensure qstart < qend
            if qstart > qend:
                qstart, qend = qend, qstart

            query_intervals[qseqid].append((qstart, qend))

            # Track best hit by evalue
            if qseqid not in query_best_hit or evalue < query_best_hit[qseqid][2]:
                query_best_hit[qseqid] = (sseqid, taxid, evalue)

    # Compute merged coverage for each query
    results: Dict[str, BlastHitSummary] = {}
    for qseqid, intervals in query_intervals.items():
        _, total_bp = merge_intervals(intervals)
        qlen = query_lengths.get(qseqid, 0)
        coverage = (total_bp / qlen) if qlen > 0 else 0.0

        best_subject, best_taxid, best_evalue = query_best_hit.get(qseqid, ("", None, float("inf")))

        results[qseqid] = BlastHitSummary(
            contig_name=qseqid,
            total_coverage=coverage,
            best_hit_subject=best_subject,
            best_hit_taxid=best_taxid,
            best_hit_evalue=best_evalue,
            total_aligned_bp=total_bp,
        )

    return results
