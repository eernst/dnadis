#!/usr/bin/env python3
"""
Sequence utilities for final_finalizer.

Contains functions for FASTA reading, writing, and DNA sequence manipulation.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.utils.io_utils import open_maybe_gzip
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("sequence_utils")


# ----------------------------
# FASTA sequence utilities
# ----------------------------
_COMPLEMENT = str.maketrans("ATCGatcgNnRrYySsWwKkMmBbDdHhVv", "TAGCtagcNnYyRrSsWwMmKkVvHhDdBb")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def read_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """
    Return dict: seqname -> length for a FASTA that may be gzipped.
    """
    lengths: dict[str, int] = {}
    name: Optional[str] = None
    seqlen = 0

    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seqlen
                name = line[1:].strip().split()[0]
                seqlen = 0
            else:
                seqlen += len(line.strip())

    if name is not None:
        lengths[name] = seqlen

    return lengths


def read_fasta_sequences(fasta_path: Path, warn_size_gb: float = 2.0) -> Dict[str, str]:
    """Read FASTA file and return dict of name -> sequence.

    Note: This loads the entire file into memory. For very large genomes,
    consider using write_filtered_fasta() which streams data instead.

    Args:
        fasta_path: Path to FASTA file (plain or gzipped)
        warn_size_gb: Warn if file exceeds this size in GB (default 2.0)

    Returns:
        Dict mapping sequence name to sequence string
    """
    # Check file size and warn for large files
    try:
        file_size_gb = fasta_path.stat().st_size / (1024**3)
        if file_size_gb > warn_size_gb:
            logger.warning(
                f"Loading {file_size_gb:.1f} GB FASTA into memory. "
                "This may use significant RAM."
            )
    except OSError:
        pass  # File might be gzipped, size check less meaningful

    sequences: Dict[str, str] = {}
    name: Optional[str] = None
    seq_parts: List[str] = []

    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(seq_parts)
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if name is not None:
        sequences[name] = "".join(seq_parts)

    return sequences


def write_fasta(sequences: Dict[str, str], output_path: Path, wrap: int = 80) -> None:
    """Write sequences to FASTA file with line wrapping."""
    with output_path.open("w") as out:
        for name, seq in sequences.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i + wrap] + "\n")


def write_filtered_fasta(
    input_fasta: Path,
    output_fasta: Path,
    include_contigs: Set[str],
    wrap: int = 80,
) -> None:
    """Write a filtered FASTA containing only specified contigs.

    Streams through the input to avoid loading entire genome into memory.
    """
    with open_maybe_gzip(input_fasta, "rt") as fh, output_fasta.open("w") as out:
        writing = False
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                writing = name in include_contigs
                if writing:
                    out.write(line)
            elif writing:
                out.write(line)


def is_hifiasm_circular(contig_name: str) -> bool:
    """Check if contig name matches hifiasm circular pattern (ptg*c)."""
    return bool(re.match(r"ptg\d+c", contig_name))


def calculate_gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as fraction (0.0-1.0). Returns 0.0 for empty sequences
        or sequences with only ambiguous bases.
    """
    seq_upper = seq.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    # Only count unambiguous bases (A, T, G, C)
    total = (
        seq_upper.count("A")
        + seq_upper.count("T")
        + seq_upper.count("G")
        + seq_upper.count("C")
    )
    if total == 0:
        return 0.0
    return gc_count / total


def calculate_gc_content_fasta(fasta_path: Path) -> Dict[str, float]:
    """Calculate GC content for all sequences in a FASTA file.

    Streams through the file to avoid loading entire genome into memory.

    Args:
        fasta_path: Path to FASTA file (plain or gzipped)

    Returns:
        Dict mapping sequence name to GC content (0.0-1.0)
    """
    gc_contents: Dict[str, float] = {}
    name: Optional[str] = None
    gc_count = 0
    total_count = 0

    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence's GC content
                if name is not None:
                    gc_contents[name] = (gc_count / total_count) if total_count > 0 else 0.0
                name = line[1:].strip().split()[0]
                gc_count = 0
                total_count = 0
            else:
                seq = line.strip().upper()
                gc_count += seq.count("G") + seq.count("C")
                total_count += seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")

    # Save last sequence
    if name is not None:
        gc_contents[name] = (gc_count / total_count) if total_count > 0 else 0.0

    return gc_contents


def calculate_gc_stats(gc_contents: Dict[str, float]) -> tuple[float, float]:
    """Calculate mean and standard deviation of GC contents.

    Args:
        gc_contents: Dict of sequence name to GC content

    Returns:
        Tuple of (mean, std_dev). Returns (0.0, 0.0) if empty.
    """
    if not gc_contents:
        return 0.0, 0.0

    values = list(gc_contents.values())
    n = len(values)
    mean = sum(values) / n

    if n < 2:
        return mean, 0.0

    variance = sum((x - mean) ** 2 for x in values) / (n - 1)
    std_dev = variance ** 0.5

    return mean, std_dev


# ----------------------------
# GC content caching
# ----------------------------


def write_gc_content_tsv(
    tsv_path: Path,
    metadata_path: Path,
    ref_gc: Dict[str, float],
    qry_gc: Dict[str, float],
    ref_path: Path,
    qry_path: Path,
) -> None:
    """Write GC content cache file with metadata for validation.

    Args:
        tsv_path: Path for TSV output
        metadata_path: Path for JSON metadata
        ref_gc: Reference GC content dict
        qry_gc: Query GC content dict
        ref_path: Path to reference FASTA (for cache validation)
        qry_path: Path to query FASTA (for cache validation)
    """
    # Write TSV with source column
    with tsv_path.open("w") as f:
        f.write("source\tseqid\tgc_content\n")
        for seqid, gc in ref_gc.items():
            f.write(f"reference\t{seqid}\t{gc:.6f}\n")
        for seqid, gc in qry_gc.items():
            f.write(f"query\t{seqid}\t{gc:.6f}\n")

    # Write metadata for cache validation
    ref_resolved = ref_path.resolve()
    qry_resolved = qry_path.resolve()

    metadata = {
        "ref_path": str(ref_resolved),
        "ref_name": ref_path.name,
        "ref_mtime": ref_resolved.stat().st_mtime if ref_resolved.exists() else 0,
        "ref_size": ref_resolved.stat().st_size if ref_resolved.exists() else 0,
        "query_path": str(qry_resolved),
        "query_name": qry_path.name,
        "query_mtime": qry_resolved.stat().st_mtime if qry_resolved.exists() else 0,
        "query_size": qry_resolved.stat().st_size if qry_resolved.exists() else 0,
        "ref_seqs": len(ref_gc),
        "query_seqs": len(qry_gc),
    }

    with metadata_path.open("w") as f:
        json.dump(metadata, f, indent=2)


def read_gc_content_tsv(tsv_path: Path) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Read cached GC content from TSV file.

    Args:
        tsv_path: Path to GC content TSV

    Returns:
        Tuple of (ref_gc, qry_gc) dicts
    """
    ref_gc: Dict[str, float] = {}
    qry_gc: Dict[str, float] = {}

    with tsv_path.open() as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue
            source, seqid, gc_str = parts
            gc = float(gc_str)
            if source == "reference":
                ref_gc[seqid] = gc
            elif source == "query":
                qry_gc[seqid] = gc

    return ref_gc, qry_gc


def check_gc_cache(
    tsv_path: Path,
    metadata_path: Path,
    ref_path: Path,
    qry_path: Path,
) -> bool:
    """Check if GC content cache is valid.

    Validates that:
    - TSV and metadata files exist
    - Reference and query paths match
    - File modification times haven't changed
    - File sizes haven't changed

    Args:
        tsv_path: Path to GC content TSV
        metadata_path: Path to metadata JSON
        ref_path: Current reference FASTA path
        qry_path: Current query FASTA path

    Returns:
        True if cache is valid and can be reused
    """
    if not tsv_path.exists() or not metadata_path.exists():
        return False

    try:
        with metadata_path.open() as f:
            metadata = json.load(f)
    except (json.JSONDecodeError, OSError):
        return False

    ref_resolved = ref_path.resolve()
    qry_resolved = qry_path.resolve()

    # Check reference file
    if metadata.get("ref_path") != str(ref_resolved):
        logger.info(
            f"GC cache invalid: ref path changed "
            f"(cached={metadata.get('ref_name', 'unknown')}, current={ref_path.name})"
        )
        return False

    try:
        ref_stat = ref_resolved.stat()
        if metadata.get("ref_mtime", 0) != ref_stat.st_mtime:
            logger.info(f"GC cache invalid: ref file modified ({ref_path.name})")
            return False
        if metadata.get("ref_size", 0) != ref_stat.st_size:
            logger.info(f"GC cache invalid: ref file size changed ({ref_path.name})")
            return False
    except OSError:
        return False

    # Check query file
    if metadata.get("query_path") != str(qry_resolved):
        logger.info(
            f"GC cache invalid: query path changed "
            f"(cached={metadata.get('query_name', 'unknown')}, current={qry_path.name})"
        )
        return False

    try:
        qry_stat = qry_resolved.stat()
        if metadata.get("query_mtime", 0) != qry_stat.st_mtime:
            logger.info(f"GC cache invalid: query file modified ({qry_path.name})")
            return False
        if metadata.get("query_size", 0) != qry_stat.st_size:
            logger.info(f"GC cache invalid: query file size changed ({qry_path.name})")
            return False
    except OSError:
        return False

    return True
