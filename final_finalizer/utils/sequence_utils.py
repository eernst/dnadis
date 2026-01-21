#!/usr/bin/env python3
"""
Sequence utilities for final_finalizer.

Contains functions for FASTA reading, writing, and DNA sequence manipulation.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Set

from final_finalizer.utils.io_utils import open_maybe_gzip


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
    import sys

    # Check file size and warn for large files
    try:
        file_size_gb = fasta_path.stat().st_size / (1024**3)
        if file_size_gb > warn_size_gb:
            print(
                f"[warn] Loading {file_size_gb:.1f} GB FASTA into memory. "
                "This may use significant RAM.",
                file=sys.stderr,
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
