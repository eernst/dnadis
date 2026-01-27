#!/usr/bin/env python3
"""
Telomere detection module for final_finalizer.

Detects telomeric repeats at contig ends to identify assembly-complete chromosomes.
Telomere presence indicates the contig has actual chromosome ends, which takes
precedence over reference coverage for full-length classification.

Telomere motifs:
- Plants (default): TTTAGGG / CCCTAAA
- Vertebrates: TTAGGG / CCCTAA
- Custom via --telomere-motif

Detection algorithm:
1. Search first/last window_size bp of each contig
2. Look for tandem telomere repeats (min_repeats consecutive matches)
3. Return TelomereResult with 5' and 3' telomere status
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional

from final_finalizer.models import TelomereResult
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.sequence_utils import read_fasta_sequences, reverse_complement

logger = get_logger("telomere")

# Default telomere motifs
DEFAULT_TELOMERE_MOTIF = "TTTAGGG"  # Plant telomere (forward strand)
VERTEBRATE_TELOMERE_MOTIF = "TTAGGG"  # Vertebrate telomere


def _count_telomere_repeats(
    sequence: str,
    motif: str,
    min_repeats: int = 3,
) -> int:
    """Count tandem telomere repeats in a sequence.

    Looks for consecutive occurrences of the telomere motif.

    Args:
        sequence: DNA sequence to search (uppercase)
        motif: Telomere motif to search for (e.g., "TTTAGGG")
        min_repeats: Minimum consecutive repeats to count as telomeric

    Returns:
        Number of consecutive telomere repeats found, or 0 if < min_repeats
    """
    if not sequence or not motif:
        return 0

    sequence = sequence.upper()
    motif = motif.upper()
    rc_motif = reverse_complement(motif)

    # Build regex for tandem repeats (both forward and reverse complement)
    # Match at least min_repeats consecutive copies
    fwd_pattern = f"({re.escape(motif)}){{{min_repeats},}}"
    rev_pattern = f"({re.escape(rc_motif)}){{{min_repeats},}}"

    max_count = 0

    # Search for forward motif repeats
    for match in re.finditer(fwd_pattern, sequence):
        count = len(match.group(0)) // len(motif)
        max_count = max(max_count, count)

    # Search for reverse complement motif repeats
    for match in re.finditer(rev_pattern, sequence):
        count = len(match.group(0)) // len(rc_motif)
        max_count = max(max_count, count)

    return max_count if max_count >= min_repeats else 0


def detect_telomeres_single(
    sequence: str,
    motif: str = DEFAULT_TELOMERE_MOTIF,
    window_size: int = 10000,
    min_repeats: int = 3,
) -> TelomereResult:
    """Detect telomeres at both ends of a single contig.

    Args:
        sequence: Full contig sequence
        motif: Telomere repeat motif (e.g., "TTTAGGG")
        window_size: Size of window at each end to search (bp)
        min_repeats: Minimum consecutive repeats to call telomere

    Returns:
        TelomereResult with detection status for both ends
    """
    if not sequence:
        return TelomereResult(
            has_5p_telomere=False,
            has_3p_telomere=False,
            telomere_5p_count=0,
            telomere_3p_count=0,
        )

    seq_len = len(sequence)
    window = min(window_size, seq_len // 2)  # Don't overlap windows

    # Extract 5' and 3' ends
    five_prime = sequence[:window]
    three_prime = sequence[-window:] if seq_len > window else sequence[window:]

    # Count telomere repeats at each end
    count_5p = _count_telomere_repeats(five_prime, motif, min_repeats)
    count_3p = _count_telomere_repeats(three_prime, motif, min_repeats)

    return TelomereResult(
        has_5p_telomere=count_5p >= min_repeats,
        has_3p_telomere=count_3p >= min_repeats,
        telomere_5p_count=count_5p,
        telomere_3p_count=count_3p,
    )


def detect_telomeres(
    query_fasta: Path,
    contig_names: Optional[set] = None,
    motif: str = DEFAULT_TELOMERE_MOTIF,
    window_size: int = 10000,
    min_repeats: int = 3,
) -> Dict[str, TelomereResult]:
    """Detect telomeres in all contigs from a FASTA file.

    Args:
        query_fasta: Path to query FASTA file
        contig_names: Optional set of contig names to analyze (None = all)
        motif: Telomere repeat motif
        window_size: Size of window at each end to search (bp)
        min_repeats: Minimum consecutive repeats to call telomere

    Returns:
        Dict mapping contig name -> TelomereResult
    """
    logger.info(f"Detecting telomeres (motif={motif}, window={window_size}bp, min_repeats={min_repeats})")

    sequences = read_fasta_sequences(query_fasta)
    results: Dict[str, TelomereResult] = {}

    n_with_5p = 0
    n_with_3p = 0
    n_with_both = 0

    for contig_name, sequence in sequences.items():
        if contig_names is not None and contig_name not in contig_names:
            continue

        result = detect_telomeres_single(
            sequence=sequence,
            motif=motif,
            window_size=window_size,
            min_repeats=min_repeats,
        )
        results[contig_name] = result

        if result.has_5p_telomere:
            n_with_5p += 1
        if result.has_3p_telomere:
            n_with_3p += 1
        if result.has_5p_telomere and result.has_3p_telomere:
            n_with_both += 1

    n_total = len(results)
    if n_total > 0:
        logger.info(
            f"Telomere detection: {n_with_5p}/{n_total} with 5' telomere, "
            f"{n_with_3p}/{n_total} with 3' telomere, "
            f"{n_with_both}/{n_total} with both"
        )

    return results
