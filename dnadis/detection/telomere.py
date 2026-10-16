#!/usr/bin/env python3
"""
Telomere detection module for dnadis.

Detects telomeric repeats at contig ends to identify assembly-complete chromosomes.
Telomere presence indicates the contig has actual chromosome ends, which takes
precedence over reference coverage for full-length classification.

Telomere motifs:
- Plants (default): TTTAGGG / CCCTAAA
- Vertebrates: TTAGGG / CCCTAA
- Custom via --telomere-motif

Detection algorithm:
1. Search first/last window_size bp of each contig
2. Find runs of tandem telomere repeats (min_repeats consecutive matches)
3. Call a telomere when the runs total at least min_bp; arrays are
   kilobases long, while a few repeats occur by chance in other simple
   sequence near contig ends
4. Return TelomereResult with 5' and 3' telomere status
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional, Tuple

from dnadis.models import TelomereResult
from dnadis.utils.logging_config import get_logger
from dnadis.utils.sequence_utils import read_fasta_sequences, reverse_complement

logger = get_logger("telomere")

# Default telomere motifs
DEFAULT_TELOMERE_MOTIF = "TTTAGGG"  # Plant telomere (forward strand)
VERTEBRATE_TELOMERE_MOTIF = "TTAGGG"  # Vertebrate telomere


def _telomere_runs(
    sequence: str,
    motif: str,
    min_repeats: int = 3,
) -> Tuple[int, int]:
    """Measure tandem telomere repeat runs in a sequence.

    A run is at least ``min_repeats`` consecutive copies of the motif or its
    reverse complement.  Arrays carry variant repeats that break exact runs,
    so the bp of all runs is summed rather than taking the longest one.

    Args:
        sequence: DNA sequence to search
        motif: Telomere motif to search for (e.g., "TTTAGGG")
        min_repeats: Minimum consecutive repeats in a run

    Returns:
        Tuple of (total bp in runs, longest run in repeats)
    """
    if not sequence or not motif:
        return 0, 0

    sequence = sequence.upper()
    motif = motif.upper()
    rc_motif = reverse_complement(motif)

    total_bp = 0
    max_count = 0
    for m in (motif, rc_motif):
        for match in re.finditer(f"(?:{re.escape(m)}){{{min_repeats},}}", sequence):
            run_bp = match.end() - match.start()
            total_bp += run_bp
            max_count = max(max_count, run_bp // len(m))

    return total_bp, max_count


def detect_telomeres_single(
    sequence: str,
    motif: str = DEFAULT_TELOMERE_MOTIF,
    window_size: int = 10000,
    min_repeats: int = 3,
    min_bp: int = 200,
) -> TelomereResult:
    """Detect telomeres at both ends of a single contig.

    Args:
        sequence: Full contig sequence
        motif: Telomere repeat motif (e.g., "TTTAGGG")
        window_size: Size of window at each end to search (bp)
        min_repeats: Minimum consecutive repeats in a run
        min_bp: Minimum total bp of runs at an end to call a telomere

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

    bp_5p, count_5p = _telomere_runs(five_prime, motif, min_repeats)
    bp_3p, count_3p = _telomere_runs(three_prime, motif, min_repeats)

    return TelomereResult(
        has_5p_telomere=bp_5p >= min_bp,
        has_3p_telomere=bp_3p >= min_bp,
        telomere_5p_count=count_5p,
        telomere_3p_count=count_3p,
        telomere_5p_bp=bp_5p,
        telomere_3p_bp=bp_3p,
    )


def detect_telomeres(
    query_fasta: Path,
    contig_names: Optional[set] = None,
    motif: str = DEFAULT_TELOMERE_MOTIF,
    window_size: int = 10000,
    min_repeats: int = 3,
    min_bp: int = 200,
) -> Dict[str, TelomereResult]:
    """Detect telomeres in all contigs from a FASTA file.

    Args:
        query_fasta: Path to query FASTA file
        contig_names: Optional set of contig names to analyze (None = all)
        motif: Telomere repeat motif
        window_size: Size of window at each end to search (bp)
        min_repeats: Minimum consecutive repeats in a run
        min_bp: Minimum total bp of runs at an end to call a telomere

    Returns:
        Dict mapping contig name -> TelomereResult
    """
    logger.info(
        f"Detecting telomeres (motif={motif}, window={window_size}bp, "
        f"min_repeats={min_repeats}, min_bp={min_bp})"
    )

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
            min_bp=min_bp,
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
