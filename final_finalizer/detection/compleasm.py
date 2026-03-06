#!/usr/bin/env python3
"""Compleasm (BUSCO completeness) evaluation for final_finalizer.

Runs compleasm on classified FASTA subsets to provide reference-independent
gene completeness metrics based on BUSCO orthologs.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Optional

from final_finalizer.models import CompleasmResult
from final_finalizer.utils.io_utils import file_exists_and_valid, have_exe
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("compleasm")


def parse_compleasm_summary(summary_path: Path) -> Optional[CompleasmResult]:
    """Parse a compleasm summary.txt file into a CompleasmResult.

    Expected format::

        ## lineage: eukaryota_odb12
        S:85.27%, 110
        D:6.98%, 9
        F:3.88%, 5
        I:0.78%, 1
        M:3.10%, 4
        N:129

    Args:
        summary_path: Path to compleasm summary.txt output.

    Returns:
        CompleasmResult or None if the file cannot be parsed.
    """
    if not file_exists_and_valid(summary_path):
        logger.warning(f"Compleasm summary not found: {summary_path}")
        return None

    text = summary_path.read_text()

    # Parse lineage
    lineage_match = re.search(r"## lineage:\s*(\S+)", text)
    lineage = lineage_match.group(1) if lineage_match else "unknown"

    # Parse category lines: X:pct%, count
    categories = {}
    for code in ("S", "D", "F", "I", "M"):
        match = re.search(rf"{code}:([\d.]+)%,\s*(\d+)", text)
        if match:
            categories[code] = (float(match.group(1)), int(match.group(2)))
        else:
            categories[code] = (0.0, 0)

    # Parse N (total)
    n_match = re.search(r"N:(\d+)", text)
    n_total = int(n_match.group(1)) if n_match else 0

    return CompleasmResult(
        lineage=lineage,
        n_total=n_total,
        n_single=categories["S"][1],
        n_duplicated=categories["D"][1],
        n_fragmented=categories["F"][1],
        n_interspersed=categories["I"][1],
        n_missing=categories["M"][1],
        pct_single=categories["S"][0],
        pct_duplicated=categories["D"][0],
        pct_fragmented=categories["F"][0],
        pct_interspersed=categories["I"][0],
        pct_missing=categories["M"][0],
        summary_path=summary_path,
    )


def run_compleasm(
    fasta: Path,
    output_dir: Path,
    lineage: str,
    threads: int,
    library_path: Optional[str] = None,
    compleasm_exe: Optional[str] = None,
    resource_spec=None,
) -> Optional[CompleasmResult]:
    """Run compleasm on a FASTA file and return parsed results.

    Args:
        fasta: Input FASTA file.
        output_dir: Output directory for compleasm.
        lineage: BUSCO lineage name (e.g., "eukaryota", "viridiplantae").
        threads: Number of threads.
        library_path: Path to pre-downloaded lineage files (optional).
        compleasm_exe: Path to compleasm executable (e.g., from a separate
            conda environment). If None, uses ``compleasm`` from PATH.
        resource_spec: ResourceSpec for SLURM submission (unused locally).

    Returns:
        CompleasmResult or None if compleasm is unavailable or fails.
    """
    exe = compleasm_exe or "compleasm"
    if compleasm_exe:
        if not Path(compleasm_exe).is_file():
            logger.warning(f"compleasm not found at {compleasm_exe}, skipping BUSCO evaluation")
            return None
    elif not have_exe("compleasm"):
        logger.warning("compleasm not found in PATH, skipping BUSCO evaluation")
        return None

    if not file_exists_and_valid(fasta):
        logger.info(f"Skipping compleasm: input FASTA empty or missing ({fasta.name})")
        return None

    summary_path = output_dir / "summary.txt"

    # Re-use cached result if available
    if file_exists_and_valid(summary_path):
        logger.info(f"Reusing cached compleasm result: {summary_path}")
        return parse_compleasm_summary(summary_path)

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        exe, "run",
        "-a", str(fasta),
        "-o", str(output_dir),
        "-l", lineage,
        "-t", str(threads),
    ]
    if library_path:
        cmd.extend(["-L", library_path])

    logger.info(f"Running compleasm ({lineage}) on {fasta.name}")

    err_path = output_dir / "compleasm.err"
    try:
        with err_path.open("w") as err_fh:
            ret = subprocess.run(cmd, capture_output=False, stderr=err_fh, stdout=subprocess.DEVNULL, check=False)
        if ret.returncode != 0:
            logger.warning(f"compleasm failed (exit {ret.returncode}); see {err_path}")
            return None
    except Exception as e:
        logger.warning(f"compleasm error: {e}")
        return None

    return parse_compleasm_summary(summary_path)
