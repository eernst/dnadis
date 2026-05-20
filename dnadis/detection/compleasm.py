#!/usr/bin/env python3
"""Compleasm (BUSCO completeness) evaluation for dnadis.

Runs compleasm on classified FASTA subsets to provide reference-independent
gene completeness metrics based on BUSCO orthologs.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Optional

from collections import defaultdict
from typing import Iterable, List, Tuple

from dnadis.models import CompleasmResult, ContigClassification
from dnadis.utils.io_utils import file_exists_and_valid, have_exe
from dnadis.utils.logging_config import get_logger

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

    lineage_subdir = summary_path.parent / lineage if lineage != "unknown" else None
    full_table_path: Optional[Path] = None
    translated_protein_path: Optional[Path] = None
    if lineage_subdir is not None:
        candidate_full = lineage_subdir / "full_table_busco_format.tsv"
        if candidate_full.exists():
            full_table_path = candidate_full
        candidate_prot = lineage_subdir / "translated_protein.fasta"
        if candidate_prot.exists():
            translated_protein_path = candidate_prot

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
        full_table_path=full_table_path,
        translated_protein_path=translated_protein_path,
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

    # When using --compleasm-path, prepend its parent directory to PATH so
    # compleasm can find sibling tools (hmmsearch, miniprot) from the same
    # conda environment.
    env = None
    if compleasm_exe:
        env = os.environ.copy()
        env["PATH"] = str(Path(compleasm_exe).resolve().parent) + os.pathsep + env.get("PATH", "")

    err_path = output_dir / "compleasm.err"
    try:
        with err_path.open("w") as err_fh:
            ret = subprocess.run(cmd, capture_output=False, stderr=err_fh, stdout=subprocess.DEVNULL, check=False, env=env)
        if ret.returncode != 0:
            logger.warning(f"compleasm failed (exit {ret.returncode}); see {err_path}")
            return None
    except Exception as e:
        logger.warning(f"compleasm error: {e}")
        return None

    return parse_compleasm_summary(summary_path)


def _write_derived_summary(
    summary_path: Path,
    lineage: str,
    n_single: int,
    n_duplicated: int,
    n_fragmented: int,
    n_missing: int,
    n_total: int,
) -> None:
    """Write a compleasm-format summary.txt for a derived subset.

    Shape matches compleasm's own output so it parses cleanly with
    :func:`parse_compleasm_summary` and so the same downstream report
    paths consume it without special casing.  The ``I`` (interspaced)
    line is always emitted as zero because compleasm's
    ``full_table_busco_format.tsv`` does not propagate that status.
    """
    if n_total <= 0:
        pct = {"S": 0.0, "D": 0.0, "F": 0.0, "I": 0.0, "M": 0.0}
    else:
        pct = {
            "S": 100.0 * n_single / n_total,
            "D": 100.0 * n_duplicated / n_total,
            "F": 100.0 * n_fragmented / n_total,
            "I": 0.0,
            "M": 100.0 * n_missing / n_total,
        }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        f"## lineage: {lineage}\n"
        f"S:{pct['S']:.2f}%, {n_single}\n"
        f"D:{pct['D']:.2f}%, {n_duplicated}\n"
        f"F:{pct['F']:.2f}%, {n_fragmented}\n"
        f"I:{pct['I']:.2f}%, 0\n"
        f"M:{pct['M']:.2f}%, {n_missing}\n"
        f"N:{n_total}\n"
    )


def derive_chrs_non_chrs_compleasm(
    full_result: CompleasmResult,
    classifications: Iterable[ContigClassification],
    chrs_out_dir: Path,
    non_chrs_out_dir: Path,
) -> Tuple[Optional[CompleasmResult], Optional[CompleasmResult]]:
    """Derive (compleasm_chrs, compleasm_non_chrs) from a full-assembly run.

    Splits ``full_result.full_table_path`` by joining each row's contig
    against ``classifications`` (chrom_assigned → chrs; everything else
    → non_chrs).  For each split, writes a compleasm-format
    ``summary.txt`` and a filtered ``full_table_busco_format.tsv`` under
    ``{out_dir}/{lineage}/`` so the user can browse the split directly.

    Status counts in derived summaries follow the same rule compleasm
    applies to the whole assembly, but on the per-split row subset:

    * 0 rows for this BUSCO in the split → ``Missing``
    * 1 non-Fragmented row → ``Complete``
    * ≥ 2 non-Fragmented rows → ``Duplicated``
    * only Fragmented rows → ``Fragmented``

    Returns ``(chrs_result, non_chrs_result)``.  Either can be ``None``
    if the full-table cannot be read.
    """
    if full_result is None or full_result.full_table_path is None:
        return None, None
    if not file_exists_and_valid(full_result.full_table_path):
        return None, None

    contig_to_split: dict = {}
    for clf in classifications:
        cat = "chrs" if clf.classification == "chrom_assigned" else "non_chrs"
        contig_to_split[clf.original_name] = cat
        contig_to_split[clf.new_name] = cat

    per_split_status: dict = {"chrs": defaultdict(list), "non_chrs": defaultdict(list)}
    raw_lines: dict = {"chrs": [], "non_chrs": []}

    header = (
        "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength"
    )

    with full_result.full_table_path.open() as fh:
        for line in fh:
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue
            cols = stripped.split("\t")
            if len(cols) < 2:
                continue
            bid, status = cols[0], cols[1]
            if status == "Missing":
                continue
            contig = cols[2] if len(cols) >= 3 else ""
            split = contig_to_split.get(contig)
            if split is None:
                # Hit on a contig dnadis did not classify (e.g. compleasm
                # ran on a FASTA whose contig didn't survive into the
                # per-contig classification list).  Skip rather than
                # misattribute.
                continue
            per_split_status[split][bid].append(status)
            raw_lines[split].append(stripped)

    def _build(split: str, out_dir: Path) -> Optional[CompleasmResult]:
        bid_to_statuses = per_split_status[split]
        n_S = n_D = n_F = 0
        for statuses in bid_to_statuses.values():
            non_frag = [s for s in statuses if s != "Fragmented"]
            if not non_frag:
                n_F += 1
            elif len(non_frag) == 1:
                n_S += 1
            else:
                n_D += 1
        n_present = n_S + n_D + n_F
        n_missing = max(0, full_result.n_total - n_present)

        out_dir.mkdir(parents=True, exist_ok=True)
        lineage_subdir = out_dir / full_result.lineage
        lineage_subdir.mkdir(parents=True, exist_ok=True)
        summary_path = out_dir / "summary.txt"
        ft_path = lineage_subdir / "full_table_busco_format.tsv"

        _write_derived_summary(
            summary_path,
            lineage=full_result.lineage,
            n_single=n_S,
            n_duplicated=n_D,
            n_fragmented=n_F,
            n_missing=n_missing,
            n_total=full_result.n_total,
        )
        with ft_path.open("w") as out_fh:
            out_fh.write(header + "\n")
            for line in raw_lines[split]:
                out_fh.write(line + "\n")
        return parse_compleasm_summary(summary_path)

    return _build("chrs", chrs_out_dir), _build("non_chrs", non_chrs_out_dir)
