#!/usr/bin/env python3
"""Plotting functions using R/ggplot2."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, Optional

from final_finalizer.utils.logging_config import get_logger

logger = get_logger("plotting")


def have_rscript() -> bool:
    try:
        subprocess.run(["Rscript", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception:
        return False


def _read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8") as fh:
        return fh.read()


def _esc(s) -> str:
    """Escape path for R script (convert backslashes to forward slashes)."""
    return str(s).replace("\\", "/")


def run_plot(
    summary_tsv: Path,
    ref_lengths_tsv: Path,
    segments_tsv: Path,
    chain_summary_tsv: Path,
    macro_blocks_tsv: Path,
    outprefix: Path,
    chr_like_minlen: int,
    plot_title_suffix: str,
    plot_html: bool,
):
    if not have_rscript():
        logger.warning("--plot specified, but Rscript not found in PATH; skipping plot.")
        return

    plot_pdf = Path(str(outprefix) + ".chromosome_overview.pdf")
    plot_html_path = Path(str(outprefix) + ".chromosome_overview.html")
    r_script_path = Path(str(outprefix) + ".chromosome_overview.R")

    # Template expected alongside this Python file
    tmpl_path = Path(__file__).resolve().with_name("visual_reporting.tmpl.R")
    if not tmpl_path.exists():
        raise FileNotFoundError(
            f"Missing R template: {tmpl_path}\n"
            "Create it next to the Python script (see instructions)."
        )

    tmpl = _read_text(tmpl_path)

    # Very simple placeholder replacement (safe + robust)
    def esc(s: Path | str) -> str:
        # ensure backslashes don't accidentally escape quotes on Windows
        return str(s).replace("\\", "/")

    filled = (
        tmpl.replace("__SUMMARY__", esc(summary_tsv))
        .replace("__REF__", esc(ref_lengths_tsv))
        .replace("__SEGMENTS__", esc(segments_tsv))
        .replace("__EVIDENCE__", esc(chain_summary_tsv))
        .replace("__MACRO__", esc(macro_blocks_tsv))
        .replace("__OUTPDF__", esc(plot_pdf))
        .replace("__OUTHTML__", esc(plot_html_path))
        .replace("__PLOTHTML__", "TRUE" if plot_html else "FALSE")
        .replace("__CHRLIKE__", str(int(chr_like_minlen)))
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    logger.info(f"Running Rscript to generate plot: {plot_pdf}")
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; plot not generated.")
    else:
        logger.done(f"Plot written to: {plot_pdf}")
        if plot_html:
            logger.done(f"Plot written to: {plot_html_path}")


def run_depth_plot(
    summary_tsv: Path,
    outprefix: Path,
    plot_title_suffix: str,
    plot_html: bool,
) -> None:
    """Generate depth overview plot if depth data is available.

    Creates:
    - {outprefix}.depth_overview.pdf: Static PDF with depth visualizations
    - {outprefix}.depth_overview.html: Interactive HTML (if plot_html=True)

    Args:
        summary_tsv: Path to contig_summary.tsv with depth columns
        outprefix: Output prefix for generated files
        plot_title_suffix: Suffix for plot titles
        plot_html: Whether to generate interactive HTML version
    """
    if not have_rscript():
        logger.warning("--plot specified, but Rscript not found in PATH; skipping depth plot.")
        return

    # Check if summary file has depth data
    try:
        with summary_tsv.open("r") as fh:
            header = fh.readline().strip().split("\t")
            if "depth_mean" not in header:
                logger.info("No depth data in summary; skipping depth plot.")
                return
    except Exception as e:
        logger.warning(f"Could not check summary file for depth data: {e}")
        return

    plot_pdf = Path(str(outprefix) + ".depth_overview.pdf")
    plot_html_path = Path(str(outprefix) + ".depth_overview.html")
    r_script_path = Path(str(outprefix) + ".depth_overview.R")

    # Template expected alongside this Python file
    tmpl_path = Path(__file__).resolve().with_name("depth_overview.tmpl.R")
    if not tmpl_path.exists():
        logger.warning(f"Missing R template for depth plots: {tmpl_path}")
        return

    tmpl = _read_text(tmpl_path)

    filled = (
        tmpl.replace("__SUMMARY__", _esc(summary_tsv))
        .replace("__OUTPDF__", _esc(plot_pdf))
        .replace("__OUTHTML__", _esc(plot_html_path))
        .replace("__PLOTHTML__", "TRUE" if plot_html else "FALSE")
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    logger.info(f"Running Rscript to generate depth plot: {plot_pdf}")
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; depth plot not generated.")
    else:
        logger.done(f"Depth plot written to: {plot_pdf}")
        if plot_html:
            logger.done(f"Depth plot written to: {plot_html_path}")


def run_contaminant_plot(
    contaminants_tsv: Path,
    outprefix: Path,
    plot_title_suffix: str,
    plot_html: bool,
) -> None:
    """Generate contamination treemap showing phylogenetic breakdown.

    Creates a treemap showing the taxonomic hierarchy as nested rectangles
    (Kingdom > Family > Genus > Species) with area proportional to total
    contamination span bp.

    Note: Coverage filtering is applied upstream in Python before writing the TSV.

    Args:
        contaminants_tsv: Path to contaminants.tsv with taxonomic lineage
        outprefix: Output prefix for generated files
        plot_title_suffix: Suffix for plot titles
        plot_html: Whether to generate interactive HTML version
    """
    if not have_rscript():
        logger.warning("Rscript not found in PATH; skipping contaminant plot.")
        return

    # Check if file exists and has content
    if not contaminants_tsv.exists():
        logger.warning(f"Contaminants TSV not found: {contaminants_tsv}")
        return

    # Check for empty file (just header or no data)
    try:
        with contaminants_tsv.open("r") as fh:
            lines = fh.readlines()
            if len(lines) <= 1:
                logger.info("No contaminants in TSV; skipping contaminant plot.")
                return
    except Exception as e:
        logger.warning(f"Could not read contaminants TSV: {e}")
        return

    plot_pdf = Path(str(outprefix) + ".contaminant_treemap.pdf")
    plot_html_path = Path(str(outprefix) + ".contaminant_treemap.html")
    r_script_path = Path(str(outprefix) + ".contaminant_treemap.R")

    # Template expected alongside this Python file
    tmpl_path = Path(__file__).resolve().with_name("contaminant_treemap.tmpl.R")
    if not tmpl_path.exists():
        logger.warning(f"Missing R template for contaminant plots: {tmpl_path}")
        return

    tmpl = _read_text(tmpl_path)

    filled = (
        tmpl.replace("__CONTAMINANTS_TSV__", _esc(contaminants_tsv))
        .replace("__OUTPDF__", _esc(plot_pdf))
        .replace("__OUTHTML__", _esc(plot_html_path))
        .replace("__PLOTHTML__", "TRUE" if plot_html else "FALSE")
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    logger.info(f"Running Rscript to generate contaminant treemap: {plot_pdf}")
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; contaminant plot not generated.")
    else:
        logger.done(f"Contaminant treemap written to: {plot_pdf}")


def run_contaminant_bandage_plot(
    contaminants_tsv: Path,
    outprefix: Path,
    plot_title_suffix: str,
) -> None:
    """Generate Bandage-like contamination plot showing individual contigs.

    Shows individual contigs as geometric shapes:
    - Circular contigs (name ending in "c"): Drawn as rings/donuts
    - Linear contigs (name ending in "l"): Drawn as pills (stadium shapes)

    Size uses power-scale areas so small contigs remain visible.
    Contigs are colored by family and ordered by decreasing depth.

    Note: Coverage filtering is applied upstream in Python before writing the TSV.
    Requires depth data (--reads) for depth-ordered layout.

    Args:
        contaminants_tsv: Path to contaminants.tsv with taxonomic lineage
        outprefix: Output prefix for generated files
        plot_title_suffix: Suffix for plot titles
    """
    if not have_rscript():
        logger.warning("Rscript not found; skipping contaminant bandage plot.")
        return

    if not contaminants_tsv.exists():
        return

    # Check for empty file and depth data availability
    try:
        with contaminants_tsv.open("r") as fh:
            lines = fh.readlines()
            if len(lines) <= 1:
                return
            header = lines[0].strip().split("\t")
            if "depth_mean" not in header:
                logger.info("No depth data in contaminants TSV; skipping bandage plot (requires --reads).")
                return
    except Exception:
        return

    plot_pdf = Path(str(outprefix) + ".contaminant_bandage.pdf")
    r_script_path = Path(str(outprefix) + ".contaminant_bandage.R")

    tmpl_path = Path(__file__).resolve().with_name("contaminant_bandage.tmpl.R")
    if not tmpl_path.exists():
        logger.warning(f"Missing R template: {tmpl_path}")
        return

    tmpl = _read_text(tmpl_path)
    filled = (
        tmpl.replace("__CONTAMINANTS_TSV__", _esc(contaminants_tsv))
        .replace("__OUTPDF__", _esc(plot_pdf))
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    logger.info(f"Running Rscript to generate contaminant bandage plot: {plot_pdf}")
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; bandage plot not generated.")
    else:
        logger.done(f"Contaminant bandage plot written to: {plot_pdf}")
