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
    rdna_annotations_tsv: Optional[Path] = None,
    agp_tsv: Optional[Path] = None,
):
    if not have_rscript():
        logger.warning("--plot specified, but Rscript not found in PATH; skipping plot.")
        return

    plot_pdf = Path(str(outprefix) + ".chromosome_overview.pdf")
    plot_html_path = Path(str(outprefix) + ".chromosome_overview.html")
    r_script_path = Path(str(outprefix) + ".chromosome_overview.R")

    # Template expected alongside this Python file
    tmpl_path = Path(__file__).resolve().with_name("chromosome_overview.tmpl.R")
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

    # Handle optional rdna_annotations_tsv
    rdna_file_str = esc(rdna_annotations_tsv) if rdna_annotations_tsv and rdna_annotations_tsv.exists() else ""
    agp_file_str = esc(agp_tsv) if agp_tsv and agp_tsv.exists() else ""

    filled = (
        tmpl.replace("__SUMMARY__", esc(summary_tsv))
        .replace("__REF__", esc(ref_lengths_tsv))
        .replace("__SEGMENTS__", esc(segments_tsv))
        .replace("__EVIDENCE__", esc(chain_summary_tsv))
        .replace("__MACRO__", esc(macro_blocks_tsv))
        .replace("__RDNA_ANNOTATIONS__", rdna_file_str)
        .replace("__AGP__", agp_file_str)
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
    ref_lengths_tsv: Path,
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
        ref_lengths_tsv: Path to ref_lengths.tsv for subgenome colors
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
        .replace("__REF__", _esc(ref_lengths_tsv))
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


def run_contaminant_table(
    contaminants_tsv: Path,
    outprefix: Path,
    plot_title_suffix: str,
    top_n: int = 10,
) -> None:
    """Generate gt table showing top contaminants ranked by abundance.

    Creates an HTML table using the gt package showing:
    - Rank by abundance (depth × length) or total length if no depth data
    - Domain indicator (colored badge)
    - Family (colored left border)
    - Binomial species name (italic)
    - Contig count
    - Total Mb with inline bar + spread showing individual contig sizes
    - Mean depth with inline bar + spread (only if depth data available)
    - Mean coverage percentage
    - Circular indicator (●/◐/○)

    The spread values highlight large contigs (potential complete genomes)
    and depth variability across contigs of the same species.

    Note: HTML-only output (CSS gradient bars don't render to PDF).
    Depth column is excluded if no depth data available.

    Args:
        contaminants_tsv: Path to contaminants.tsv with taxonomic lineage
        outprefix: Output prefix for generated files
        plot_title_suffix: Suffix for plot titles
        top_n: Number of top contaminants to display (default 10)
    """
    if not have_rscript():
        logger.warning("Rscript not found; skipping contaminant table.")
        return

    if not contaminants_tsv.exists():
        return

    # Check for empty file
    try:
        with contaminants_tsv.open("r") as fh:
            lines = fh.readlines()
            if len(lines) <= 1:
                return
    except Exception:
        return

    plot_html = Path(str(outprefix) + ".contaminant_table.html")
    r_script_path = Path(str(outprefix) + ".contaminant_table.R")

    tmpl_path = Path(__file__).resolve().with_name("contaminant_table.tmpl.R")
    if not tmpl_path.exists():
        logger.warning(f"Missing R template: {tmpl_path}")
        return

    tmpl = _read_text(tmpl_path)
    filled = (
        tmpl.replace("__CONTAMINANTS_TSV__", _esc(contaminants_tsv))
        .replace("__OUTHTML__", _esc(plot_html))
        .replace("__TOP_N__", str(top_n))
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    logger.info(f"Running Rscript to generate contaminant table: {plot_html}")
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; contaminant table not generated.")
    else:
        if plot_html.exists():
            logger.done(f"Contaminant table written to: {plot_html}")


def run_classification_summary_bar(
    summary_tsv: Path,
    outprefix: Path,
    plot_title_suffix: str,
    plot_html: bool = False,
    assembly_name: str = "",
    reference_name: str = "",
) -> None:
    """Generate 100% stacked bar chart showing classification proportions.

    Creates a horizontal bar chart with:
    - Each segment representing a classification category
    - Labels below showing classification name, contig count (%), and total Mbp (%)
    - Optional interactive HTML version with tooltips

    Args:
        summary_tsv: Path to contig_summary.tsv
        outprefix: Output prefix for generated files
        plot_title_suffix: Suffix for plot titles
        plot_html: Whether to generate HTML version
        assembly_name: Optional assembly name for subtitle
        reference_name: Optional reference name for subtitle
    """
    if not have_rscript():
        logger.warning("Rscript not found; skipping classification summary bar.")
        return

    if not summary_tsv.exists():
        return

    plot_pdf = Path(str(outprefix) + ".classification_summary_bar.pdf")
    plot_html_path = Path(str(outprefix) + ".classification_summary_bar.html")
    r_script_path = Path(str(outprefix) + ".classification_summary_bar.R")

    tmpl_path = Path(__file__).resolve().with_name("classification_summary_bar.tmpl.R")
    if not tmpl_path.exists():
        logger.warning(f"Missing R template: {tmpl_path}")
        return

    tmpl = _read_text(tmpl_path)
    filled = (
        tmpl.replace("__SUMMARY__", _esc(summary_tsv))
        .replace("__OUTPDF__", _esc(plot_pdf))
        .replace("__OUTHTML__", _esc(plot_html_path))
        .replace("__PLOTHTML__", "TRUE" if plot_html else "FALSE")
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
        .replace("__ASMNAME__", str(assembly_name).replace('"', '\\"'))
        .replace("__REFNAME__", str(reference_name).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    logger.info(f"Running Rscript to generate classification summary bar: {plot_pdf}")
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; classification summary bar not generated.")
    else:
        logger.done(f"Classification summary bar written to: {plot_pdf}")
        if plot_html and plot_html_path.exists():
            logger.done(f"Classification summary bar HTML written to: {plot_html_path}")


