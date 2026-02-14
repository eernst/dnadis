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
    rdna_arrays_tsv: Optional[Path] = None,
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


def _have_rmarkdown() -> bool:
    """Check if rmarkdown and pandoc are available for rendering .Rmd files."""
    try:
        result = subprocess.run(
            ["Rscript", "-e", "if (!requireNamespace('rmarkdown', quietly=TRUE)) quit(status=1)"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        if result.returncode != 0:
            return False
        result2 = subprocess.run(
            ["Rscript", "-e", "if (!rmarkdown::pandoc_available()) quit(status=1)"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return result2.returncode == 0
    except Exception:
        return False


def run_unified_report(
    summary_tsv: Path,
    ref_lengths_tsv: Path,
    segments_tsv: Path,
    chain_summary_tsv: Path,
    macro_blocks_tsv: Path,
    outprefix: Path,
    chr_like_minlen: int,
    plot_title_suffix: str,
    synteny_mode: str,
    assembly_name: str = "",
    reference_name: str = "",
    rdna_annotations_tsv: Optional[Path] = None,
    rdna_arrays_tsv: Optional[Path] = None,
    contaminants_tsv: Optional[Path] = None,
    agp_tsv: Optional[Path] = None,
    top_n_contaminants: int = 10,
) -> bool:
    """Generate unified HTML report with all plots and summary tables.

    All plots (chromosome overview, classification bar, depth overview,
    contaminant table) are built natively within the Rmd template — no
    external .rds files needed.  Individual PDFs are also exported from
    within the Rmd via ggsave().

    Returns True if the report was generated successfully.
    """
    if not have_rscript():
        logger.warning("Rscript not found; skipping unified report.")
        return False

    if not _have_rmarkdown():
        logger.warning(
            "rmarkdown or pandoc not available; skipping unified report. "
            "Install with: conda install -c conda-forge r-rmarkdown pandoc"
        )
        return False

    tmpl_path = Path(__file__).resolve().with_name("unified_report.tmpl.Rmd")
    if not tmpl_path.exists():
        logger.warning(f"Missing Rmd template: {tmpl_path}")
        return False

    report_rmd = Path(str(outprefix) + ".unified_report.Rmd").resolve()
    report_html = Path(str(outprefix) + ".unified_report.html").resolve()

    tmpl = _read_text(tmpl_path)

    # rmarkdown::render() changes cwd to the Rmd directory, so all paths
    # must be absolute to avoid "file not found" errors.
    def abs_esc(p: Path) -> str:
        return _esc(p.resolve())

    # Optional TSV paths (absolute, empty string if missing)
    rdna_str = abs_esc(rdna_annotations_tsv) if rdna_annotations_tsv and rdna_annotations_tsv.exists() else ""
    rdna_arrays_str = abs_esc(rdna_arrays_tsv) if rdna_arrays_tsv and rdna_arrays_tsv.exists() else ""
    contam_str = abs_esc(contaminants_tsv) if contaminants_tsv and contaminants_tsv.exists() else ""
    agp_str = abs_esc(agp_tsv) if agp_tsv and agp_tsv.exists() else ""

    filled = (
        tmpl.replace("__SUMMARY__", abs_esc(summary_tsv))
        .replace("__REF__", abs_esc(ref_lengths_tsv))
        .replace("__SEGMENTS__", abs_esc(segments_tsv))
        .replace("__EVIDENCE__", abs_esc(chain_summary_tsv))
        .replace("__MACRO__", abs_esc(macro_blocks_tsv))
        .replace("__RDNA_ANNOTATIONS__", rdna_str)
        .replace("__RDNA_ARRAYS__", rdna_arrays_str)
        .replace("__CONTAMINANTS_TSV__", contam_str)
        .replace("__AGP__", agp_str)
        .replace("__ASMNAME__", str(assembly_name).replace('"', '\\"'))
        .replace("__REFNAME__", str(reference_name).replace('"', '\\"'))
        .replace("__SYNTENY_MODE__", str(synteny_mode))
        .replace("__CHRLIKE__", str(int(chr_like_minlen)))
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
        .replace("__OUTPREFIX__", abs_esc(outprefix))
        .replace("__TOP_N__", str(int(top_n_contaminants)))
    )

    with report_rmd.open("w", encoding="utf-8") as fh:
        fh.write(filled)

    logger.info(f"Rendering unified report: {report_html}")
    try:
        subprocess.run(
            ["Rscript", "-e", f"rmarkdown::render('{_esc(report_rmd)}')"],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; unified report not generated.")
        return False
    else:
        if report_html.exists():
            logger.done(f"Unified report written to: {report_html}")
            return True
        return False


