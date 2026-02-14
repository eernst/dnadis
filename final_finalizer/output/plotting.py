#!/usr/bin/env python3
"""Plotting functions using R/ggplot2."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

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
        logger.error(
            "--plot requires Rscript in PATH. "
            "Install with: conda install -c conda-forge r-base"
        )
        return False

    if not _have_rmarkdown():
        logger.error(
            "--plot requires rmarkdown and pandoc. "
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


