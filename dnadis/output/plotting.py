#!/usr/bin/env python3
"""Plotting functions using R/ggplot2."""
from __future__ import annotations

import functools
import subprocess
from pathlib import Path
from typing import List, Optional

from dnadis.utils.logging_config import get_logger

logger = get_logger("plotting")

# --- Module-level constants ---
_REPORTS_DIR = Path(__file__).resolve().parent / "reports"
_COMMON_R = _REPORTS_DIR / "common.R"
_COMMON_CSS = _REPORTS_DIR / "common.css"


@functools.lru_cache(maxsize=None)
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


def _abs_esc(p: Path) -> str:
    """Resolve path to absolute and escape for R."""
    return _esc(p.resolve())


@functools.lru_cache(maxsize=None)
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


def _check_r_dependencies() -> bool:
    """Return True if Rscript and rmarkdown+pandoc are available."""
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
    return True


def run_assembly_report(
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
    compleasm_chrs_summary: Optional[Path] = None,
    compleasm_non_chrs_summary: Optional[Path] = None,
    top_n_contaminants: int = 10,
    rearrangements_tsv: Optional[Path] = None,
    self_contained: bool = False,
) -> bool:
    """Generate unified HTML report with all plots and summary tables.

    All plots (chromosome overview, classification bar, depth overview,
    contaminant table) are built natively within the Rmd template — no
    external .rds files needed.  Individual PDFs are also exported from
    within the Rmd via ggsave().

    Returns True if the report was generated successfully.
    """
    if not _check_r_dependencies():
        return False

    tmpl_path = _REPORTS_DIR / "assembly_report.tmpl.Rmd"
    if not tmpl_path.exists():
        logger.warning(f"Missing Rmd template: {tmpl_path}")
        return False

    report_rmd = Path(str(outprefix) + ".assembly_report.Rmd").resolve()
    report_html = Path(str(outprefix) + ".assembly_report.html").resolve()

    tmpl = _read_text(tmpl_path)

    # Optional TSV paths (absolute, empty string if missing)
    rdna_str = _abs_esc(rdna_annotations_tsv) if rdna_annotations_tsv and rdna_annotations_tsv.exists() else ""
    rdna_arrays_str = _abs_esc(rdna_arrays_tsv) if rdna_arrays_tsv and rdna_arrays_tsv.exists() else ""
    contam_str = _abs_esc(contaminants_tsv) if contaminants_tsv and contaminants_tsv.exists() else ""
    agp_str = _abs_esc(agp_tsv) if agp_tsv and agp_tsv.exists() else ""
    compleasm_chrs_str = _abs_esc(compleasm_chrs_summary) if compleasm_chrs_summary and compleasm_chrs_summary.exists() else ""
    compleasm_non_str = _abs_esc(compleasm_non_chrs_summary) if compleasm_non_chrs_summary and compleasm_non_chrs_summary.exists() else ""
    rearr_str = _abs_esc(rearrangements_tsv) if rearrangements_tsv and rearrangements_tsv.exists() else ""

    filled = (
        tmpl.replace("__SUMMARY__", _abs_esc(summary_tsv))
        .replace("__REF__", _abs_esc(ref_lengths_tsv))
        .replace("__SEGMENTS__", _abs_esc(segments_tsv))
        .replace("__EVIDENCE__", _abs_esc(chain_summary_tsv))
        .replace("__MACRO__", _abs_esc(macro_blocks_tsv))
        .replace("__RDNA_ANNOTATIONS__", rdna_str)
        .replace("__RDNA_ARRAYS__", rdna_arrays_str)
        .replace("__CONTAMINANTS_TSV__", contam_str)
        .replace("__AGP__", agp_str)
        .replace("__COMPLEASM_CHRS__", compleasm_chrs_str)
        .replace("__COMPLEASM_NONCHRS__", compleasm_non_str)
        .replace("__REARRANGEMENTS__", rearr_str)
        .replace("__ASMNAME__", str(assembly_name).replace('"', '\\"'))
        .replace("__REFNAME__", str(reference_name).replace('"', '\\"'))
        .replace("__SYNTENY_MODE__", str(synteny_mode))
        .replace("__CHRLIKE__", str(int(chr_like_minlen)))
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
        .replace("__OUTPREFIX__", _abs_esc(outprefix))
        .replace("__TOP_N__", str(int(top_n_contaminants)))
        .replace("__COMMON_R__", _esc(_COMMON_R))
        .replace("__COMMON_CSS__", _esc(_COMMON_CSS))
        .replace("__SELF_CONTAINED__", "true" if self_contained else "false")
    )

    with report_rmd.open("w", encoding="utf-8") as fh:
        fh.write(filled)

    logger.info(f"Rendering assembly report: {report_html}")
    try:
        subprocess.run(
            ["Rscript", "-e", f"rmarkdown::render('{_esc(report_rmd)}', quiet = TRUE)"],
            check=True, capture_output=True, text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; assembly report not generated.")
        if e.stderr:
            for line in e.stderr.strip().splitlines()[-20:]:
                logger.warning(f"  R: {line}")
        return False
    else:
        if report_html.exists():
            logger.done(f"Assembly report written to: {report_html}")
            return True
        return False


def run_comparison_report(
    comparison_tsv: Path,
    completeness_tsv: Path,
    ref_lengths_tsv: Path,
    assembly_results: List,
    outprefix: Path,
    chr_like_minlen: int,
    synteny_mode: str,
    reference_name: str = "",
    pairwise_pairs: Optional[List[tuple]] = None,
    self_contained: bool = False,
    assembly_sort_order: str = "input",
) -> bool:
    """Generate cross-assembly comparison HTML report.

    Follows the same pattern as :func:`run_assembly_report`: loads the
    ``comparison_report.tmpl.Rmd`` template, substitutes placeholders,
    and calls ``rmarkdown::render()``.

    Args:
        comparison_tsv: Path to comparison_summary.tsv.
        completeness_tsv: Path to chromosome_completeness.tsv.
        ref_lengths_tsv: Path to reference lengths TSV.
        assembly_results: List of :class:`AssemblyResult` objects.
        outprefix: Output prefix for the comparison report.
        chr_like_minlen: Minimum length threshold for chromosome-like contigs.
        synteny_mode: Synteny mode used ("protein" or "nucleotide").
        reference_name: Reference name for the report header.
        pairwise_pairs: List of (pair_name, tsv_path) tuples for pairwise
            macro_blocks (nucleotide mode only). Each pair_name is
            "{left_asm}_vs_{right_asm}".
        assembly_sort_order: Assembly ordering in comparison report.
            "input" preserves FOFN/directory order; "identity" sorts by
            descending median sequence identity vs reference.

    Returns:
        True if the report was generated successfully.
    """
    if not _check_r_dependencies():
        return False

    tmpl_path = _REPORTS_DIR / "comparison_report.tmpl.Rmd"
    if not tmpl_path.exists():
        logger.warning(f"Missing Rmd template: {tmpl_path}")
        return False

    report_rmd = Path(str(outprefix) + ".comparison_report.Rmd").resolve()
    report_html = Path(str(outprefix) + ".comparison_report.html").resolve()

    tmpl = _read_text(tmpl_path)

    # Build semicolon-separated lists of per-assembly TSV paths and names
    per_asm_tsvs = ";".join(
        _abs_esc(r.summary_tsv) for r in assembly_results
    )
    asm_names = ";".join(r.assembly_name for r in assembly_results)

    # Per-assembly macro_blocks TSVs (for ref→asm ribbons in synteny plot)
    per_asm_macro_tsvs = ";".join(
        _abs_esc(r.macro_blocks_tsv) for r in assembly_results
    )

    # Per-assembly rDNA arrays TSVs (for rDNA comparison section)
    per_asm_rdna_tsvs = ";".join(
        _abs_esc(r.rdna_arrays_tsv) if r.rdna_arrays_tsv and r.rdna_arrays_tsv.exists() else ""
        for r in assembly_results
    )

    # Per-assembly contaminants TSVs (for contamination Top Taxa section)
    per_asm_contam_tsvs = ";".join(
        _abs_esc(r.contaminants_tsv) if r.contaminants_tsv and r.contaminants_tsv.exists() else ""
        for r in assembly_results
    )

    # Per-assembly rearrangement TSVs (for chromosome summary section)
    per_asm_rearr_tsvs = ";".join(
        _abs_esc(r.rearrangements_tsv) if r.rearrangements_tsv and r.rearrangements_tsv.exists() else ""
        for r in assembly_results
    )

    # Per-assembly AGP files (for scaffold join indicators in synteny plot)
    per_asm_agp_tsvs = ";".join(
        _abs_esc(r.agp_tsv) if r.agp_tsv and r.agp_tsv.exists() else ""
        for r in assembly_results
    )

    # Pairwise macro_blocks TSVs and pair labels (for asm→asm ribbons)
    # Names and paths are kept synchronized as (pair_name, tsv_path) tuples
    if pairwise_pairs:
        pw_macro_str = ";".join(_abs_esc(tsv) for _name, tsv in pairwise_pairs)
        pw_names_str = ";".join(name for name, _tsv in pairwise_pairs)
    else:
        pw_macro_str = ""
        pw_names_str = ""

    filled = (
        tmpl.replace("__COMPARISON_TSV__", _abs_esc(comparison_tsv))
        .replace("__COMPLETENESS_TSV__", _abs_esc(completeness_tsv))
        .replace("__REF_LENGTHS_TSV__", _abs_esc(ref_lengths_tsv))
        .replace("__PER_ASM_TSVS__", per_asm_tsvs)
        .replace("__ASM_NAMES__", asm_names)
        .replace("__PER_ASM_MACRO_TSVS__", per_asm_macro_tsvs)
        .replace("__PER_ASM_RDNA_TSVS__", per_asm_rdna_tsvs)
        .replace("__PER_ASM_CONTAM_TSVS__", per_asm_contam_tsvs)
        .replace("__PER_ASM_REARR_TSVS__", per_asm_rearr_tsvs)
        .replace("__PER_ASM_AGP_TSVS__", per_asm_agp_tsvs)
        .replace("__PAIRWISE_MACRO_TSVS__", pw_macro_str)
        .replace("__PAIRWISE_NAMES__", pw_names_str)
        .replace("__REFNAME__", str(reference_name).replace('"', '\\"'))
        .replace("__SYNTENY_MODE__", str(synteny_mode))
        .replace("__CHRLIKE__", str(int(chr_like_minlen)))
        .replace("__OUTPREFIX__", _abs_esc(outprefix))
        .replace("__COMMON_R__", _esc(_COMMON_R))
        .replace("__COMMON_CSS__", _esc(_COMMON_CSS))
        .replace("__SELF_CONTAINED__", "true" if self_contained else "false")
        .replace("__ASM_SORT_ORDER__", str(assembly_sort_order))
    )

    with report_rmd.open("w", encoding="utf-8") as fh:
        fh.write(filled)

    logger.info(f"Rendering comparison report: {report_html}")
    try:
        subprocess.run(
            ["Rscript", "-e", f"rmarkdown::render('{_esc(report_rmd)}', quiet = TRUE)"],
            check=True, capture_output=True, text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.warning(f"Rscript failed with code {e.returncode}; comparison report not generated.")
        if e.stderr:
            for line in e.stderr.strip().splitlines()[-20:]:
                logger.warning(f"  R: {line}")
        return False
    else:
        if report_html.exists():
            logger.done(f"Comparison report written to: {report_html}")
            return True
        return False
