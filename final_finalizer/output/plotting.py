#!/usr/bin/env python3
"""Plotting functions using R/ggplot2."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional


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
        print("[warn] --plot specified, but Rscript not found in PATH; skipping plot.", file=sys.stderr)
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

    print(f"[info] Running Rscript to generate plot: {plot_pdf}", file=sys.stderr)
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[warn] Rscript failed with code {e.returncode}; plot not generated.", file=sys.stderr)
    else:
        print(f"[done] Plot written to: {plot_pdf}", file=sys.stderr)
        if plot_html:
            print(f"[done] Plot written to: {plot_html_path}", file=sys.stderr)


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
        print("[warn] --plot specified, but Rscript not found in PATH; skipping depth plot.", file=sys.stderr)
        return

    # Check if summary file has depth data
    try:
        with summary_tsv.open("r") as fh:
            header = fh.readline().strip().split("\t")
            if "depth_mean" not in header:
                print("[info] No depth data in summary; skipping depth plot.", file=sys.stderr)
                return
    except Exception as e:
        print(f"[warn] Could not check summary file for depth data: {e}", file=sys.stderr)
        return

    plot_pdf = Path(str(outprefix) + ".depth_overview.pdf")
    plot_html_path = Path(str(outprefix) + ".depth_overview.html")
    r_script_path = Path(str(outprefix) + ".depth_overview.R")

    # Template expected alongside this Python file
    tmpl_path = Path(__file__).resolve().with_name("depth_overview.tmpl.R")
    if not tmpl_path.exists():
        print(f"[warn] Missing R template for depth plots: {tmpl_path}", file=sys.stderr)
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

    print(f"[info] Running Rscript to generate depth plot: {plot_pdf}", file=sys.stderr)
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[warn] Rscript failed with code {e.returncode}; depth plot not generated.", file=sys.stderr)
    else:
        print(f"[done] Depth plot written to: {plot_pdf}", file=sys.stderr)
        if plot_html:
            print(f"[done] Depth plot written to: {plot_html_path}", file=sys.stderr)
