#!/usr/bin/env python3
"""
Bundle the HTML reports from a dnadis output folder into a compressed tarball.

Point this at a dnadis output directory and it collects every ``*.html`` report
(the top-level comparison report plus each per-assembly report) into a single
``.tar.gz``. The per-assembly subdirectory layout is preserved, so the relative
links from the comparison report to the individual assembly reports keep
working when the archive is extracted.

Only HTML files are included; all other outputs (FASTA, TSV, PDF, logs, etc.)
are left out to keep the archive small and shareable.

Usage:
    python scripts/bundle_reports.py OUTPUT_DIR
    python scripts/bundle_reports.py OUTPUT_DIR -o reports.tar.gz
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path


def find_html_reports(output_dir: Path) -> list[Path]:
    """Return all HTML files under output_dir, sorted for stable archives."""
    return sorted(output_dir.rglob("*.html"))


def bundle_reports(output_dir: Path, archive_path: Path) -> list[Path]:
    """Write an HTML-only tarball, preserving paths relative to output_dir.

    Returns the list of files added to the archive.
    """
    html_files = find_html_reports(output_dir)
    if not html_files:
        raise FileNotFoundError(f"No .html reports found under {output_dir}")

    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "w:gz") as tar:
        for html in html_files:
            # arcname is relative to the output folder, so the comparison
            # report's links into the per-assembly subdirs still resolve.
            tar.add(html, arcname=str(html.relative_to(output_dir)))
    return html_files


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Bundle the HTML reports from a dnadis output folder into a "
        "compressed tarball, preserving per-assembly subdirectories.",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="dnadis output directory containing the HTML reports",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Path for the tarball (default: <output_dir>.reports.tar.gz "
        "alongside the output folder)",
    )
    args = parser.parse_args()

    output_dir: Path = args.output_dir
    if not output_dir.is_dir():
        print(f"[error] not a directory: {output_dir}", file=sys.stderr)
        return 1

    output_dir = output_dir.resolve()
    archive_path: Path = args.output
    if archive_path is None:
        archive_path = output_dir.parent / f"{output_dir.name}.reports.tar.gz"
    archive_path = archive_path.resolve()

    try:
        html_files = bundle_reports(output_dir, archive_path)
    except FileNotFoundError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    for html in html_files:
        print(f"  + {html.relative_to(output_dir)}")
    size_mb = archive_path.stat().st_size / (1024 * 1024)
    print(f"[done] wrote {len(html_files)} report(s) to {archive_path} ({size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
