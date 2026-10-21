#!/usr/bin/env python3
"""
Bundle the HTML reports from a dnadis output folder into a compressed tarball.

Point this at a dnadis output directory and it collects every ``*.html`` report
(the top-level comparison report plus each per-assembly report) into a single
``.tar.gz``. The per-assembly subdirectory layout is preserved, so the relative
links from the comparison report to the individual assembly reports keep
working when the archive is extracted.

By default only HTML files are included; all other outputs (FASTA, TSV, PDF,
logs, etc.) are left out to keep the archive small and shareable. With
``--tsv``, the TSV tables dnadis writes are added as well: those at the top
level, and in each subdirectory the ones named after it (``NAME/NAME.*.tsv``,
e.g. per-assembly tables and ``reference/reference.*.tsv``). Intermediate TSVs
in deeper folders (pairwise caches, classification work directories) are not.

Usage:
    python scripts/bundle_reports.py OUTPUT_DIR
    python scripts/bundle_reports.py OUTPUT_DIR -o reports.tar.gz
    python scripts/bundle_reports.py OUTPUT_DIR --tsv
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path


def find_html_reports(output_dir: Path) -> list[Path]:
    """Return all HTML files under output_dir, sorted for stable archives."""
    return sorted(output_dir.rglob("*.html"))


def find_tsv_outputs(output_dir: Path) -> list[Path]:
    """Return dnadis TSV tables: top-level ones and NAME/NAME.*.tsv in each subdirectory."""
    tsvs = list(output_dir.glob("*.tsv"))
    for sub in output_dir.iterdir():
        if sub.is_dir():
            tsvs += sub.glob(f"{sub.name}.tsv")
            tsvs += sub.glob(f"{sub.name}.*.tsv")
    return sorted(set(tsvs))


def bundle_reports(output_dir: Path, archive_path: Path, include_tsv: bool = False) -> list[Path]:
    """Write a tarball of the HTML reports (and optionally TSV tables), preserving
    paths relative to output_dir.

    Returns the list of files added to the archive.
    """
    files = find_html_reports(output_dir)
    if include_tsv:
        files += find_tsv_outputs(output_dir)
    if not files:
        kinds = ".html reports or .tsv tables" if include_tsv else ".html reports"
        raise FileNotFoundError(f"No {kinds} found under {output_dir}")

    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "w:gz") as tar:
        for path in files:
            # arcname is relative to the output folder, so the comparison
            # report's links into the per-assembly subdirs still resolve.
            tar.add(path, arcname=str(path.relative_to(output_dir)))
    return files


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
        help="Path for the tarball (default: <output_dir>.reports.tar.gz, or "
        "<output_dir>.reports_tsv.tar.gz with --tsv, alongside the output folder)",
    )
    parser.add_argument(
        "--tsv",
        action="store_true",
        help="Also include the TSV tables: top-level ones and NAME/NAME.*.tsv in "
        "each subdirectory (per-assembly tables, reference/reference.*.tsv)",
    )
    args = parser.parse_args()

    output_dir: Path = args.output_dir
    if not output_dir.is_dir():
        print(f"[error] not a directory: {output_dir}", file=sys.stderr)
        return 1

    output_dir = output_dir.resolve()
    archive_path: Path = args.output
    if archive_path is None:
        suffix = "reports_tsv" if args.tsv else "reports"
        archive_path = output_dir.parent / f"{output_dir.name}.{suffix}.tar.gz"
    archive_path = archive_path.resolve()

    try:
        files = bundle_reports(output_dir, archive_path, include_tsv=args.tsv)
    except FileNotFoundError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    for path in files:
        print(f"  + {path.relative_to(output_dir)}")
    n_html = sum(p.suffix == ".html" for p in files)
    tsv_note = f" and {len(files) - n_html} TSV table(s)" if args.tsv else ""
    size_mb = archive_path.stat().st_size / (1024 * 1024)
    print(f"[done] wrote {n_html} report(s){tsv_note} to {archive_path} ({size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
