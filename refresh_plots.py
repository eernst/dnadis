#!/usr/bin/env python3
"""Re-render R plot scripts from current templates and optionally re-run them.

Usage:
    ./refresh_plots.py /path/to/run_folder              # refresh + run changed
    ./refresh_plots.py /path/to/run_folder --force      # refresh + run all
    ./refresh_plots.py /path/to/run_folder --dry-run    # show what would change
    ./refresh_plots.py /path/to/run_folder --no-run     # refresh scripts only
    ./refresh_plots.py /path/to/run_folder --only bar   # only matching scripts
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Map generated R script suffix → template filename
SCRIPT_TO_TEMPLATE = {
    "assembly_report": "assembly_report.tmpl.Rmd",
    "comparison_report": "comparison_report.tmpl.Rmd",
    # Legacy name for backward compat with existing output dirs
    "unified_report": "assembly_report.tmpl.Rmd",
}

# Regex to find all __PLACEHOLDER__ tokens in a template
PLACEHOLDER_RE = re.compile(r"__([A-Z_]+)__")

# Which TSV file suffix provides the value for each placeholder
# (relative to the output prefix)
PLACEHOLDER_TSV = {
    "__SUMMARY__": ".contig_summary.tsv",
    "__REF__": ".ref_lengths.tsv",
    "__SEGMENTS__": ".segments.tsv",
    "__EVIDENCE__": ".evidence_summary.tsv",
    "__MACRO__": ".macro_blocks.tsv",
    "__RDNA_ANNOTATIONS__": ".rdna_annotations.tsv",
    "__AGP__": ".scaffolded.agp",
    "__CONTAMINANTS_TSV__": ".contaminants.tsv",
}

# Output file suffixes per script type
SCRIPT_OUTPUT_SUFFIXES = {
    "assembly_report": (None, ".assembly_report.html"),
    "comparison_report": (None, ".comparison_report.html"),
    # Legacy name for backward compat with existing output dirs
    "unified_report": (None, ".unified_report.html"),
}


def find_generated_scripts(run_dir: Path) -> list[tuple[str, Path]]:
    """Find generated .R and .Rmd scripts in run_dir and return (suffix, path) pairs."""
    found = []
    for r_file in sorted(list(run_dir.glob("*.R")) + list(run_dir.glob("*.Rmd"))):
        for suffix, tmpl_name in SCRIPT_TO_TEMPLATE.items():
            ext = ".Rmd" if tmpl_name.endswith(".Rmd") else ".R"
            if r_file.name.endswith(f".{suffix}{ext}"):
                found.append((suffix, r_file))
                break
    return found


def detect_prefix(run_dir: Path, scripts: list[tuple[str, Path]]) -> Path:
    """Detect the output prefix from a generated R script filename.

    E.g. if the script is 'la8230_vs_diploids.chromosome_overview.R'
    the prefix is run_dir / 'la8230_vs_diploids'.
    """
    for suffix, script_path in scripts:
        stem = script_path.name
        tmpl_name = SCRIPT_TO_TEMPLATE[suffix]
        ext = ".Rmd" if tmpl_name.endswith(".Rmd") else ".R"
        expected_end = f".{suffix}{ext}"
        if stem.endswith(expected_end):
            prefix_name = stem[: -len(expected_end)]
            return run_dir / prefix_name
    sys.exit("Error: could not detect output prefix from R script names")


def extract_non_path_placeholders(
    template_text: str, generated_text: str
) -> dict[str, str]:
    """Extract non-path placeholder values (PLOTHTML, CHRLIKE, SUFFIX, TOP_N)
    by matching context around each placeholder in the generated script.
    """
    # Only extract placeholders that aren't TSV paths or output paths
    path_placeholders = set(PLACEHOLDER_TSV.keys()) | {
        "__OUTPDF__", "__OUTHTML__",
    }

    values: dict[str, str] = {}
    for match in PLACEHOLDER_RE.finditer(template_text):
        placeholder = match.group(0)
        if placeholder in values or placeholder in path_placeholders:
            continue

        start, end = match.start(), match.end()

        # Get context on the same line before the placeholder
        line_start = template_text.rfind("\n", 0, start) + 1
        context_before = template_text[line_start:start]

        # Get context on the same line after the placeholder
        line_end = template_text.find("\n", end)
        if line_end == -1:
            line_end = len(template_text)
        context_after = template_text[end:line_end]

        # Build regex: literal context before + capture group + literal context after
        pattern = re.escape(context_before) + r"(.*?)" + re.escape(context_after)
        m = re.search(pattern, generated_text)
        if m:
            values[placeholder] = m.group(1)
        else:
            # Fallback: try without trailing context
            pattern_no_trail = re.escape(context_before) + r"(.*)"
            m2 = re.search(pattern_no_trail, generated_text)
            if m2:
                values[placeholder] = m2.group(1).rstrip()

    return values


def esc(s: Path | str) -> str:
    """Escape path for R (forward slashes)."""
    return str(s).replace("\\", "/")


def build_placeholder_values(
    suffix: str,
    prefix: Path,
    template_text: str,
    generated_text: str,
    template_dir: Path | None = None,
) -> dict[str, str]:
    """Build placeholder → value mapping from TSVs in run_dir and the generated script."""
    values: dict[str, str] = {}

    # Path-based placeholders: resolve from prefix + known TSV suffixes
    for placeholder, tsv_suffix in PLACEHOLDER_TSV.items():
        if placeholder in template_text:
            tsv_path = Path(str(prefix) + tsv_suffix)
            # For optional TSVs, use empty string if file doesn't exist
            # (R templates handle missing files gracefully)
            if tsv_path.exists():
                values[placeholder] = esc(tsv_path)
            else:
                values[placeholder] = ""

    # Output file placeholders
    pdf_suffix, html_suffix = SCRIPT_OUTPUT_SUFFIXES[suffix]
    if "__OUTPDF__" in template_text and pdf_suffix:
        values["__OUTPDF__"] = esc(str(prefix) + pdf_suffix)
    if "__OUTHTML__" in template_text and html_suffix:
        values["__OUTHTML__"] = esc(str(prefix) + html_suffix)

    # Output prefix placeholder
    if "__OUTPREFIX__" in template_text:
        values["__OUTPREFIX__"] = esc(prefix)

    # Shared R and CSS files (resolved relative to template directory)
    if "__COMMON_R__" in template_text and template_dir:
        values["__COMMON_R__"] = esc(template_dir / "common.R")
    if "__COMMON_CSS__" in template_text and template_dir:
        values["__COMMON_CSS__"] = esc(template_dir / "common.css")

    # Non-path placeholders: extract from generated script (PLOTHTML, CHRLIKE, SUFFIX, TOP_N)
    non_path = extract_non_path_placeholders(template_text, generated_text)
    values.update(non_path)

    # Cross-reference shared placeholders from sibling scripts in the same run
    # folder.  This is necessary because some scripts (e.g., the unified report)
    # may have been generated with stale/default values during development, while
    # the standalone scripts were generated with the correct pipeline values.
    missing = {f"__{p}__" for p in PLACEHOLDER_RE.findall(template_text)} - set(values.keys())
    needs_crossref = missing | {
        p for p in values
        if values[p].startswith("__") or (p == "__CHRLIKE__" and values.get(p, "") in ("", "50000"))
    }
    if needs_crossref:
        # Try extracting from the chromosome_overview.R script first (most
        # likely to have authoritative CHRLIKE, SUFFIX values)
        for sibling_suffix in ("chromosome_overview", "classification_summary_bar", "depth_overview"):
            if sibling_suffix == suffix:
                continue
            tmpl_name = SCRIPT_TO_TEMPLATE.get(sibling_suffix)
            if not tmpl_name:
                continue
            ext = ".Rmd" if tmpl_name.endswith(".Rmd") else ".R"
            sibling_path = Path(str(prefix) + f".{sibling_suffix}{ext}")
            if not sibling_path.exists():
                continue
            sibling_tmpl = (template_dir / tmpl_name) if template_dir else None
            if sibling_tmpl and sibling_tmpl.exists():
                sibling_tmpl_text = sibling_tmpl.read_text(encoding="utf-8")
                sibling_gen_text = sibling_path.read_text(encoding="utf-8")
                sibling_vals = extract_non_path_placeholders(sibling_tmpl_text, sibling_gen_text)
                for p in list(needs_crossref):
                    if p in sibling_vals and sibling_vals[p] and not sibling_vals[p].startswith("__"):
                        values[p] = sibling_vals[p]
                        needs_crossref.discard(p)
            if not needs_crossref:
                break

    # Defaults for placeholders that couldn't be extracted or were extracted
    # as literal placeholder strings (from a previously broken generated script)
    defaults = {
        "__PLOTHTML__": "TRUE",
        "__TOP_N__": "10",
    }
    for placeholder, default in defaults.items():
        if placeholder in template_text:
            current = values.get(placeholder, "")
            # Use default if missing or if the "value" is still a placeholder
            if not current or current.startswith("__"):
                values[placeholder] = default

    return values


def refresh_script(
    suffix: str,
    generated_path: Path,
    template_dir: Path,
    prefix: Path,
    dry_run: bool = False,
) -> bool:
    """Re-render a single generated R script from the current template.

    Returns True if the script was updated (or would be updated in dry-run mode).
    """
    template_name = SCRIPT_TO_TEMPLATE[suffix]
    template_path = template_dir / template_name

    if not template_path.exists():
        print(f"  SKIP {suffix}: template not found at {template_path}", file=sys.stderr)
        return False

    template_text = template_path.read_text(encoding="utf-8")
    generated_text = generated_path.read_text(encoding="utf-8")

    # Build placeholder values from TSVs in run_dir + generated script
    values = build_placeholder_values(suffix, prefix, template_text, generated_text,
                                       template_dir=template_dir)

    # Check we found all placeholders
    all_placeholders = {f"__{p}__" for p in PLACEHOLDER_RE.findall(template_text)}
    missing = all_placeholders - set(values.keys())
    if missing:
        print(f"  WARN {suffix}: could not extract values for: {', '.join(sorted(missing))}", file=sys.stderr)

    # Apply values to fresh template
    filled = template_text
    for placeholder, value in values.items():
        filled = filled.replace(placeholder, value)

    # Check if anything changed
    if filled == generated_text:
        print(f"  {suffix}: unchanged")
        return False

    if dry_run:
        print(f"  {suffix}: would update ({len(values)} placeholders)")
        return True

    generated_path.write_text(filled, encoding="utf-8")
    print(f"  {suffix}: updated ({len(values)} placeholders)")
    return True


def run_script(script_path: Path) -> bool:
    """Run a generated R script with Rscript, or render an Rmd with rmarkdown."""
    is_rmd = script_path.suffix.lower() == ".rmd"
    if is_rmd:
        esc_path = str(script_path).replace("\\", "/")
        cmd = ["Rscript", "-e", f"rmarkdown::render('{esc_path}')"]
        label = f"rmarkdown::render({script_path.name})"
    else:
        cmd = ["Rscript", str(script_path)]
        label = f"Rscript {script_path.name}"

    print(f"  Running: {label} ... ", end="", flush=True)
    result = subprocess.run(
        cmd,
        cwd=str(script_path.parent),
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        print("ok")
        return True
    else:
        print(f"FAILED (exit {result.returncode})")
        # Show last few lines of stderr for debugging
        stderr_lines = result.stderr.strip().splitlines()
        for line in stderr_lines[-10:]:
            print(f"    {line}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Refresh R plot scripts from current templates and optionally re-run them.",
    )
    parser.add_argument(
        "run_dir",
        type=Path,
        help="Path to a final_finalizer output/run folder",
    )
    parser.add_argument(
        "--template-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "final_finalizer" / "output" / "reports",
        help="Path to directory containing .tmpl.Rmd files (default: auto-detect from script location)",
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Show what would change without writing files",
    )
    parser.add_argument(
        "--no-run",
        action="store_true",
        help="Refresh R scripts but don't execute them",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        metavar="PATTERN",
        help="Only refresh scripts whose suffix contains one of these substrings "
             "(e.g., 'bar' matches classification_summary_bar)",
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Re-run scripts even if unchanged",
    )
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        sys.exit(f"Error: not a directory: {run_dir}")

    template_dir = args.template_dir.resolve()
    if not template_dir.is_dir():
        sys.exit(f"Error: template directory not found: {template_dir}")

    scripts = find_generated_scripts(run_dir)
    if not scripts:
        sys.exit(f"No generated .R scripts found in {run_dir}")

    # Detect output prefix from script filenames
    prefix = detect_prefix(run_dir, scripts)
    print(f"Run folder: {run_dir}")
    print(f"Prefix:     {prefix}")
    print(f"Templates:  {template_dir}")

    # Filter by --only patterns
    if args.only:
        scripts = [
            (suffix, path) for suffix, path in scripts
            if any(pat in suffix for pat in args.only)
        ]
        if not scripts:
            sys.exit(f"No scripts matched --only {' '.join(args.only)}")

    print(f"Found {len(scripts)} script(s):")
    print()

    updated = []
    for suffix, script_path in scripts:
        changed = refresh_script(suffix, script_path, template_dir, prefix, dry_run=args.dry_run)
        if changed:
            updated.append((suffix, script_path))

    if args.dry_run:
        print(f"\nDry run: {len(updated)} script(s) would be updated.")
        return

    if args.no_run:
        return

    # Determine which scripts to run
    if args.force:
        to_run = scripts
    elif updated:
        to_run = updated
    else:
        print("\nNo scripts changed. Use --force to re-run anyway.")
        return

    print(f"\nRunning {len(to_run)} script(s):")
    for suffix, script_path in to_run:
            run_script(script_path)


if __name__ == "__main__":
    main()
