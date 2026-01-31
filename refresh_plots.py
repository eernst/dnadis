#!/usr/bin/env python3
"""Re-render R plot scripts from current templates and optionally re-run them.

Usage:
    ./refresh_plots.py /path/to/run_folder              # refresh + run all
    ./refresh_plots.py /path/to/run_folder --dry-run     # show what would change
    ./refresh_plots.py /path/to/run_folder --no-run      # refresh scripts only
    ./refresh_plots.py /path/to/run_folder --only bar    # refresh only matching scripts
    ./refresh_plots.py /path/to/run_folder --only bar contigs  # multiple filters
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Map generated R script suffix → template filename
SCRIPT_TO_TEMPLATE = {
    "chromosome_overview": "chromosome_overview.tmpl.R",
    "depth_overview": "depth_overview.tmpl.R",
    "contaminant_treemap": "contaminant_treemap.tmpl.R",
    "contaminant_table": "contaminant_table.tmpl.R",
    "classification_summary_bar": "classification_summary_bar.tmpl.R",
    "classification_summary_contigs": "classification_summary_contigs.tmpl.R",
}

# Regex to find all __PLACEHOLDER__ tokens in a template
PLACEHOLDER_RE = re.compile(r"__([A-Z_]+)__")


def find_generated_scripts(run_dir: Path) -> list[tuple[str, Path]]:
    """Find generated .R scripts in run_dir and return (suffix, path) pairs."""
    found = []
    for r_file in sorted(run_dir.glob("*.R")):
        for suffix in SCRIPT_TO_TEMPLATE:
            if r_file.name.endswith(f".{suffix}.R"):
                found.append((suffix, r_file))
                break
    return found


def extract_placeholder_values(
    template_text: str, generated_text: str
) -> dict[str, str]:
    """Extract placeholder values by matching context around each placeholder.

    For each __FOO__ in the template, takes the surrounding text (same line,
    before and after the placeholder) and searches for that context in the
    generated script to recover the substituted value.
    """
    values: dict[str, str] = {}

    for match in PLACEHOLDER_RE.finditer(template_text):
        placeholder = match.group(0)
        if placeholder in values:
            continue  # already extracted from first occurrence

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
        # Use non-greedy match for the value
        pattern = re.escape(context_before) + r"(.+?)" + re.escape(context_after)
        m = re.search(pattern, generated_text)
        if m:
            values[placeholder] = m.group(1)
        else:
            # Fallback: try without trailing context (handles trailing whitespace diffs)
            pattern_no_trail = re.escape(context_before) + r"(.+)"
            m2 = re.search(pattern_no_trail, generated_text)
            if m2:
                values[placeholder] = m2.group(1).rstrip()

    return values


def refresh_script(
    suffix: str,
    generated_path: Path,
    template_dir: Path,
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

    # Extract placeholder values from the existing generated script
    values = extract_placeholder_values(template_text, generated_text)

    # Check we found all placeholders
    all_placeholders = set(PLACEHOLDER_RE.findall(template_text))
    found_keys = {k.strip("_") for k in values}  # strip __ wrappers
    missing = {f"__{p}__" for p in all_placeholders} - set(values.keys())
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
        print(f"  {suffix}: would update ({len(values)} placeholders extracted)")
        return True

    generated_path.write_text(filled, encoding="utf-8")
    print(f"  {suffix}: updated ({len(values)} placeholders)")
    return True


def run_script(script_path: Path) -> bool:
    """Run a generated R script with Rscript."""
    print(f"  Running: Rscript {script_path.name} ... ", end="", flush=True)
    result = subprocess.run(
        ["Rscript", str(script_path)],
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
        default=Path(__file__).resolve().parent / "final_finalizer" / "output",
        help="Path to directory containing .tmpl.R files (default: auto-detect from script location)",
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

    # Filter by --only patterns
    if args.only:
        scripts = [
            (suffix, path) for suffix, path in scripts
            if any(pat in suffix for pat in args.only)
        ]
        if not scripts:
            sys.exit(f"No scripts matched --only {' '.join(args.only)}")

    print(f"Run folder: {run_dir}")
    print(f"Templates:  {template_dir}")
    print(f"Found {len(scripts)} script(s):")
    print()

    updated = []
    for suffix, script_path in scripts:
        changed = refresh_script(suffix, script_path, template_dir, dry_run=args.dry_run)
        if changed:
            updated.append((suffix, script_path))

    if args.dry_run:
        print(f"\nDry run: {len(updated)} script(s) would be updated.")
        return

    if not updated and not args.no_run:
        # Even if nothing changed, user might want to re-run
        print("\nNo scripts changed. Use Rscript directly to re-run existing scripts.")
        return

    if updated and not args.no_run:
        print(f"\nRe-running {len(updated)} updated script(s):")
        for suffix, script_path in updated:
            run_script(script_path)


if __name__ == "__main__":
    main()
