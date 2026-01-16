#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def parse_conda_list(text: str) -> dict[str, str]:
    versions: dict[str, str] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line.strip())
        if len(parts) < 2:
            continue
        name, version = parts[0], parts[1]
        versions[name] = version
    return versions


def render_versions(packages: list[str], versions: dict[str, str]) -> str:
    lines = []
    for pkg in packages:
        ver = versions.get(pkg, "unknown")
        lines.append(f"- {pkg}: {ver}")
    return "\n".join(lines)


def update_readme(readme_path: Path, packages: list[str], conda_list_text: str) -> bool:
    text = readme_path.read_text(encoding="utf-8")
    versions = parse_conda_list(conda_list_text)
    block = render_versions(packages, versions)

    pattern = re.compile(
        r"(<!-- conda-versions-start -->)(.*?)(<!-- conda-versions-end -->)",
        re.DOTALL,
    )
    match = pattern.search(text)
    if not match:
        raise ValueError("Unable to find conda versions markers in README.md")
    replacement = f"{match.group(1)}\n{block}\n{match.group(3)}"
    replaced = pattern.sub(replacement, text, count=1)
    if replaced != text:
        readme_path.write_text(replaced, encoding="utf-8")
        return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description="Update README with conda package versions.")
    parser.add_argument("--readme", type=Path, default=Path("README.md"))
    parser.add_argument("--packages", nargs="+", required=True)
    parser.add_argument("--conda-list", type=Path, required=True)
    args = parser.parse_args()

    conda_list_text = args.conda_list.read_text(encoding="utf-8")
    changed = update_readme(args.readme, args.packages, conda_list_text)
    if changed:
        print(f"Updated {args.readme} conda versions block")
    else:
        print(f"No change needed in {args.readme}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
