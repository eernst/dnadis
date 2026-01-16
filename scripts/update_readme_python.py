#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path


def update_readme(readme_path: Path, python_version: str) -> bool:
    text = readme_path.read_text(encoding="utf-8")
    pattern = re.compile(r"(conda create -n final_finalizer python=)(\d+\.\d+)")
    replaced, n = pattern.subn(rf"\g<1>{python_version}", text)
    if n == 0:
        raise ValueError("Unable to find conda create line with python= in README.md")
    if replaced != text:
        readme_path.write_text(replaced, encoding="utf-8")
        return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description="Update README conda python version.")
    parser.add_argument("--readme", type=Path, default=Path("README.md"))
    parser.add_argument("--python", required=True, help="Python version (e.g. 3.11)")
    args = parser.parse_args()

    changed = update_readme(args.readme, args.python)
    if changed:
        print(f"Updated {args.readme} to python={args.python}")
    else:
        print(f"No change needed in {args.readme}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
