#!/usr/bin/env python3
"""
Multi-assembly input resolution for dnadis.

Provides functions for reading assembly lists from FOFN (file-of-filenames)
TSV files or scanning directories for FASTA files.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Tuple

from dnadis.utils.logging_config import get_logger

logger = get_logger("multi_assembly")

# Type alias: (assembly_path, name, reads_path_or_None)
AssemblySpec = Tuple[Path, str, Optional[Path]]

_FASTA_EXTENSIONS = {
    ".fasta", ".fa", ".fna",
    ".fasta.gz", ".fa.gz", ".fna.gz",
}


def _has_fasta_extension(path: Path) -> bool:
    """Check whether *path* has a recognised FASTA extension."""
    name = path.name.lower()
    for ext in _FASTA_EXTENSIONS:
        if name.endswith(ext):
            return True
    return False


def _strip_fasta_extension(name: str) -> str:
    """Strip FASTA (and optional .gz) extension from a filename."""
    lower = name.lower()
    for ext in sorted(_FASTA_EXTENSIONS, key=len, reverse=True):
        if lower.endswith(ext):
            return name[: len(name) - len(ext)]
    return name


def parse_fofn(fofn_path: Path) -> List[AssemblySpec]:
    """Parse a file-of-filenames TSV listing assemblies.

    A header row is optional.  When present, ``path`` is required and
    ``name`` / ``reads`` are optional, in any column order.  When absent,
    columns are taken positionally: column 0 = path (required), column 1
    = name (optional), column 2 = reads (optional).

    A header is detected when the first non-blank, non-comment line's
    first column is exactly ``path`` (case-insensitive).  Anything else
    is treated as the first data row and the file is read headerless.

    Args:
        fofn_path: Path to the TSV file.

    Returns:
        List of (assembly_path, name, reads_path) tuples.

    Raises:
        ValueError: On duplicate names or missing files.
    """
    fofn_path = Path(fofn_path).resolve()
    fofn_dir = fofn_path.parent

    # Read all data-bearing lines (skip blanks and comments) up front so
    # we can decide whether the first one is a header.
    raw_lines: List[Tuple[int, str]] = []
    with open(fofn_path) as fh:
        for lineno, line in enumerate(fh, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            raw_lines.append((lineno, stripped))

    if not raw_lines:
        raise ValueError(f"FOFN file is empty: {fofn_path}")

    first_cols = [c.strip() for c in raw_lines[0][1].split("\t")]
    has_header = bool(first_cols) and first_cols[0].lower() == "path"

    if has_header:
        col_idx = {c.lower(): i for i, c in enumerate(first_cols)}
        path_idx = col_idx["path"]
        name_idx = col_idx.get("name")
        reads_idx = col_idx.get("reads")
        data_rows = raw_lines[1:]
    else:
        # Headerless: positional columns.  name/reads are read when the
        # row has enough columns and skipped otherwise.
        path_idx = 0
        name_idx = 1
        reads_idx = 2
        data_rows = raw_lines

    results: List[AssemblySpec] = []
    seen_names: dict[str, int] = {}

    for lineno, line in data_rows:
        fields = line.split("\t")
        if len(fields) <= path_idx:
            raise ValueError(
                f"FOFN line {lineno}: not enough columns "
                f"(expected >= {path_idx + 1}, got {len(fields)})"
            )

        # Resolve assembly path (relative to FOFN directory)
        raw_path = fields[path_idx].strip()
        asm_path = Path(raw_path)
        if not asm_path.is_absolute():
            asm_path = (fofn_dir / asm_path).resolve()
        else:
            asm_path = asm_path.resolve()
        if not asm_path.exists():
            raise ValueError(f"FOFN line {lineno}: assembly not found: {asm_path}")
        if not asm_path.is_file():
            raise ValueError(f"FOFN line {lineno}: not a file: {asm_path}")

        # Name
        if name_idx is not None and len(fields) > name_idx:
            name = fields[name_idx].strip()
        else:
            name = ""
        if not name:
            name = _strip_fasta_extension(asm_path.name)

        # Reads
        reads_path: Optional[Path] = None
        if reads_idx is not None and len(fields) > reads_idx:
            reads_raw = fields[reads_idx].strip()
            if reads_raw:
                reads_path = Path(reads_raw)
                if not reads_path.is_absolute():
                    reads_path = (fofn_dir / reads_path).resolve()
                else:
                    reads_path = reads_path.resolve()
                if not reads_path.exists():
                    raise ValueError(
                        f"FOFN line {lineno}: reads file not found: {reads_path}"
                    )

        # Track duplicates
        if name in seen_names:
            raise ValueError(
                f"FOFN has duplicate assembly name '{name}' "
                f"(lines {seen_names[name]} and {lineno})"
            )
        seen_names[name] = lineno

        results.append((asm_path, name, reads_path))

    if not results:
        raise ValueError(f"FOFN contains no assembly entries: {fofn_path}")

    return results


def scan_assembly_dir(dir_path: Path) -> List[Tuple[Path, str]]:
    """Scan a directory for FASTA files.

    Derives short names by stripping FASTA extensions.  If name collisions
    occur after stripping, the full filename stem is kept.

    Args:
        dir_path: Directory to scan.

    Returns:
        List of (assembly_path, name) tuples, sorted by name.

    Raises:
        ValueError: If no FASTA files found or the path is not a directory.
    """
    dir_path = Path(dir_path).resolve()
    if not dir_path.is_dir():
        raise ValueError(f"Not a directory: {dir_path}")

    fasta_files: List[Path] = []
    for entry in sorted(dir_path.iterdir()):
        if entry.is_file() and _has_fasta_extension(entry):
            fasta_files.append(entry)

    if not fasta_files:
        raise ValueError(f"No FASTA files found in directory: {dir_path}")

    # Derive names
    names = [_strip_fasta_extension(f.name) for f in fasta_files]

    # Handle collisions: if duplicates exist, use full filename stem instead
    if len(set(names)) < len(names):
        names = [f.stem for f in fasta_files]
        # If still duplicate (unlikely), append numeric suffix
        if len(set(names)) < len(names):
            seen: dict[str, int] = {}
            unique_names: list[str] = []
            for n in names:
                if n in seen:
                    seen[n] += 1
                    unique_names.append(f"{n}_{seen[n]}")
                else:
                    seen[n] = 0
                    unique_names.append(n)
            names = unique_names

    pairs = list(zip(fasta_files, names))
    pairs.sort(key=lambda x: x[1])
    return pairs


def resolve_assemblies(args) -> List[AssemblySpec]:
    """Resolve the list of assemblies from CLI arguments.

    Dispatches to :func:`parse_fofn` or :func:`scan_assembly_dir` based on
    which argument was provided.  For directory mode, ``--reads`` from the
    CLI applies as the global reads for all assemblies.  For FOFN mode,
    per-assembly reads from the FOFN take priority; ``--reads`` is the
    fallback for entries without a reads column/value.

    Args:
        args: Parsed argparse namespace (must have ``fofn``, ``assembly_dir``,
              and ``reads`` attributes).

    Returns:
        List of (assembly_path, name, reads_path) tuples.
    """
    global_reads: Optional[Path] = getattr(args, "reads", None)

    if getattr(args, "fofn", None):
        assemblies = parse_fofn(args.fofn)
        # Apply global --reads as fallback for entries without reads
        if global_reads:
            assemblies = [
                (p, n, r if r is not None else global_reads)
                for p, n, r in assemblies
            ]
        return assemblies

    if getattr(args, "assembly_dir", None):
        pairs = scan_assembly_dir(args.assembly_dir)
        return [(p, n, global_reads) for p, n in pairs]

    raise RuntimeError("resolve_assemblies called without --fofn or --assembly-dir")
