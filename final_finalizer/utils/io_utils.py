#!/usr/bin/env python3
"""
I/O utilities for final_finalizer.

Contains functions for file handling, gzip operations, and subprocess management.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str):
    """
    Open plain or gzipped text/binary files based on suffix.
    mode should be compatible with open()/gzip.open().
    """
    if path.suffix == ".gz" or path.suffix == ".bgz":
        return gzip.open(path, mode)
    return open(path, mode)


def have_exe(exe: str) -> bool:
    """Check if an executable exists in PATH using shutil.which()."""
    return shutil.which(exe) is not None


def run_to_gzip(cmd: list[str], gz_out: Path, err_path: Path) -> None:
    """
    Run a command, streaming stdout into a gzipped file. Stderr -> err_path.
    """
    err_path.parent.mkdir(parents=True, exist_ok=True)

    with err_path.open("wb") as err_fh, gzip.open(gz_out, "wb") as gz_fh:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err_fh)
        try:
            assert proc.stdout is not None
            for chunk in iter(lambda: proc.stdout.read(1 << 20), b""):
                gz_fh.write(chunk)
        finally:
            try:
                if proc.stdout:
                    proc.stdout.close()
            finally:
                ret = proc.wait()

    if ret != 0:
        raise RuntimeError(f"Command failed with return code {ret} (see {err_path})")


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], int]:
    """
    Merge closed-open intervals [(s,e),...] where s<e.
    Returns merged list and total merged length.
    """
    if not intervals:
        return [], 0
    intervals = sorted(intervals)
    merged: list[tuple[int, int]] = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    total = sum(e - s for s, e in merged)
    return merged, total
