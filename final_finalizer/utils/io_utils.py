#!/usr/bin/env python3
"""
I/O utilities for final_finalizer.

Contains functions for file handling, gzip operations, and subprocess management.
"""
from __future__ import annotations

import gzip
import shutil
import stat as stat_module
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


def file_exists_and_valid(path: Path, min_size: int = 1) -> bool:
    """Check if a file exists, is readable, and has at least min_size bytes.

    This uses a single stat() call to avoid TOCTOU race conditions and detect
    corrupted/truncated files.

    Args:
        path: File path to check
        min_size: Minimum file size in bytes (default 1, i.e., non-empty)

    Returns:
        True if file exists, is a regular file, and meets size requirement
    """
    try:
        st = path.stat()
        if not stat_module.S_ISREG(st.st_mode):
            return False
        if st.st_size < min_size:
            return False
        return True
    except (OSError, PermissionError, FileNotFoundError):
        # File doesn't exist, was deleted, or we don't have permission
        return False


def run_to_gzip(cmd: list[str], gz_out: Path, err_path: Path, timeout: int = 7200) -> None:
    """Run a command, streaming stdout into a gzipped file. Stderr -> err_path.

    Args:
        cmd: Command and arguments to run
        gz_out: Output path for gzipped stdout
        err_path: Path for stderr output
        timeout: Maximum runtime in seconds (default 2 hours)

    Raises:
        RuntimeError: If command fails or times out
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
                try:
                    ret = proc.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                    raise RuntimeError(
                        f"Command timed out after {timeout}s: {' '.join(cmd)}\n"
                        f"Error log: {err_path}"
                    )

    if ret != 0:
        # Read error log preview for better error messages
        error_preview = ""
        try:
            if err_path.exists():
                with err_path.open("r", errors="replace") as ef:
                    lines = ef.readlines()[:10]
                    error_preview = "".join(lines)
        except Exception:
            pass

        raise RuntimeError(
            f"Command failed with return code {ret}\n"
            f"Command: {' '.join(cmd)}\n"
            f"Error log: {err_path}\n"
            f"Preview:\n{error_preview}"
        )


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
