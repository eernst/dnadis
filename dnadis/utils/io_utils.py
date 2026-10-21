#!/usr/bin/env python3
"""
I/O utilities for dnadis.

Contains functions for file handling, gzip operations, and subprocess management.
"""
from __future__ import annotations

import gzip
import os
import re
import shutil
import stat as stat_module
import struct
import subprocess
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Optional, Union


_FILENAME_UNSAFE_RE = re.compile(r"[\s/\\]+")


def sanitize_filename_component(name: str, fallback: str = "comparison") -> str:
    """Return a filename-safe version of ``name``.

    Replaces runs of whitespace and path separators with a single ``_``,
    strips leading/trailing ``_``, and falls back to ``fallback`` if the
    result is empty.  Other characters are left intact so users can keep
    Unicode or punctuation that downstream tools handle (rmarkdown, in
    particular, mishandles spaces in input filenames — those are the
    ones we have to fix).
    """
    cleaned = _FILENAME_UNSAFE_RE.sub("_", name).strip("_")
    return cleaned or fallback



def open_maybe_gzip(path: Path, mode: str):
    """
    Open plain or gzipped text/binary files based on suffix.
    mode should be compatible with open()/gzip.open().
    """
    if is_gzip_path(path):
        return gzip.open(path, mode)
    return open(path, mode)


# BGZF (blocked gzip, as written by htslib's bgzip): a series of gzip members
# of at most 64 KiB uncompressed each, carrying their compressed size in a
# 'BC' extra field, terminated by an empty EOF member.  Plain gzip readers
# decompress it as one stream; samtools faidx can index it for random access.
_BGZF_BLOCK_SIZE = 0xFF00  # uncompressed bytes per block (htslib default)
_BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def _bgzf_block(data: bytes, level: int) -> bytes:
    """Compress *data* (<= 64 KiB) into one BGZF block."""
    comp = zlib.compressobj(level, zlib.DEFLATED, -15)
    cdata = comp.compress(data) + comp.flush()
    header = struct.pack(
        "<4BI2BH2BHH",
        0x1F, 0x8B, 8, 4,     # magic, deflate, FEXTRA
        0,                    # mtime 0 keeps output reproducible
        0, 0xFF,              # XFL, OS unknown
        6,                    # XLEN
        ord("B"), ord("C"), 2,
        len(cdata) + 25,      # BSIZE = total block size - 1
    )
    trailer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + cdata + trailer


class BgzfWriter:
    """Write text to a BGZF-compressed file, compressing blocks in parallel.

    zlib releases the GIL, so blocks compress concurrently on *threads*
    threads.  Usable as a context manager; the EOF block is written on close.
    """

    def __init__(self, path: Path, threads: int = 1, level: int = 6) -> None:
        self._fh = open(path, "wb")
        self._level = level
        self._buf = bytearray()
        self._threads = max(1, threads)
        self._batch = _BGZF_BLOCK_SIZE * 16 * self._threads
        self._pool: Optional[ThreadPoolExecutor] = (
            ThreadPoolExecutor(max_workers=self._threads) if self._threads > 1 else None
        )

    def write(self, text: Union[str, bytes]) -> None:
        self._buf += text.encode() if isinstance(text, str) else text
        if len(self._buf) >= self._batch:
            self._flush(final=False)

    def _flush(self, final: bool) -> None:
        n = len(self._buf) if final else len(self._buf) - len(self._buf) % _BGZF_BLOCK_SIZE
        chunks: List[bytes] = [
            bytes(self._buf[i:i + _BGZF_BLOCK_SIZE]) for i in range(0, n, _BGZF_BLOCK_SIZE)
        ]
        del self._buf[:n]
        if self._pool is not None:
            blocks = self._pool.map(lambda c: _bgzf_block(c, self._level), chunks)
        else:
            blocks = (_bgzf_block(c, self._level) for c in chunks)
        for block in blocks:
            self._fh.write(block)

    def close(self) -> None:
        if self._fh.closed:
            return
        try:
            self._flush(final=True)
            self._fh.write(_BGZF_EOF)
        finally:
            self._fh.close()
            if self._pool is not None:
                self._pool.shutdown()

    def __enter__(self) -> "BgzfWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def is_gzip_path(path: Path) -> bool:
    """True if *path* names a gzip/BGZF file by suffix."""
    return path.suffix in (".gz", ".bgz")


def open_text_writer(path: Path, threads: int = 1):
    """Open *path* for text writing: BGZF if it ends in .gz/.bgz, else plain."""
    if is_gzip_path(path):
        return BgzfWriter(path, threads=threads)
    return open(path, "w")


def symlink_if_needed(link: Path, target: Path) -> bool:
    """Make *link* a relative symlink to *target* (same directory).

    Replaces whatever is at *link* unless it already points at *target*.
    Returns True if the link was (re)created.
    """
    try:
        if link.is_symlink() and os.readlink(link) == target.name:
            return False
    except OSError:
        pass
    tmp = link.with_name(link.name + ".tmp")
    if tmp.is_symlink() or tmp.exists():
        tmp.unlink()
    os.symlink(target.name, tmp)
    tmp.replace(link)
    return True


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


def output_is_current(output: Path, *inputs: Path) -> bool:
    """True if *output* is a valid file at least as new as every input.

    Used to reuse a cached tool output only when its inputs have not been
    rewritten since.  A missing input counts as changed.
    """
    if not file_exists_and_valid(output):
        return False
    try:
        out_mtime = output.stat().st_mtime_ns
        return all(Path(p).stat().st_mtime_ns <= out_mtime for p in inputs)
    except OSError:
        return False


def write_text_if_changed(path: Path, text: str) -> bool:
    """Write *text* to *path* unless the file already holds exactly that text.

    Leaving unchanged files untouched keeps their mtimes, so downstream
    outputs checked with :func:`output_is_current` stay reusable.  Returns
    True if the file was written.
    """
    try:
        if path.read_text() == text:
            return False
    except (OSError, UnicodeDecodeError):
        pass
    path.write_text(text)
    return True


def run_to_gzip(cmd: list[str], gz_out: Path, err_path: Path, timeout: int = 7200) -> None:
    """Run a command, streaming stdout into a gzipped file. Stderr -> err_path.

    Writes to a temporary file first, then atomically renames on success.
    This prevents truncated gzip files from being left on disk when the
    process is killed (e.g. SLURM timeout, OOM).

    Args:
        cmd: Command and arguments to run
        gz_out: Output path for gzipped stdout
        err_path: Path for stderr output
        timeout: Maximum runtime in seconds (default 2 hours)

    Raises:
        RuntimeError: If command fails or times out
    """
    err_path.parent.mkdir(parents=True, exist_ok=True)
    gz_out.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = gz_out.with_suffix(gz_out.suffix + ".tmp")

    try:
        with err_path.open("wb") as err_fh, gzip.open(tmp_out, "wb") as gz_fh:
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

        tmp_out.rename(gz_out)
    except BaseException:
        tmp_out.unlink(missing_ok=True)
        raise


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
