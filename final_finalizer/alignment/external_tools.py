#!/usr/bin/env python3
"""
External alignment tool runners for final_finalizer.

Contains functions to run minimap2, mm2plus, gffread, and miniprot.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import List, Optional

from final_finalizer.utils.io_utils import have_exe, run_to_gzip


# ----------------------------
# Minimap2/mm2plus detection and helpers
# ----------------------------
# mm2plus is a drop-in replacement for minimap2 with additional features.
# We prefer mm2plus if available, otherwise fall back to minimap2.
_MINIMAP2_EXE: Optional[str] = None


def get_minimap2_exe() -> Optional[str]:
    """Detect and cache the minimap2/mm2plus executable.

    Prefers mm2plus if available, falls back to minimap2.
    Returns None if neither is found.
    """
    global _MINIMAP2_EXE
    if _MINIMAP2_EXE is None:
        if have_exe("mm2plus"):
            _MINIMAP2_EXE = "mm2plus"
        elif have_exe("minimap2"):
            _MINIMAP2_EXE = "minimap2"
        else:
            _MINIMAP2_EXE = ""  # Mark as checked but not found
    return _MINIMAP2_EXE if _MINIMAP2_EXE else None


def run_minimap2(
    ref: Path,
    qry: Path,
    paf_out: Path,
    threads: int,
    preset: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    gzip_output: bool = False,
    err_path: Optional[Path] = None,
) -> bool:
    """Run minimap2/mm2plus alignment.

    Uses mm2plus if available, otherwise minimap2. Both are functionally
    equivalent for our use cases.

    Args:
        ref: Reference FASTA path
        qry: Query FASTA path
        paf_out: Output PAF path (will be gzipped if gzip_output=True)
        threads: Number of threads
        preset: Minimap2 preset (e.g., "asm5", "asm20")
        extra_args: Additional command-line arguments
        gzip_output: If True, compress output with gzip
        err_path: Path for stderr output (optional)

    Returns:
        True if alignment succeeded, False otherwise
    """
    # Check if output already exists
    if paf_out.exists():
        print(f"[info] PAF exists, reusing: {paf_out}", file=sys.stderr)
        return True

    mapper = get_minimap2_exe()
    if not mapper:
        print("[warn] Neither mm2plus nor minimap2 found in PATH", file=sys.stderr)
        return False

    # Build command
    cmd = [mapper, "-t", str(threads)]
    if preset:
        cmd += ["-x", preset]
    if extra_args:
        cmd.extend(extra_args)

    if gzip_output:
        # Output to stdout for gzip compression
        cmd += [str(ref), str(qry)]
    else:
        # Output directly to file
        cmd += ["-o", str(paf_out), str(ref), str(qry)]

    print(f"[info] Running {mapper} -> {paf_out}", file=sys.stderr)

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if gzip_output:
            run_to_gzip(cmd, paf_out, err_path or Path("/dev/null"))
        else:
            if err_path:
                with err_path.open("wb") as err_fh:
                    ret = subprocess.run(cmd, stderr=err_fh, check=False)
            else:
                ret = subprocess.run(cmd, stderr=subprocess.DEVNULL, check=False)
            if ret.returncode != 0:
                print(f"[warn] {mapper} failed with return code {ret.returncode}", file=sys.stderr)
                return False
    except Exception as e:
        print(f"[warn] {mapper} failed: {e}", file=sys.stderr)
        return False

    return True


def run_minimap2_synteny(
    ref: Path,
    qry: Path,
    paf_gz_out: Path,
    threads: int,
    preset: Optional[str],
    k: Optional[int],
    w: Optional[int],
    err_path: Path,
) -> None:
    """Run minimap2/mm2plus with parameterization tuned for cross-species synteny alignments.

    This uses specialized parameters (-c -f1 -B2 -O1,24 -s 10 --max-chain-skip=80 -z 2000,200)
    for genome-wide synteny analysis. Output is gzipped PAF.
    """
    if paf_gz_out.exists():
        print(f"[info] PAF.gz exists, reusing: {paf_gz_out}", file=sys.stderr)
        return

    mapper = get_minimap2_exe()
    if not mapper:
        raise RuntimeError("Neither mm2plus nor minimap2 found in PATH")

    cmd = [
        mapper,
        "-c",
        "-f1",
        "-B2",
        "-O1,24",
        "-s", "10",
        "--max-chain-skip=80",
        "-z", "2000,200",
        "-t", str(threads),
    ]
    if preset:
        cmd += ["-x", preset]
    if k is not None:
        cmd += ["-k", str(k)]
    if w is not None:
        cmd += ["-w", str(w)]
    cmd += [str(ref), str(qry)]

    print(
        f"[info] Running {mapper} -> {paf_gz_out} (stderr: {err_path}): " + " ".join(cmd),
        file=sys.stderr,
    )
    run_to_gzip(cmd, paf_gz_out, err_path)


def run_gffread_extract_proteins(
    gffread_exe: str,
    ref_fasta: Path,
    ref_gff3: Path,
    proteins_faa: Path,
    err_path: Path,
) -> None:
    """
    Use gffread to extract protein sequences from reference genome + GFF3.
      gffread -g ref.fa -y proteins.fa ref.gff3
    """
    if proteins_faa.exists():
        print(f"[info] Proteins FASTA exists, reusing: {proteins_faa}", file=sys.stderr)
        return

    if not have_exe(gffread_exe):
        raise RuntimeError(f"gffread executable not found/usable: {gffread_exe}")

    cmd = [gffread_exe, "-g", str(ref_fasta), "-y", str(proteins_faa), str(ref_gff3)]
    print(f"[info] Extracting proteins with gffread (stderr: {err_path}): " + " ".join(cmd), file=sys.stderr)

    err_path.parent.mkdir(parents=True, exist_ok=True)
    with err_path.open("wb") as err_fh:
        ret = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=err_fh)

    if ret != 0:
        raise RuntimeError(f"gffread failed with return code {ret} (see {err_path})")

    if not proteins_faa.exists() or proteins_faa.stat().st_size == 0:
        raise RuntimeError(f"gffread produced empty proteins FASTA: {proteins_faa}")


def run_miniprot(
    miniprot_exe: str,
    query_fasta: Path,
    proteins_faa: Path,
    paf_gz_out: Path,
    threads: int,
    extra_args: str,
    err_path: Path,
) -> None:
    """
    Align proteins to query genome with miniprot.
    Typical usage:
      miniprot -t N <target_genome.fa> <query_proteins.faa> > out.paf
    Here target = query genome, query = reference proteins.
    Output is gzipped PAF.
    """
    if paf_gz_out.exists():
        print(f"[info] miniprot PAF.gz exists, reusing: {paf_gz_out}", file=sys.stderr)
        return

    if not have_exe(miniprot_exe):
        raise RuntimeError(f"miniprot executable not found/usable: {miniprot_exe}")

    cmd = [miniprot_exe, "-t", str(threads)]
    if extra_args:
        cmd += extra_args.strip().split()
    cmd += [str(query_fasta), str(proteins_faa)]

    print(
        f"[info] Running miniprot -> {paf_gz_out} (stderr: {err_path}): " + " ".join(cmd),
        file=sys.stderr,
    )
    run_to_gzip(cmd, paf_gz_out, err_path)
