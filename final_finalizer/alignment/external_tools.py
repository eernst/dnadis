#!/usr/bin/env python3
"""
External alignment tool runners for final_finalizer.

Contains functions to run minimap2, mm2plus, gffread, and miniprot.
"""
from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional

from final_finalizer.utils.io_utils import file_exists_and_valid, have_exe, run_to_gzip
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("external_tools")


# Pattern to detect dangerous shell metacharacters that could enable command injection
_DANGEROUS_SHELL_CHARS = re.compile(r'[;&|`$(){}[\]<>!\\\'\"#]')


def validate_extra_args(args_str: str, arg_name: str = "extra_args") -> str:
    """Validate extra command-line arguments for shell safety.

    Checks for dangerous shell metacharacters that could enable command injection.
    This is a defense-in-depth measure since shlex.split() is used, but provides
    an additional layer of protection and clearer error messages.

    Args:
        args_str: The argument string to validate
        arg_name: Name of the argument for error messages

    Returns:
        The validated and stripped argument string

    Raises:
        ValueError: If dangerous shell metacharacters are detected
    """
    if not args_str:
        return ""
    args_str = args_str.strip()
    if _DANGEROUS_SHELL_CHARS.search(args_str):
        raise ValueError(f"Invalid characters in {arg_name}: shell metacharacters not allowed")
    return args_str


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
    # Check if output already exists and is valid
    if file_exists_and_valid(paf_out):
        logger.info(f"PAF exists, reusing: {paf_out}")
        return True

    mapper = get_minimap2_exe()
    if not mapper:
        logger.warning("Neither mm2plus nor minimap2 found in PATH")
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

    logger.info(f"Running {mapper} -> {paf_out}")

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
                logger.warning(f"{mapper} failed with return code {ret.returncode}")
                return False
    except Exception as e:
        logger.warning(f"{mapper} failed: {e}")
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
    permissive: bool = False,
) -> None:
    """Run minimap2/mm2plus for whole-genome synteny alignment (nucleotide mode).

    This function is used exclusively for --synteny-mode nucleotide, which performs
    chromosome-scale structural composition analysis via whole-genome nucleotide
    alignment. The permissive chaining parameters create megabase-scale synteny
    blocks suitable for chromosome classification and architecture visualization.

    Args:
        ref: Reference genome FASTA
        qry: Query assembly FASTA
        paf_gz_out: Output PAF file (will be gzipped)
        threads: Number of threads for minimap2
        preset: Minimap2 preset (e.g., asm20 for cross-species, asm5 for same-species)
        k: K-mer size override (optional)
        w: Minimizer window size override (optional)
        err_path: Path for stderr output
        permissive: If True, use permissive chaining (default for nucleotide mode).
            If False, use conservative parameters (legacy, not used in final_finalizer).

    Permissive mode (default for nucleotide mode):
        Parameters: -c -f1 -B2 -O1,24 -s 10 --max-chain-skip=300 -z 10000,1000 -r 50000
        - Chains through repetitive regions and small gaps (kb-scale)
        - Creates continuous megabase-scale blocks for chromosome architecture analysis
        - Suitable for both within-species and cross-species comparisons
        - Balanced by downstream filtering (identity, length, gate thresholds)

    Conservative mode (legacy, not used):
        Parameters: -c -f1 -B2 -O1,24 -s 10 --max-chain-skip=80 -z 2000,200
        - Fragments alignments at ambiguous repeats
        - Preserves fine-scale structural variant details
        - Not suitable for chromosome classification (too fragmented)

    Output:
        Gzipped PAF file with whole-genome alignment records

    Raises:
        RuntimeError: If neither mm2plus nor minimap2 is available
    """
    if file_exists_and_valid(paf_gz_out):
        logger.info(f"PAF.gz exists, reusing: {paf_gz_out}")
        return

    mapper = get_minimap2_exe()
    if not mapper:
        raise RuntimeError("Neither mm2plus nor minimap2 found in PATH")

    if permissive:
        # Permissive chaining for within-species structural analysis
        cmd = [
            mapper,
            "-c",  # Generate CIGAR
            "-f1",  # Filter: keep all alignments
            "-B2",  # Mismatch penalty (low for nearly identical sequences)
            "-O1,24",  # Gap open penalty (low extension, high open)
            "-s", "10",  # Minimal chaining score
            "--max-chain-skip=300",  # Allow skipping many seeds (chain through repeats)
            "-z", "10000,1000",  # High Z-drop (tolerate long gaps in chain)
            "-r", "50000",  # Large bandwidth (50kb - tolerate structural variations)
            "-t", str(threads),
        ]
    else:
        # Conservative chaining for cross-species comparisons
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

    mode = "permissive" if permissive else "conservative"
    logger.info(f"Running {mapper} ({mode} mode) -> {paf_gz_out} (stderr: {err_path})")
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
    if file_exists_and_valid(proteins_faa):
        logger.info(f"Proteins FASTA exists, reusing: {proteins_faa}")
        return

    if not have_exe(gffread_exe):
        raise RuntimeError(f"gffread executable not found/usable: {gffread_exe}")

    cmd = [gffread_exe, "-g", str(ref_fasta), "-y", str(proteins_faa), str(ref_gff3)]
    logger.info(f"Extracting proteins with gffread (stderr: {err_path}): " + " ".join(cmd))

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
    if file_exists_and_valid(paf_gz_out):
        logger.info(f"miniprot PAF.gz exists, reusing: {paf_gz_out}")
        return

    if not have_exe(miniprot_exe):
        raise RuntimeError(f"miniprot executable not found/usable: {miniprot_exe}")

    cmd = [miniprot_exe, "-t", str(threads)]
    if extra_args:
        # Validate for shell safety and use shlex.split for proper argument parsing
        validated = validate_extra_args(extra_args, "--miniprot-args")
        cmd += shlex.split(validated)
    cmd += [str(query_fasta), str(proteins_faa)]

    logger.info(f"Running miniprot -> {paf_gz_out} (stderr: {err_path}): " + " ".join(cmd))
    run_to_gzip(cmd, paf_gz_out, err_path)


def run_cdhit_est(
    input_fasta: Path,
    output_prefix: Path,
    identity: float = 0.95,
    threads: int = 1,
    err_path: Optional[Path] = None,
) -> bool:
    """Run cd-hit-est to cluster nucleotide sequences.

    Args:
        input_fasta: Input FASTA file with sequences to cluster
        output_prefix: Output prefix (creates prefix and prefix.clstr)
        identity: Sequence identity threshold (0.0-1.0) [0.95]
        threads: Number of threads
        err_path: Path for stderr output

    Returns:
        True if successful, False if cd-hit-est not available
    """
    if not have_exe("cd-hit-est"):
        return False

    if file_exists_and_valid(output_prefix):
        logger.info(f"cd-hit-est output exists, reusing: {output_prefix}")
        return True

    # Choose word size based on identity threshold (cd-hit-est requirement)
    if identity >= 0.90:
        word_size = 8
    elif identity >= 0.88:
        word_size = 7
    elif identity >= 0.85:
        word_size = 6
    elif identity >= 0.80:
        word_size = 5
    else:
        word_size = 4

    cmd = [
        "cd-hit-est",
        "-i", str(input_fasta),
        "-o", str(output_prefix),
        "-c", str(identity),
        "-n", str(word_size),
        "-T", str(threads),
        "-M", "0",  # Unlimited memory
        "-d", "0",  # Full sequence name in output
    ]

    logger.info(f"Running cd-hit-est (identity={identity}) -> {output_prefix}")

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if err_path:
            with err_path.open("wb") as err_fh:
                ret = subprocess.run(cmd, stderr=err_fh, stdout=subprocess.DEVNULL, check=False)
        else:
            ret = subprocess.run(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, check=False)
        if ret.returncode != 0:
            logger.warning(f"cd-hit-est failed with return code {ret.returncode}")
            return False
    except Exception as e:
        logger.warning(f"cd-hit-est failed: {e}")
        return False

    return True


def run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    threads: int = 1,
    err_path: Optional[Path] = None,
) -> bool:
    """Run MAFFT multiple sequence alignment.

    Args:
        input_fasta: Input FASTA with sequences to align
        output_fasta: Output aligned FASTA
        threads: Number of threads
        err_path: Path for stderr output

    Returns:
        True if successful, False if MAFFT not available
    """
    if not have_exe("mafft"):
        return False

    if file_exists_and_valid(output_fasta):
        logger.info(f"MAFFT output exists, reusing: {output_fasta}")
        return True

    cmd = [
        "mafft",
        "--auto",
        "--thread", str(threads),
        str(input_fasta),
    ]

    logger.info(f"Running MAFFT -> {output_fasta}")

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if err_path:
            with output_fasta.open("w") as out_fh, err_path.open("wb") as err_fh:
                ret = subprocess.run(cmd, stdout=out_fh, stderr=err_fh, check=False)
        else:
            with output_fasta.open("w") as out_fh:
                ret = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.DEVNULL, check=False)
        if ret.returncode != 0:
            logger.warning(f"MAFFT failed with return code {ret.returncode}")
            return False
    except Exception as e:
        logger.warning(f"MAFFT failed: {e}")
        return False

    return True
