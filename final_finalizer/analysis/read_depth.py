#!/usr/bin/env python3
"""
Read depth analysis for final_finalizer.

Provides functions for:
- Detecting read format (FASTQ, unaligned BAM, aligned BAM/CRAM)
- Aligning reads to assembly with minimap2
- Running mosdepth for depth calculation
- Parsing mosdepth output and computing per-contig depth statistics
"""
from __future__ import annotations

import gzip
import shutil
import statistics
import subprocess
import sys
from enum import Enum
from pathlib import Path
from typing import Dict, Optional

from final_finalizer.models import DepthStats
from final_finalizer.utils.io_utils import file_exists_and_valid, have_exe


class ReadFormat(Enum):
    """Detected format of input reads."""
    FASTQ = "fastq"
    FASTQ_GZ = "fastq_gz"
    UNALIGNED_BAM = "unaligned_bam"
    UNALIGNED_CRAM = "unaligned_cram"
    ALIGNED_BAM = "aligned_bam"
    ALIGNED_CRAM = "aligned_cram"


def detect_read_format(reads_path: Path) -> ReadFormat:
    """Detect the format of input reads.

    Args:
        reads_path: Path to reads file

    Returns:
        ReadFormat enum indicating detected format

    Raises:
        ValueError: If format cannot be determined
    """
    suffix = reads_path.suffix.lower()
    suffixes = [s.lower() for s in reads_path.suffixes]

    # FASTQ detection
    if suffix in (".fq", ".fastq"):
        return ReadFormat.FASTQ
    if ".fq.gz" in "".join(suffixes) or ".fastq.gz" in "".join(suffixes):
        return ReadFormat.FASTQ_GZ
    if suffix == ".gz" and len(suffixes) >= 2:
        pre_suffix = suffixes[-2]
        if pre_suffix in (".fq", ".fastq"):
            return ReadFormat.FASTQ_GZ

    # BAM/CRAM detection
    if suffix == ".bam":
        if is_aligned_bam(reads_path):
            return ReadFormat.ALIGNED_BAM
        return ReadFormat.UNALIGNED_BAM

    if suffix == ".cram":
        if is_aligned_bam(reads_path):
            return ReadFormat.ALIGNED_CRAM
        return ReadFormat.UNALIGNED_CRAM

    raise ValueError(f"Cannot determine format of reads file: {reads_path}")


def is_aligned_bam(bam_path: Path) -> bool:
    """Check if a BAM/CRAM file contains aligned reads.

    Uses samtools view -H to check for @SQ header lines (reference sequences).
    Aligned BAM/CRAM files have @SQ headers; unaligned (uBAM) files do not.

    Args:
        bam_path: Path to BAM or CRAM file

    Returns:
        True if file contains aligned reads, False if unaligned
    """
    if not have_exe("samtools"):
        print("[warn] samtools not found; assuming BAM is unaligned", file=sys.stderr)
        return False

    try:
        result = subprocess.run(
            ["samtools", "view", "-H", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode != 0:
            print(f"[warn] samtools view -H failed for {bam_path}; assuming unaligned", file=sys.stderr)
            return False

        # Check for @SQ lines (reference sequences) - indicates aligned BAM
        for line in result.stdout.split("\n"):
            if line.startswith("@SQ"):
                return True
        return False

    except subprocess.TimeoutExpired:
        print(f"[warn] samtools timeout checking {bam_path}; assuming unaligned", file=sys.stderr)
        return False
    except Exception as e:
        print(f"[warn] Error checking BAM alignment status: {e}; assuming unaligned", file=sys.stderr)
        return False


def get_minimap2_preset(reads_type: str) -> str:
    """Get minimap2 preset for read type.

    Args:
        reads_type: One of "hifi", "ont", "sr"

    Returns:
        Minimap2 -x preset string
    """
    presets = {
        "hifi": "lr:hqae",  # Long reads, high quality, error rate < 1%
        "ont": "map-ont",   # ONT reads, standard accuracy
        "sr": "sr",         # Short reads (Illumina)
    }
    return presets.get(reads_type, "lr:hqae")


def align_reads_to_assembly(
    reads: Path,
    assembly: Path,
    output_bam: Path,
    threads: int,
    reads_type: str,
    work_dir: Path,
    err_path: Optional[Path] = None,
) -> bool:
    """Align reads to assembly using minimap2 and sort with samtools.

    Args:
        reads: Path to reads (FASTQ, FASTQ.gz, or unaligned BAM/CRAM)
        assembly: Path to assembly FASTA
        output_bam: Path for output sorted BAM
        threads: Number of threads
        reads_type: Read type for preset selection ("hifi", "ont", "sr")
        work_dir: Working directory for intermediate files
        err_path: Path for error log (optional)

    Returns:
        True if alignment succeeded, False otherwise
    """
    if file_exists_and_valid(output_bam, min_size=100):
        print(f"[info] Aligned BAM exists, reusing: {output_bam}", file=sys.stderr)
        return True

    work_dir.mkdir(parents=True, exist_ok=True)

    # Check for required tools
    mapper = None
    for exe in ["mm2plus", "minimap2"]:
        if have_exe(exe):
            mapper = exe
            break

    if not mapper:
        print("[error] Neither mm2plus nor minimap2 found in PATH", file=sys.stderr)
        return False

    if not have_exe("samtools"):
        print("[error] samtools not found in PATH", file=sys.stderr)
        return False

    preset = get_minimap2_preset(reads_type)
    err_file = err_path or (work_dir / "alignment.err")

    # Detect read format to handle different inputs
    read_format = detect_read_format(reads)

    # Build minimap2 command
    # -a outputs SAM, -x preset, -t threads
    mm2_cmd = [mapper, "-a", "-x", preset, "-t", str(threads), str(assembly)]

    if read_format in (ReadFormat.UNALIGNED_BAM, ReadFormat.UNALIGNED_CRAM):
        # For unaligned BAM/CRAM, use samtools fastq to convert to FASTQ on the fly
        # This avoids creating large intermediate files
        mm2_cmd.append("-")  # Read from stdin
    else:
        mm2_cmd.append(str(reads))

    # samtools sort command
    sort_cmd = ["samtools", "sort", "-@", str(max(1, threads // 2)), "-o", str(output_bam), "-"]

    print(f"[info] Aligning reads with {mapper} -x {preset}", file=sys.stderr)

    try:
        err_file.parent.mkdir(parents=True, exist_ok=True)
        with err_file.open("wb") as err_fh:
            if read_format in (ReadFormat.UNALIGNED_BAM, ReadFormat.UNALIGNED_CRAM):
                # Pipeline: samtools fastq | minimap2 | samtools sort
                fastq_cmd = ["samtools", "fastq", "-@", "2", str(reads)]

                p_fastq = subprocess.Popen(fastq_cmd, stdout=subprocess.PIPE, stderr=err_fh)
                p_mm2 = subprocess.Popen(mm2_cmd, stdin=p_fastq.stdout, stdout=subprocess.PIPE, stderr=err_fh)
                if p_fastq.stdout:
                    p_fastq.stdout.close()
                p_sort = subprocess.Popen(sort_cmd, stdin=p_mm2.stdout, stderr=err_fh)
                if p_mm2.stdout:
                    p_mm2.stdout.close()

                # Wait for pipeline to complete
                p_sort.wait()
                p_mm2.wait()
                p_fastq.wait()

                if p_sort.returncode != 0 or p_mm2.returncode != 0:
                    print(f"[error] Alignment pipeline failed (see {err_file})", file=sys.stderr)
                    return False
            else:
                # Pipeline: minimap2 | samtools sort
                p_mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=err_fh)
                p_sort = subprocess.Popen(sort_cmd, stdin=p_mm2.stdout, stderr=err_fh)
                if p_mm2.stdout:
                    p_mm2.stdout.close()

                p_sort.wait()
                p_mm2.wait()

                if p_sort.returncode != 0 or p_mm2.returncode != 0:
                    print(f"[error] Alignment pipeline failed (see {err_file})", file=sys.stderr)
                    return False

        # Index the BAM file
        idx_cmd = ["samtools", "index", "-@", str(threads), str(output_bam)]
        subprocess.run(idx_cmd, check=True, stderr=subprocess.DEVNULL)

        print(f"[done] Alignment complete: {output_bam}", file=sys.stderr)
        return True

    except Exception as e:
        print(f"[error] Alignment failed: {e}", file=sys.stderr)
        return False


def run_mosdepth(
    bam_path: Path,
    output_prefix: Path,
    threads: int,
    window_size: int = 1000,
    err_path: Optional[Path] = None,
) -> bool:
    """Run mosdepth to calculate depth statistics.

    Args:
        bam_path: Path to sorted, indexed BAM file
        output_prefix: Output prefix for mosdepth files
        threads: Number of threads
        window_size: Window size for depth calculation (default 1000)
        err_path: Path for error log (optional)

    Returns:
        True if mosdepth succeeded, False otherwise
    """
    regions_file = Path(str(output_prefix) + ".regions.bed.gz")
    if file_exists_and_valid(regions_file, min_size=50):
        print(f"[info] mosdepth output exists, reusing: {regions_file}", file=sys.stderr)
        return True

    if not have_exe("mosdepth"):
        print("[error] mosdepth not found in PATH", file=sys.stderr)
        return False

    # mosdepth command: window-based depth
    cmd = [
        "mosdepth",
        "--by", str(window_size),
        "-t", str(threads),
        "--no-per-base",  # Don't output per-base depth (saves space/time)
        str(output_prefix),
        str(bam_path),
    ]

    print(f"[info] Running mosdepth with {window_size}bp windows", file=sys.stderr)

    try:
        err_file = err_path or Path(str(output_prefix) + ".mosdepth.err")
        err_file.parent.mkdir(parents=True, exist_ok=True)

        with err_file.open("wb") as err_fh:
            result = subprocess.run(cmd, stderr=err_fh, check=False)

        if result.returncode != 0:
            print(f"[error] mosdepth failed (see {err_file})", file=sys.stderr)
            return False

        print(f"[done] mosdepth complete: {regions_file}", file=sys.stderr)
        return True

    except Exception as e:
        print(f"[error] mosdepth failed: {e}", file=sys.stderr)
        return False


def parse_mosdepth_regions(
    regions_bed_gz: Path,
    contig_lengths: Dict[str, int],
) -> Dict[str, DepthStats]:
    """Parse mosdepth regions.bed.gz output and compute per-contig depth stats.

    Args:
        regions_bed_gz: Path to mosdepth *.regions.bed.gz file
        contig_lengths: Dictionary mapping contig names to lengths

    Returns:
        Dictionary mapping contig names to DepthStats objects
    """
    # Collect depths per contig (list of window depths)
    contig_depths: Dict[str, list] = {name: [] for name in contig_lengths.keys()}
    contig_windows: Dict[str, list] = {name: [] for name in contig_lengths.keys()}

    try:
        with gzip.open(regions_bed_gz, "rt") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue

                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                depth = float(parts[3])
                window_size = end - start

                if chrom in contig_depths:
                    contig_depths[chrom].append(depth)
                    contig_windows[chrom].append((start, end, depth, window_size))

    except Exception as e:
        print(f"[warn] Error parsing mosdepth output {regions_bed_gz}: {e}", file=sys.stderr)
        return {}

    # Compute stats per contig
    result: Dict[str, DepthStats] = {}

    for contig, depths in contig_depths.items():
        if not depths:
            # No coverage data - use zeros
            result[contig] = DepthStats(
                mean_depth=0.0,
                median_depth=0.0,
                std_depth=0.0,
                min_depth=0.0,
                max_depth=0.0,
                breadth_1x=0.0,
                breadth_10x=0.0,
            )
            continue

        windows = contig_windows[contig]
        total_bp = sum(ws for _, _, _, ws in windows)
        contig_len = contig_lengths.get(contig, total_bp)

        # Weighted mean by window size
        weighted_sum = sum(d * ws for _, _, d, ws in windows)
        mean_depth = weighted_sum / total_bp if total_bp > 0 else 0.0

        # Median (simple approach - use raw depths, which is approximately weighted)
        median_depth = statistics.median(depths) if depths else 0.0

        # Standard deviation
        std_depth = statistics.stdev(depths) if len(depths) > 1 else 0.0

        # Min/Max
        min_depth = min(depths)
        max_depth = max(depths)

        # Breadth calculations: fraction of bases with >= threshold coverage
        bp_1x = sum(ws for _, _, d, ws in windows if d >= 1.0)
        bp_10x = sum(ws for _, _, d, ws in windows if d >= 10.0)

        breadth_1x = bp_1x / contig_len if contig_len > 0 else 0.0
        breadth_10x = bp_10x / contig_len if contig_len > 0 else 0.0

        result[contig] = DepthStats(
            mean_depth=mean_depth,
            median_depth=median_depth,
            std_depth=std_depth,
            min_depth=min_depth,
            max_depth=max_depth,
            breadth_1x=min(1.0, breadth_1x),  # Cap at 1.0
            breadth_10x=min(1.0, breadth_10x),
        )

    return result


def calculate_depth_metrics(
    reads: Path,
    assembly: Path,
    contig_lengths: Dict[str, int],
    work_dir: Path,
    threads: int,
    reads_type: str = "hifi",
    window_size: int = 1000,
) -> Dict[str, DepthStats]:
    """Calculate read depth metrics for an assembly.

    This is the main entry point for depth analysis. It handles:
    1. Format detection for input reads
    2. Alignment if reads are not pre-aligned
    3. Running mosdepth for depth calculation
    4. Parsing results into per-contig DepthStats

    Args:
        reads: Path to reads file (FASTQ, BAM, CRAM)
        assembly: Path to assembly FASTA
        contig_lengths: Dictionary of contig names to lengths
        work_dir: Working directory for intermediate files
        threads: Number of threads to use
        reads_type: Read type for minimap2 preset ("hifi", "ont", "sr")
        window_size: Window size for mosdepth (default 1000)

    Returns:
        Dictionary mapping contig names to DepthStats objects
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Detect read format
    try:
        read_format = detect_read_format(reads)
        print(f"[info] Detected read format: {read_format.value}", file=sys.stderr)
    except ValueError as e:
        print(f"[error] {e}", file=sys.stderr)
        return {}

    # Determine if we need to align
    if read_format in (ReadFormat.ALIGNED_BAM, ReadFormat.ALIGNED_CRAM):
        # Already aligned - use directly
        bam_path = reads
        print("[info] Using pre-aligned BAM/CRAM", file=sys.stderr)

        # Check if index exists, create if needed
        bai_path = Path(str(reads) + ".bai")
        csi_path = Path(str(reads) + ".csi")
        if not bai_path.exists() and not csi_path.exists():
            print("[info] Indexing BAM file", file=sys.stderr)
            if have_exe("samtools"):
                subprocess.run(
                    ["samtools", "index", "-@", str(threads), str(reads)],
                    check=True,
                    stderr=subprocess.DEVNULL,
                )
    else:
        # Need to align
        bam_path = work_dir / "aligned_reads.sorted.bam"
        err_path = work_dir / "alignment.err"

        success = align_reads_to_assembly(
            reads=reads,
            assembly=assembly,
            output_bam=bam_path,
            threads=threads,
            reads_type=reads_type,
            work_dir=work_dir,
            err_path=err_path,
        )

        if not success:
            print("[error] Read alignment failed", file=sys.stderr)
            return {}

    # Run mosdepth
    mosdepth_prefix = work_dir / "depth"
    err_path = work_dir / "mosdepth.err"

    success = run_mosdepth(
        bam_path=bam_path,
        output_prefix=mosdepth_prefix,
        threads=threads,
        window_size=window_size,
        err_path=err_path,
    )

    if not success:
        print("[error] mosdepth failed", file=sys.stderr)
        return {}

    # Parse mosdepth output
    regions_file = Path(str(mosdepth_prefix) + ".regions.bed.gz")
    if not regions_file.exists():
        print(f"[error] mosdepth output not found: {regions_file}", file=sys.stderr)
        return {}

    depth_stats = parse_mosdepth_regions(regions_file, contig_lengths)
    print(f"[info] Computed depth stats for {len(depth_stats)} contigs", file=sys.stderr)

    return depth_stats
