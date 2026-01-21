#!/usr/bin/env python3
"""
Read depth analysis for final_finalizer.

Provides functions for:
- Detecting read format (FASTQ, unaligned BAM, aligned BAM/CRAM)
- Downsampling reads to target coverage (seqtk for FASTQ, samtools for BAM)
- Aligning reads to assembly with minimap2/mm2plus
- Running mosdepth for depth calculation
- Parsing mosdepth output and computing per-contig depth statistics

The depth analysis pipeline:
1. Auto-detect read format (FASTQ/BAM/CRAM, aligned vs unaligned)
2. If not pre-aligned:
   a. Estimate total read bases and coverage
   b. Downsample to target coverage if needed (default 20X)
   c. Align with minimap2 using appropriate preset
3. Run mosdepth to compute window-based depth metrics
4. Parse results into DepthStats per contig (mean, median, std, breadth)

Depth metrics are useful for:
- Quality assessment (breadth of coverage)
- Copy number evaluation (relative depth)
- Identifying misassemblies (depth discontinuities)
"""
from __future__ import annotations

import gzip
import random
import shutil
import statistics
import subprocess
import sys
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Tuple

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
        reads_type: One of "hifi_onthq", "ont", "sr"

    Returns:
        Minimap2 -x preset string

    Notes:
        - "hifi_onthq" uses lr:hqae (long reads, high quality, all-vs-all, end-to-end)
          Suitable for PacBio HiFi, Duplex, and ONT Q20+ with error rate < 1%
        - "ont" uses map-ont for standard-accuracy ONT reads
        - "sr" uses sr for Illumina short reads
    """
    presets = {
        "hifi_onthq": "lr:hqae",  # Long reads, high quality, error rate < 1% (HiFi, Duplex, ONT Q20+)
        "ont": "map-ont",         # ONT reads, standard accuracy
        "sr": "sr",               # Short reads (Illumina)
    }
    return presets.get(reads_type, "lr:hqae")


def estimate_read_bases_fastq(reads_path: Path) -> int:
    """Estimate total bases in a FASTQ file using seqtk fqchk.

    Falls back to sampling-based estimation if seqtk is not available.

    Args:
        reads_path: Path to FASTQ or FASTQ.gz file

    Returns:
        Estimated total bases in the file
    """
    if have_exe("seqtk"):
        try:
            # seqtk fqchk outputs stats including total bases
            # Format: min_len: X; max_len: Y; avg_len: Z; N reads; Q bases
            result = subprocess.run(
                ["seqtk", "fqchk", str(reads_path)],
                capture_output=True,
                text=True,
                timeout=300,  # 5 min timeout for large files
            )
            if result.returncode == 0:
                # Parse the first line which contains summary stats
                # Example: "min_len: 1000; max_len: 50000; avg_len: 15234.5; 1234567 reads; 18789012345 bases"
                for line in result.stdout.split("\n"):
                    if "bases" in line.lower():
                        # Try to extract the number before "bases"
                        parts = line.split(";")
                        for part in parts:
                            if "bases" in part.lower():
                                # Extract number from string like "18789012345 bases"
                                tokens = part.strip().split()
                                for i, tok in enumerate(tokens):
                                    if "bases" in tok.lower() and i > 0:
                                        try:
                                            return int(tokens[i - 1])
                                        except ValueError:
                                            pass
                                    # Also try format "N bases"
                                    try:
                                        return int(tok)
                                    except ValueError:
                                        continue
        except subprocess.TimeoutExpired:
            print("[warn] seqtk fqchk timed out, falling back to sampling", file=sys.stderr)
        except Exception as e:
            print(f"[warn] seqtk fqchk failed: {e}, falling back to sampling", file=sys.stderr)

    # Fallback: sample first N reads and extrapolate
    return _estimate_bases_by_sampling_fastq(reads_path)


def _estimate_bases_by_sampling_fastq(reads_path: Path, sample_reads: int = 10000) -> int:
    """Estimate total bases by sampling first N reads from FASTQ.

    Args:
        reads_path: Path to FASTQ or FASTQ.gz file
        sample_reads: Number of reads to sample

    Returns:
        Estimated total bases
    """
    total_bases = 0
    read_count = 0
    line_count = 0
    bytes_read = 0

    try:
        opener = gzip.open if str(reads_path).endswith(".gz") else open
        with opener(reads_path, "rt") as fh:
            for line in fh:
                bytes_read += len(line)
                line_count += 1
                # FASTQ format: line 2 of every 4 is the sequence
                if line_count % 4 == 2:
                    total_bases += len(line.strip())
                    read_count += 1
                    if read_count >= sample_reads:
                        break

        if read_count == 0:
            return 0

        avg_read_len = total_bases / read_count

        # Estimate total reads by file size ratio
        # For gzipped files, use compressed file size and estimate compression ratio
        file_size = reads_path.stat().st_size
        if str(reads_path).endswith(".gz"):
            # Estimate compression ratio (~3-4x for FASTQ)
            # Use bytes_read (uncompressed) to estimate reads, then scale by compression
            bytes_per_read_uncompressed = bytes_read / read_count
            # Typical gzip compression for FASTQ is ~3.5x
            compression_ratio = 3.5
            estimated_total_reads = (file_size * compression_ratio) / bytes_per_read_uncompressed
        else:
            bytes_per_read = bytes_read / read_count
            estimated_total_reads = file_size / bytes_per_read

        estimated_total_bases = int(estimated_total_reads * avg_read_len)

        print(f"[info] Estimated {estimated_total_reads:.0f} reads, {estimated_total_bases:,} bases (sampling-based)", file=sys.stderr)
        return estimated_total_bases

    except Exception as e:
        print(f"[warn] Could not estimate read bases: {e}", file=sys.stderr)
        return 0


def estimate_read_bases_bam(bam_path: Path) -> int:
    """Estimate total bases in a BAM file using samtools.

    Args:
        bam_path: Path to BAM or CRAM file

    Returns:
        Estimated total bases
    """
    if not have_exe("samtools"):
        print("[warn] samtools not found, cannot estimate BAM bases", file=sys.stderr)
        return 0

    try:
        # Use samtools stats for accurate base count
        result = subprocess.run(
            ["samtools", "stats", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=600,  # 10 min timeout
        )
        if result.returncode == 0:
            for line in result.stdout.split("\n"):
                # Look for "SN	bases mapped:	NNNN" or "SN	total length:	NNNN"
                if line.startswith("SN\ttotal length:"):
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        return int(parts[2])
                # For unaligned BAM, look at raw total bases
                if line.startswith("SN\tbases mapped and QC-failed:"):
                    # Skip this one
                    continue

        # Fallback: count reads and estimate
        count_result = subprocess.run(
            ["samtools", "view", "-c", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if count_result.returncode == 0:
            read_count = int(count_result.stdout.strip())
            # Assume average read length (HiFi ~15kb, ONT ~10kb, SR ~150bp)
            # Use conservative estimate of 10kb for long reads
            avg_len = 10000
            return read_count * avg_len

    except subprocess.TimeoutExpired:
        print("[warn] samtools timed out estimating BAM bases", file=sys.stderr)
    except Exception as e:
        print(f"[warn] Could not estimate BAM bases: {e}", file=sys.stderr)

    return 0


def estimate_read_bases(reads_path: Path, read_format: ReadFormat) -> int:
    """Estimate total bases in a reads file.

    Args:
        reads_path: Path to reads file
        read_format: Detected format of the reads

    Returns:
        Estimated total bases
    """
    if read_format in (ReadFormat.FASTQ, ReadFormat.FASTQ_GZ):
        return estimate_read_bases_fastq(reads_path)
    elif read_format in (ReadFormat.UNALIGNED_BAM, ReadFormat.UNALIGNED_CRAM):
        return estimate_read_bases_bam(reads_path)
    else:
        # Aligned BAM/CRAM - shouldn't need downsampling
        return 0


def calculate_subsample_fraction(
    total_read_bases: int,
    assembly_size: int,
    target_coverage: float,
) -> float:
    """Calculate fraction to subsample reads to target coverage.

    Args:
        total_read_bases: Total bases in input reads
        assembly_size: Total assembly size in bp
        target_coverage: Target coverage (e.g., 20.0 for 20X)

    Returns:
        Fraction to subsample (0.0-1.0), or 1.0 if no subsampling needed
    """
    if total_read_bases <= 0 or assembly_size <= 0:
        return 1.0

    target_bases = target_coverage * assembly_size
    fraction = target_bases / total_read_bases

    return min(1.0, fraction)


def downsample_fastq(
    reads_path: Path,
    output_path: Path,
    fraction: float,
    seed: int = 42,
) -> bool:
    """Downsample a FASTQ file using seqtk.

    Args:
        reads_path: Path to input FASTQ or FASTQ.gz
        output_path: Path for output (will be gzipped)
        fraction: Fraction to retain (0.0-1.0)
        seed: Random seed for reproducibility

    Returns:
        True if successful, False otherwise
    """
    if not have_exe("seqtk"):
        print("[error] seqtk not found, cannot downsample FASTQ", file=sys.stderr)
        return False

    try:
        # seqtk sample -s seed input.fq fraction | gzip > output.fq.gz
        print(f"[info] Downsampling FASTQ to {fraction:.2%} with seqtk", file=sys.stderr)

        with open(output_path, "wb") as out_fh:
            seqtk_proc = subprocess.Popen(
                ["seqtk", "sample", "-s", str(seed), str(reads_path), str(fraction)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            gzip_proc = subprocess.Popen(
                ["gzip", "-c"],
                stdin=seqtk_proc.stdout,
                stdout=out_fh,
                stderr=subprocess.PIPE,
            )
            if seqtk_proc.stdout:
                seqtk_proc.stdout.close()

            gzip_proc.wait()
            seqtk_proc.wait()

            if seqtk_proc.returncode != 0 or gzip_proc.returncode != 0:
                print("[error] seqtk/gzip pipeline failed", file=sys.stderr)
                return False

        print(f"[done] Downsampled FASTQ: {output_path}", file=sys.stderr)
        return True

    except Exception as e:
        print(f"[error] FASTQ downsampling failed: {e}", file=sys.stderr)
        return False


def downsample_bam(
    bam_path: Path,
    output_path: Path,
    fraction: float,
    seed: int = 42,
    threads: int = 4,
) -> bool:
    """Downsample a BAM file using samtools.

    Args:
        bam_path: Path to input BAM or CRAM
        output_path: Path for output BAM
        fraction: Fraction to retain (0.0-1.0)
        seed: Random seed for reproducibility
        threads: Number of threads

    Returns:
        True if successful, False otherwise
    """
    if not have_exe("samtools"):
        print("[error] samtools not found, cannot downsample BAM", file=sys.stderr)
        return False

    try:
        # samtools view -s seed.fraction -b input.bam > output.bam
        # The -s option takes seed.fraction format (e.g., 42.25 for seed=42, frac=0.25)
        subsample_arg = f"{seed}.{int(fraction * 100):02d}"

        print(f"[info] Downsampling BAM to {fraction:.2%} with samtools", file=sys.stderr)

        cmd = [
            "samtools", "view",
            "-s", subsample_arg,
            "-@", str(threads),
            "-b",
            "-o", str(output_path),
            str(bam_path),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[error] samtools view failed: {result.stderr}", file=sys.stderr)
            return False

        print(f"[done] Downsampled BAM: {output_path}", file=sys.stderr)
        return True

    except Exception as e:
        print(f"[error] BAM downsampling failed: {e}", file=sys.stderr)
        return False


def downsample_reads(
    reads_path: Path,
    read_format: ReadFormat,
    output_path: Path,
    fraction: float,
    seed: int = 42,
    threads: int = 4,
) -> Tuple[bool, Path]:
    """Downsample reads to a target fraction.

    Args:
        reads_path: Path to input reads
        read_format: Detected format of the reads
        output_path: Path for output file
        fraction: Fraction to retain (0.0-1.0)
        seed: Random seed for reproducibility
        threads: Number of threads for BAM operations

    Returns:
        Tuple of (success, output_path). If fraction >= 1.0, returns original path.
    """
    if fraction >= 1.0:
        print("[info] No downsampling needed (fraction >= 1.0)", file=sys.stderr)
        return True, reads_path

    if read_format in (ReadFormat.FASTQ, ReadFormat.FASTQ_GZ):
        # Output as gzipped FASTQ
        if not str(output_path).endswith(".gz"):
            output_path = Path(str(output_path) + ".gz")
        success = downsample_fastq(reads_path, output_path, fraction, seed)
        return success, output_path if success else reads_path

    elif read_format in (ReadFormat.UNALIGNED_BAM, ReadFormat.UNALIGNED_CRAM):
        # Output as BAM
        if not str(output_path).endswith(".bam"):
            output_path = Path(str(output_path) + ".bam")
        success = downsample_bam(reads_path, output_path, fraction, seed, threads)
        return success, output_path if success else reads_path

    else:
        # Aligned BAM/CRAM - shouldn't downsample
        print("[warn] Cannot downsample aligned BAM/CRAM", file=sys.stderr)
        return True, reads_path


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
        reads_type: Read type for preset selection ("hifi_onthq", "ont", "sr")
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
    reads_type: str = "hifi_onthq",
    window_size: int = 1000,
    target_coverage: Optional[float] = 20.0,
) -> Dict[str, DepthStats]:
    """Calculate read depth metrics for an assembly.

    This is the main entry point for depth analysis. It handles:
    1. Format detection for input reads
    2. Downsampling to target coverage (if alignment needed)
    3. Alignment if reads are not pre-aligned
    4. Running mosdepth for depth calculation
    5. Parsing results into per-contig DepthStats

    Args:
        reads: Path to reads file (FASTQ, BAM, CRAM)
        assembly: Path to assembly FASTA
        contig_lengths: Dictionary of contig names to lengths
        work_dir: Working directory for intermediate files
        threads: Number of threads to use
        reads_type: Read type for minimap2 preset ("hifi_onthq", "ont", "sr")
        window_size: Window size for mosdepth (default 1000bp)
        target_coverage: Target coverage for downsampling before alignment.
            Set to None or 0 to disable downsampling. Default 20.0 (20X).
            Downsampling uses seqtk (FASTQ) or samtools (BAM/CRAM).

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
        # Already aligned - use directly (no downsampling for pre-aligned)
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
        # Need to align - consider downsampling first
        reads_to_align = reads

        # Downsample if target coverage is specified
        if target_coverage and target_coverage > 0:
            assembly_size = sum(contig_lengths.values())
            print(f"[info] Assembly size: {assembly_size:,} bp", file=sys.stderr)

            # Estimate total read bases
            print("[info] Estimating total read bases for downsampling", file=sys.stderr)
            total_read_bases = estimate_read_bases(reads, read_format)

            if total_read_bases > 0:
                current_coverage = total_read_bases / assembly_size
                print(f"[info] Estimated read coverage: {current_coverage:.1f}X", file=sys.stderr)

                fraction = calculate_subsample_fraction(
                    total_read_bases=total_read_bases,
                    assembly_size=assembly_size,
                    target_coverage=target_coverage,
                )

                if fraction < 1.0:
                    print(f"[info] Downsampling to ~{target_coverage:.0f}X coverage (fraction={fraction:.3f})", file=sys.stderr)
                    downsampled_path = work_dir / "downsampled_reads"
                    success, reads_to_align = downsample_reads(
                        reads_path=reads,
                        read_format=read_format,
                        output_path=downsampled_path,
                        fraction=fraction,
                        seed=42,
                        threads=threads,
                    )
                    if not success:
                        print("[warn] Downsampling failed, using full reads", file=sys.stderr)
                        reads_to_align = reads
                else:
                    print(f"[info] Read coverage ({current_coverage:.1f}X) below target ({target_coverage:.0f}X), no downsampling needed", file=sys.stderr)
            else:
                print("[warn] Could not estimate read bases, skipping downsampling", file=sys.stderr)

        bam_path = work_dir / "aligned_reads.sorted.bam"
        err_path = work_dir / "alignment.err"

        success = align_reads_to_assembly(
            reads=reads_to_align,
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
