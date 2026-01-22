#!/usr/bin/env python3
"""
Read depth analysis for final_finalizer.

Provides functions for:
- Detecting read format (FASTQ, unaligned BAM, aligned BAM/CRAM)
- Downsampling reads to target coverage (rasusa for FASTQ, samtools for BAM)
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
import json
import statistics
import subprocess
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Tuple

from final_finalizer.models import DepthStats
from final_finalizer.utils.io_utils import file_exists_and_valid, have_exe
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("read_depth")


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
        logger.warning("samtools not found; assuming BAM is unaligned")
        return False

    try:
        result = subprocess.run(
            ["samtools", "view", "-H", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode != 0:
            logger.warning(f"samtools view -H failed for {bam_path}; assuming unaligned")
            return False

        # Check for @SQ lines (reference sequences) - indicates aligned BAM
        for line in result.stdout.split("\n"):
            if line.startswith("@SQ"):
                return True
        return False

    except subprocess.TimeoutExpired:
        logger.warning(f"samtools timeout checking {bam_path}; assuming unaligned")
        return False
    except Exception as e:
        logger.warning(f"Error checking BAM alignment status: {e}; assuming unaligned")
        return False


VALID_READ_TYPES = ("hifi_onthq", "ont", "sr")


def get_minimap2_preset(reads_type: str) -> str:
    """Get minimap2 preset for read type.

    Args:
        reads_type: One of "hifi_onthq", "ont", "sr"

    Returns:
        Minimap2 -x preset string

    Raises:
        ValueError: If reads_type is not a valid option

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
    if reads_type not in presets:
        raise ValueError(
            f"Invalid reads_type: '{reads_type}'. Must be one of: {', '.join(VALID_READ_TYPES)}"
        )
    return presets[reads_type]


def estimate_read_bases_bam(bam_path: Path) -> int:
    """Estimate total bases in a BAM file using samtools.

    Args:
        bam_path: Path to BAM or CRAM file

    Returns:
        Estimated total bases
    """
    if not have_exe("samtools"):
        logger.warning("samtools not found, cannot estimate BAM bases")
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
        logger.warning("samtools timed out estimating BAM bases")
    except Exception as e:
        logger.warning(f"Could not estimate BAM bases: {e}")

    return 0


def estimate_read_bases(reads_path: Path, read_format: ReadFormat) -> int:
    """Estimate total bases in a reads file (BAM/CRAM only).

    Args:
        reads_path: Path to reads file
        read_format: Detected format of the reads

    Returns:
        Estimated total bases
    """
    if read_format in (ReadFormat.UNALIGNED_BAM, ReadFormat.UNALIGNED_CRAM):
        return estimate_read_bases_bam(reads_path)
    else:
        # FASTQ and aligned BAM/CRAM are handled elsewhere.
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
    target_coverage: float,
    genome_size: int,
    seed: int = 42,
) -> bool:
    """Downsample a FASTQ file using rasusa.

    Args:
        reads_path: Path to input FASTQ or FASTQ.gz
        output_path: Path for output (will be gzipped)
        target_coverage: Target coverage to retain (e.g., 20.0 for 20X)
        genome_size: Genome size in bp (used to compute target bases)
        seed: Random seed for reproducibility

    Returns:
        True if successful, False otherwise
    """
    if not have_exe("rasusa"):
        logger.error("rasusa not found, cannot downsample FASTQ")
        return False

    if genome_size <= 0 or target_coverage <= 0:
        logger.error("Invalid genome size or target coverage for downsampling")
        return False

    try:
        logger.info(f"Downsampling FASTQ to ~{target_coverage:.2f}X with rasusa")

        cmd = [
            "rasusa",
            "reads",
            "--coverage", f"{target_coverage:.6g}",
            "--genome-size", str(genome_size),
            "--seed", str(seed),
            "-o", str(output_path),
            str(reads_path)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            stderr_output = (result.stderr or result.stdout or "").strip()
            logger.error(f"rasusa failed (code {result.returncode}): {stderr_output}")
            return False

        logger.done(f"Downsampled FASTQ: {output_path}")
        return True

    except Exception as e:
        logger.error(f"FASTQ downsampling failed: {e}")
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
        logger.error("samtools not found, cannot downsample BAM")
        return False

    # Validate fraction is in valid range
    if fraction <= 0.0 or fraction > 1.0:
        logger.error(f"Invalid subsample fraction: {fraction}. Must be in (0.0, 1.0]")
        return False

    try:
        # samtools view -s SEED.FRAC where integer part is seed, decimal part is fraction
        # e.g., 42.25 means seed=42, keep 25% of reads
        # Construct by adding seed + fraction to get correct format
        subsample_value = seed + fraction
        subsample_arg = f"{subsample_value:.6f}"  # Preserve precision for small fractions

        logger.info(f"Downsampling BAM to {fraction:.2%} with samtools (-s {subsample_arg})")

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
            logger.error(f"samtools view -s {subsample_arg} failed: {result.stderr}")
            return False

        logger.done(f"Downsampled BAM: {output_path}")
        return True

    except Exception as e:
        logger.error(f"BAM downsampling failed: {e}")
        return False


def downsample_reads(
    reads_path: Path,
    read_format: ReadFormat,
    output_path: Path,
    fraction: float,
    seed: int = 42,
    threads: int = 4,
    target_coverage: Optional[float] = None,
    genome_size: Optional[int] = None,
) -> Tuple[bool, Path]:
    """Downsample reads to a target fraction/coverage.

    Args:
        reads_path: Path to input reads
        read_format: Detected format of the reads
        output_path: Path for output file
        fraction: Fraction to retain (0.0-1.0)
        seed: Random seed for reproducibility
        threads: Number of threads for BAM operations
        target_coverage: Target coverage for FASTQ downsampling (rasusa)
        genome_size: Genome size in bp for FASTQ downsampling (rasusa)

    Returns:
        Tuple of (success, output_path). For FASTQ, target_coverage/genome_size are used.
    """

    if read_format in (ReadFormat.FASTQ, ReadFormat.FASTQ_GZ):
        # Output as gzipped FASTQ
        if not str(output_path).endswith(".gz"):
            output_path = Path(str(output_path) + ".gz")
        if target_coverage is None or genome_size is None:
            logger.error("Missing target coverage or genome size for FASTQ downsampling")
            return False, reads_path
        success = downsample_fastq(reads_path, output_path, target_coverage, genome_size, seed)
        return success, output_path if success else reads_path

    elif read_format in (ReadFormat.UNALIGNED_BAM, ReadFormat.UNALIGNED_CRAM):
        if fraction >= 1.0:
            logger.info("No downsampling needed (fraction >= 1.0)")
            return True, reads_path
        # Output as BAM
        if not str(output_path).endswith(".bam"):
            output_path = Path(str(output_path) + ".bam")
        success = downsample_bam(reads_path, output_path, fraction, seed, threads)
        return success, output_path if success else reads_path

    else:
        # Aligned BAM/CRAM - shouldn't downsample
        logger.warning("Cannot downsample aligned BAM/CRAM")
        return True, reads_path


def write_alignment_metadata(
    metadata_path: Path,
    original_reads: Path,
    assembly: Path,
    reads_type: str,
    target_coverage: Optional[float],
    downsampled: bool,
    window_size: int = 1000,
) -> None:
    """Write alignment metadata for cache validation.

    Args:
        metadata_path: Path to write metadata JSON
        original_reads: Path to original reads file (before downsampling)
        assembly: Path to assembly FASTA
        reads_type: Read type used for alignment
        target_coverage: Target coverage for downsampling (None if disabled)
        downsampled: Whether reads were downsampled
        window_size: Window size for mosdepth (for cache validation)
    """
    metadata = {
        "original_reads": str(original_reads.resolve()),
        "original_reads_name": original_reads.name,
        "assembly": str(assembly.resolve()),
        "assembly_name": assembly.name,
        "reads_type": reads_type,
        "target_coverage": target_coverage,
        "downsampled": downsampled,
        "window_size": window_size,
    }
    try:
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)
    except Exception as e:
        logger.warning(f"Could not write alignment metadata: {e}")


def check_alignment_cache(
    bam_path: Path,
    metadata_path: Path,
    original_reads: Path,
    assembly: Path,
    reads_type: str,
    target_coverage: Optional[float],
) -> bool:
    """Check if cached alignment is valid for current parameters.

    Validates that:
    1. BAM file exists and is non-empty
    2. Metadata file exists and matches current parameters
    3. BAM index exists

    Args:
        bam_path: Path to aligned BAM
        metadata_path: Path to metadata JSON
        original_reads: Current original reads path
        assembly: Current assembly path
        reads_type: Current reads type
        target_coverage: Current target coverage

    Returns:
        True if cache is valid and can be reused, False otherwise
    """
    # Check BAM exists
    if not file_exists_and_valid(bam_path, min_size=100):
        return False

    # Check metadata exists
    if not metadata_path.exists():
        logger.info("No alignment metadata found, will realign")
        return False

    # Load and validate metadata
    try:
        with open(metadata_path) as f:
            metadata = json.load(f)

        # Check original reads match (by resolved path)
        cached_reads = metadata.get("original_reads", "")
        current_reads = str(original_reads.resolve())
        if cached_reads != current_reads:
            logger.info(f"Reads file changed: cached={metadata.get('original_reads_name', 'unknown')}, "
                        f"current={original_reads.name}")
            return False

        # Check assembly matches
        cached_assembly = metadata.get("assembly", "")
        current_assembly = str(assembly.resolve())
        if cached_assembly != current_assembly:
            logger.info(f"Assembly changed: cached={metadata.get('assembly_name', 'unknown')}, "
                        f"current={assembly.name}")
            return False

        # Check reads type matches
        cached_reads_type = metadata.get("reads_type", "")
        if cached_reads_type != reads_type:
            logger.info(f"Reads type changed: cached={cached_reads_type}, current={reads_type}")
            return False

        # Check target coverage matches (for downsampling consistency)
        cached_coverage = metadata.get("target_coverage")
        if cached_coverage != target_coverage:
            logger.info(f"Target coverage changed: cached={cached_coverage}, current={target_coverage}")
            return False

    except (json.JSONDecodeError, KeyError) as e:
        logger.warning(f"Invalid alignment metadata: {e}")
        return False

    # Check BAM index exists
    bai_path = Path(str(bam_path) + ".bai")
    csi_path = Path(str(bam_path) + ".csi")
    if not bai_path.exists() and not csi_path.exists():
        logger.info("BAM index missing, will create")
        # Try to create index
        if have_exe("samtools"):
            try:
                subprocess.run(
                    ["samtools", "index", str(bam_path)],
                    check=True,
                    stderr=subprocess.DEVNULL,
                )
                logger.info("Created BAM index")
            except subprocess.CalledProcessError:
                logger.warning("Could not create BAM index, will realign")
                return False
        else:
            return False

    return True


def check_depth_cache(
    regions_bed_gz: Path,
    metadata_path: Path,
    original_reads: Path,
    assembly: Path,
    reads_type: str,
    target_coverage: Optional[float],
    window_size: int,
) -> bool:
    """Check if cached depth output is valid for current parameters.

    Returns True if mosdepth regions file exists AND metadata matches all params.
    This allows skipping both alignment AND mosdepth on re-run.

    Args:
        regions_bed_gz: Path to mosdepth regions.bed.gz output file
        metadata_path: Path to metadata JSON
        original_reads: Current original reads path
        assembly: Current assembly path
        reads_type: Current reads type
        target_coverage: Current target coverage
        window_size: Current mosdepth window size

    Returns:
        True if cache is valid and can be reused, False otherwise
    """
    # Check mosdepth output exists
    if not file_exists_and_valid(regions_bed_gz, min_size=50):
        return False

    # Check metadata exists and matches
    if not metadata_path.exists():
        return False

    try:
        with open(metadata_path) as f:
            metadata = json.load(f)

        # Validate all parameters match
        if metadata.get("original_reads") != str(original_reads.resolve()):
            return False
        if metadata.get("assembly") != str(assembly.resolve()):
            return False
        if metadata.get("reads_type") != reads_type:
            return False
        if metadata.get("target_coverage") != target_coverage:
            return False
        if metadata.get("window_size") != window_size:
            logger.info(f"Window size changed: cached={metadata.get('window_size')}, current={window_size}")
            return False

        return True
    except Exception:
        return False


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
        logger.info(f"Aligned BAM exists, reusing: {output_bam}")
        return True

    work_dir.mkdir(parents=True, exist_ok=True)

    # Check for required tools
    mapper = None
    for exe in ["mm2plus", "minimap2"]:
        if have_exe(exe):
            mapper = exe
            break

    if not mapper:
        logger.error("Neither mm2plus nor minimap2 found in PATH")
        return False

    if not have_exe("samtools"):
        logger.error("samtools not found in PATH")
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

    logger.info(f"Aligning reads with {mapper} -x {preset}")

    # Track processes for cleanup
    p_fastq = None
    p_mm2 = None
    p_sort = None

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

                # Check all return codes with detailed error messages
                if p_fastq.returncode != 0:
                    logger.error(f"samtools fastq failed (code {p_fastq.returncode}, see {err_file})")
                    return False
                if p_mm2.returncode != 0:
                    logger.error(f"{mapper} failed (code {p_mm2.returncode}, see {err_file})")
                    return False
                if p_sort.returncode != 0:
                    logger.error(f"samtools sort failed (code {p_sort.returncode}, see {err_file})")
                    return False
            else:
                # Pipeline: minimap2 | samtools sort
                p_mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=err_fh)
                p_sort = subprocess.Popen(sort_cmd, stdin=p_mm2.stdout, stderr=err_fh)
                if p_mm2.stdout:
                    p_mm2.stdout.close()

                p_sort.wait()
                p_mm2.wait()

                # Check return codes with detailed error messages
                if p_mm2.returncode != 0:
                    logger.error(f"{mapper} failed (code {p_mm2.returncode}, see {err_file})")
                    return False
                if p_sort.returncode != 0:
                    logger.error(f"samtools sort failed (code {p_sort.returncode}, see {err_file})")
                    return False

        # Index the BAM file
        idx_cmd = ["samtools", "index", "-@", str(threads), str(output_bam)]
        subprocess.run(idx_cmd, check=True, stderr=subprocess.DEVNULL)

        logger.done(f"Alignment complete: {output_bam}")
        return True

    except Exception as e:
        logger.error(f"Alignment failed: {e}")
        return False

    finally:
        # Cleanup: terminate any still-running processes
        for proc in [p_fastq, p_mm2, p_sort]:
            if proc is not None and proc.poll() is None:
                try:
                    proc.terminate()
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()  # Must wait after kill to reap zombie
                except Exception:
                    try:
                        proc.kill()
                        proc.wait()
                    except Exception:
                        pass


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
        logger.info(f"mosdepth output exists, reusing: {regions_file}")
        return True

    if not have_exe("mosdepth"):
        logger.error("mosdepth not found in PATH")
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

    logger.info(f"Running mosdepth with {window_size}bp windows")

    try:
        err_file = err_path or Path(str(output_prefix) + ".mosdepth.err")
        err_file.parent.mkdir(parents=True, exist_ok=True)

        with err_file.open("wb") as err_fh:
            result = subprocess.run(cmd, stderr=err_fh, check=False)

        if result.returncode != 0:
            logger.error(f"mosdepth failed (see {err_file})")
            return False

        logger.done(f"mosdepth complete: {regions_file}")
        return True

    except Exception as e:
        logger.error(f"mosdepth failed: {e}")
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
    parse_errors = 0
    max_parse_errors = 100  # Limit error reporting

    try:
        with gzip.open(regions_bed_gz, "rt") as fh:
            for line_num, line in enumerate(fh, 1):
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue

                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    depth = float(parts[3])

                    # Validate interval
                    if start < 0 or end <= start:
                        if parse_errors < max_parse_errors:
                            logger.warning(f"Invalid BED interval at line {line_num}: start={start}, end={end}")
                        parse_errors += 1
                        continue

                    # Validate depth (can be 0, but not negative)
                    if depth < 0:
                        if parse_errors < max_parse_errors:
                            logger.warning(f"Negative depth at line {line_num}: {depth}")
                        parse_errors += 1
                        continue

                    window_size = end - start

                    if chrom in contig_depths:
                        contig_depths[chrom].append(depth)
                        contig_windows[chrom].append((start, end, depth, window_size))

                except (ValueError, IndexError) as e:
                    if parse_errors < max_parse_errors:
                        logger.warning(f"Malformed BED line {line_num}: {e}")
                    parse_errors += 1
                    continue

        if parse_errors > 0:
            logger.warning(f"Total BED parsing errors: {parse_errors}")

    except Exception as e:
        logger.warning(f"Error parsing mosdepth output {regions_bed_gz}: {e}")
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
    keep_bam: bool = False,
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
            Downsampling uses rasusa (FASTQ) or samtools (BAM/CRAM).
        keep_bam: Keep aligned BAM file after depth analysis (default: delete to save space)

    Returns:
        Dictionary mapping contig names to DepthStats objects
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Early check: if depth output already valid, skip everything
    regions_file = work_dir / "depth.regions.bed.gz"
    metadata_path = work_dir / "alignment_metadata.json"

    if check_depth_cache(regions_file, metadata_path, reads, assembly,
                         reads_type, target_coverage, window_size):
        logger.info(f"Reusing cached depth output: {regions_file}")
        return parse_mosdepth_regions(regions_file, contig_lengths)

    # Detect read format
    try:
        read_format = detect_read_format(reads)
        logger.info(f"Detected read format: {read_format.value}")
    except ValueError as e:
        logger.error(f"{e}")
        return {}

    # Determine if we need to align
    if read_format in (ReadFormat.ALIGNED_BAM, ReadFormat.ALIGNED_CRAM):
        # Already aligned - use directly (no downsampling for pre-aligned)
        bam_path = reads
        logger.info("Using pre-aligned BAM/CRAM")

        # Check if index exists, create if needed
        bai_path = Path(str(reads) + ".bai")
        csi_path = Path(str(reads) + ".csi")
        if not bai_path.exists() and not csi_path.exists():
            logger.info("Indexing BAM file")
            if have_exe("samtools"):
                subprocess.run(
                    ["samtools", "index", "-@", str(threads), str(reads)],
                    check=True,
                    stderr=subprocess.DEVNULL,
                )
    else:
        # Need to align - check cache first, then consider downsampling
        bam_path = work_dir / "aligned_reads.sorted.bam"
        metadata_path = work_dir / "alignment_metadata.json"
        err_path = work_dir / "alignment.err"

        # Check if we have a valid cached alignment from the same reads
        cache_valid = check_alignment_cache(
            bam_path=bam_path,
            metadata_path=metadata_path,
            original_reads=reads,
            assembly=assembly,
            reads_type=reads_type,
            target_coverage=target_coverage,
        )

        if cache_valid:
            logger.info(f"Reusing cached alignment: {bam_path}")
        else:
            # Need to perform alignment - consider downsampling first
            reads_to_align = reads
            downsampled = False

            # Downsample if target coverage is specified
            if target_coverage and target_coverage > 0:
                assembly_size = sum(contig_lengths.values())
                logger.info(f"Assembly size: {assembly_size:,} bp")

                if read_format in (ReadFormat.FASTQ, ReadFormat.FASTQ_GZ):
                    downsampled_path = work_dir / "downsampled_reads"
                    success, reads_to_align = downsample_reads(
                        reads_path=reads,
                        read_format=read_format,
                        output_path=downsampled_path,
                        fraction=1.0,
                        seed=42,
                        threads=threads,
                        target_coverage=target_coverage,
                        genome_size=assembly_size,
                    )
                    if not success:
                        logger.warning("Downsampling failed, using full reads")
                        reads_to_align = reads
                    else:
                        downsampled = True
                else:
                    # Estimate total read bases for BAM/CRAM downsampling.
                    logger.info("Estimating total read bases for downsampling")
                    total_read_bases = estimate_read_bases(reads, read_format)

                    if total_read_bases > 0:
                        current_coverage = total_read_bases / assembly_size
                        logger.info(f"Estimated read coverage: {current_coverage:.1f}X")

                        fraction = calculate_subsample_fraction(
                            total_read_bases=total_read_bases,
                            assembly_size=assembly_size,
                            target_coverage=target_coverage,
                        )

                        if fraction < 1.0:
                            logger.info(f"Downsampling to ~{target_coverage:.0f}X coverage (fraction={fraction:.3f})")
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
                                logger.warning("Downsampling failed, using full reads")
                                reads_to_align = reads
                            else:
                                downsampled = True
                        else:
                            logger.info(f"Read coverage ({current_coverage:.1f}X) below target ({target_coverage:.0f}X), no downsampling needed")
                    else:
                        logger.warning("Could not estimate read bases, skipping downsampling")

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
                logger.error("Read alignment failed")
                return {}

            # Write metadata for future cache validation
            write_alignment_metadata(
                metadata_path=metadata_path,
                original_reads=reads,
                assembly=assembly,
                reads_type=reads_type,
                target_coverage=target_coverage,
                downsampled=downsampled,
                window_size=window_size,
            )

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
        logger.error("mosdepth failed")
        return {}

    # Parse mosdepth output
    regions_file = Path(str(mosdepth_prefix) + ".regions.bed.gz")
    if not regions_file.exists():
        logger.error(f"mosdepth output not found: {regions_file}")
        return {}

    depth_stats = parse_mosdepth_regions(regions_file, contig_lengths)
    logger.info(f"Computed depth stats for {len(depth_stats)} contigs")

    # After successful mosdepth, cleanup BAM unless --keep-depth-bam
    # Only cleanup if we created the BAM (not pre-aligned)
    if not keep_bam and read_format not in (ReadFormat.ALIGNED_BAM, ReadFormat.ALIGNED_CRAM):
        bam_to_clean = work_dir / "aligned_reads.sorted.bam"
        if bam_to_clean.exists():
            try:
                bam_to_clean.unlink()
                bai_path = Path(str(bam_to_clean) + ".bai")
                if bai_path.exists():
                    bai_path.unlink()
                logger.info(f"Cleaned up intermediate BAM: {bam_to_clean.name}")
            except Exception as e:
                logger.warning(f"Could not cleanup BAM: {e}")

    return depth_stats
