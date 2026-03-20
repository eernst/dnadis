"""Resource estimation for distributed job submission.

Each ``estimate_*`` function returns a :class:`ResourceSpec` with cores,
memory, and wall-time sized for the given input files.  All estimates are
deliberately generous — under-requesting causes SLURM OOM kills which are
far more expensive than over-requesting.

Estimates are capped by :func:`~final_finalizer.utils.distributed.clamp_resources`
before submission.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Optional

from final_finalizer.utils.distributed import ClusterConfig, ResourceSpec, clamp_resources


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _file_size_bytes(path: Path) -> int:
    """Return file size in bytes, or 0 if the file doesn't exist."""
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _estimate_genome_bp_from_filesize(path: Path) -> int:
    """Estimate genome size in base pairs from a FASTA file size.

    For gzipped files, assumes ~3× compression ratio.
    For plain text FASTA, assumes ~50% of file size is sequence.
    """
    size = _file_size_bytes(path)
    if size == 0:
        return 0

    name = path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        # gzip typically achieves ~3x on FASTA; sequence is ~50% of decompressed
        return int(size * 3 * 0.5)
    else:
        # Plain FASTA: roughly 50% is sequence after removing headers/newlines
        return int(size * 0.5)


# Global safety multiplier for wall-time estimates.  Under-requesting causes
# SLURM TIMEOUT kills which hang the coordinator (executorlib futures never
# resolve), so we deliberately over-request.  Over-requesting only wastes
# scheduler priority, not actual compute time — jobs exit as soon as they finish.
_TIME_SAFETY_FACTOR = 3


def _scale_time(base_minutes: int, genome_bp: int, scale_bp: int = 500_000_000) -> int:
    """Scale wall-time linearly with genome size, then apply safety factor.

    *base_minutes* is for a genome of *scale_bp* bases.
    Minimum returned is *base_minutes* × ``_TIME_SAFETY_FACTOR``.
    """
    if genome_bp <= 0:
        return base_minutes * _TIME_SAFETY_FACTOR
    factor = max(1.0, genome_bp / scale_bp)
    return max(base_minutes, int(math.ceil(base_minutes * factor * _TIME_SAFETY_FACTOR)))


# ---------------------------------------------------------------------------
# Per-phase estimators
# ---------------------------------------------------------------------------
def estimate_synteny_resources(
    ref_path: Path,
    qry_path: Path,
    synteny_mode: str,
    config: ClusterConfig,
) -> ResourceSpec:
    """Estimate resources for synteny alignment (phase 2)."""
    ref_bp = _estimate_genome_bp_from_filesize(ref_path)
    qry_bp = _estimate_genome_bp_from_filesize(qry_path)

    if synteny_mode == "protein":
        # miniprot: memory dominated by query index
        mem_gb = max(2.0, qry_bp * 6 / 1e9 + 1.0)
        time_min = _scale_time(30, qry_bp)
    else:
        # minimap2 whole-genome: memory dominated by reference index
        mem_gb = max(4.0, ref_bp * 8 / 1e9 + 2.0)
        time_min = _scale_time(60, max(ref_bp, qry_bp))

    cores = min(32, config.max_threads)

    spec = ResourceSpec(
        cores=cores,
        memory_gb=mem_gb,
        time_minutes=time_min,
        job_name="synteny",
    )
    return clamp_resources(spec, config)


def estimate_blast_resources(
    query_path: Path,
    config: ClusterConfig,
    job_name: str = "blast",
) -> ResourceSpec:
    """Estimate resources for BLAST-based detection (organelle, rDNA)."""
    qry_bp = _estimate_genome_bp_from_filesize(query_path)
    time_min = _scale_time(15, qry_bp, scale_bp=200_000_000)

    spec = ResourceSpec(
        cores=min(8, config.max_threads),
        memory_gb=2.0,
        time_minutes=time_min,
        job_name=job_name,
    )
    return clamp_resources(spec, config)


def estimate_debris_resources(
    query_path: Path,
    config: ClusterConfig,
) -> ResourceSpec:
    """Estimate resources for chromosome debris detection (phase 5)."""
    qry_bp = _estimate_genome_bp_from_filesize(query_path)
    mem_gb = max(2.0, qry_bp * 8 / 1e9 + 1.0)
    time_min = _scale_time(30, qry_bp)

    spec = ResourceSpec(
        cores=min(8, config.max_threads),
        memory_gb=mem_gb,
        time_minutes=time_min,
        job_name="debris",
    )
    return clamp_resources(spec, config)


def estimate_contaminant_resources(
    centrifuger_idx: str,
    config: ClusterConfig,
) -> ResourceSpec:
    """Estimate resources for contaminant detection (phase 6)."""
    # Index size drives memory; try .1.cfr first, fall back to prefix.cfr
    idx_path = Path(centrifuger_idx + ".1.cfr")
    if not idx_path.exists():
        idx_path = Path(centrifuger_idx + ".cfr")
    idx_bytes = _file_size_bytes(idx_path)
    mem_gb = max(4.0, idx_bytes * 1.2 / 1e9 + 2.0)

    spec = ResourceSpec(
        cores=min(16, config.max_threads),
        memory_gb=mem_gb,
        time_minutes=30 * _TIME_SAFETY_FACTOR,
        job_name="contaminant",
    )
    return clamp_resources(spec, config)


def estimate_pairwise_resources(
    left_fasta: Path,
    right_fasta: Path,
    config: ClusterConfig,
) -> ResourceSpec:
    """Estimate resources for pairwise assembly-vs-assembly synteny alignment."""
    left_bp = _estimate_genome_bp_from_filesize(left_fasta)
    right_bp = _estimate_genome_bp_from_filesize(right_fasta)
    max_bp = max(left_bp, right_bp)

    # minimap2 whole-genome: memory dominated by target (left) index
    mem_gb = max(4.0, max_bp * 8 / 1e9 + 2.0)
    time_min = _scale_time(30, max_bp)
    cores = min(32, config.max_threads)

    spec = ResourceSpec(
        cores=cores,
        memory_gb=mem_gb,
        time_minutes=time_min,
        job_name="pairwise",
    )
    return clamp_resources(spec, config)


def estimate_compleasm_resources(
    fasta_path: Path,
    config: ClusterConfig,
) -> ResourceSpec:
    """Estimate resources for compleasm BUSCO evaluation (phase 17)."""
    qry_bp = _estimate_genome_bp_from_filesize(fasta_path)
    # compleasm/miniprot: memory scales modestly with genome size
    mem_gb = max(4.0, qry_bp * 4 / 1e9 + 2.0)
    time_min = _scale_time(30, qry_bp)

    spec = ResourceSpec(
        cores=min(16, config.max_threads),
        memory_gb=mem_gb,
        time_minutes=time_min,
        job_name="compleasm",
    )
    return clamp_resources(spec, config)


def estimate_depth_resources(
    reads_path: Path,
    assembly_path: Path,
    config: ClusterConfig,
) -> ResourceSpec:
    """Estimate resources for read depth analysis (phase 12)."""
    asm_bp = _estimate_genome_bp_from_filesize(assembly_path)
    reads_bytes = _file_size_bytes(reads_path)

    # Memory: index + sort buffer (8 GB sort budget to reduce NFS spill)
    mem_gb = max(32.0, asm_bp * 8 / 1e9 + 8.0)
    # Time scales with reads file size (alignment dominates)
    time_min = _scale_time(60, reads_bytes, scale_bp=2_000_000_000)

    spec = ResourceSpec(
        cores=min(32, config.max_threads),
        memory_gb=mem_gb,
        time_minutes=time_min,
        job_name="depth",
    )
    return clamp_resources(spec, config)
