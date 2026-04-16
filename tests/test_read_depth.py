#!/usr/bin/env python3
"""Unit tests for read depth analysis module."""
import gzip
import textwrap
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnadis.analysis.read_depth import (
    ReadFormat,
    detect_read_format,
    get_minimap2_preset,
    parse_mosdepth_regions,
)
from dnadis.models import DepthStats


def test_detect_read_format_fastq(tmp_path):
    """Test detection of FASTQ format."""
    # Plain FASTQ
    fq = tmp_path / "reads.fq"
    fq.write_text("@read1\nACGT\n+\nIIII\n")
    assert detect_read_format(fq) == ReadFormat.FASTQ

    fq2 = tmp_path / "reads.fastq"
    fq2.write_text("@read1\nACGT\n+\nIIII\n")
    assert detect_read_format(fq2) == ReadFormat.FASTQ


def test_detect_read_format_fastq_gz(tmp_path):
    """Test detection of gzipped FASTQ format."""
    # Gzipped FASTQ
    fq_gz = tmp_path / "reads.fq.gz"
    with gzip.open(fq_gz, "wt") as fh:
        fh.write("@read1\nACGT\n+\nIIII\n")
    assert detect_read_format(fq_gz) == ReadFormat.FASTQ_GZ

    fq_gz2 = tmp_path / "reads.fastq.gz"
    with gzip.open(fq_gz2, "wt") as fh:
        fh.write("@read1\nACGT\n+\nIIII\n")
    assert detect_read_format(fq_gz2) == ReadFormat.FASTQ_GZ


def test_get_minimap2_preset():
    """Test minimap2 preset mapping for different read types."""
    assert get_minimap2_preset("lrhq") == "lr:hqae"
    assert get_minimap2_preset("r9") == "map-ont"
    assert get_minimap2_preset("sr") == "sr"
    # Unknown type should raise ValueError
    with pytest.raises(ValueError, match="Invalid reads_type"):
        get_minimap2_preset("unknown")


def test_parse_mosdepth_regions(tmp_path):
    """Test parsing of mosdepth regions.bed.gz output."""
    regions_file = tmp_path / "depth.regions.bed.gz"

    # Create mock mosdepth output
    # Format: chrom, start, end, depth
    regions_data = textwrap.dedent("""
        contigA\t0\t1000\t50.0
        contigA\t1000\t2000\t55.0
        contigA\t2000\t3000\t45.0
        contigB\t0\t1000\t5.0
        contigB\t1000\t2000\t0.5
    """).strip()

    with gzip.open(regions_file, "wt") as fh:
        fh.write(regions_data + "\n")

    contig_lengths = {"contigA": 3000, "contigB": 2000}

    result = parse_mosdepth_regions(regions_file, contig_lengths)

    assert "contigA" in result
    assert "contigB" in result

    # Check contigA stats
    stats_a = result["contigA"]
    assert isinstance(stats_a, DepthStats)
    assert stats_a.mean_depth == 50.0  # (50*1000 + 55*1000 + 45*1000) / 3000 = 50
    assert stats_a.median_depth == 50.0
    assert stats_a.min_depth == 45.0
    assert stats_a.max_depth == 55.0
    assert stats_a.breadth_1x == 1.0  # All windows >= 1x
    assert stats_a.breadth_10x == 1.0  # All windows >= 10x

    # Check contigB stats (partial coverage)
    stats_b = result["contigB"]
    assert stats_b.min_depth == 0.5
    assert stats_b.max_depth == 5.0
    assert stats_b.breadth_1x == 0.5  # Only first window >= 1x
    assert stats_b.breadth_10x == 0.0  # No windows >= 10x


def test_parse_mosdepth_regions_empty_contig(tmp_path):
    """Test handling of contigs with no coverage data."""
    regions_file = tmp_path / "depth.regions.bed.gz"

    # Only contigA has data
    regions_data = "contigA\t0\t1000\t30.0\n"
    with gzip.open(regions_file, "wt") as fh:
        fh.write(regions_data)

    contig_lengths = {"contigA": 1000, "contigB": 500}

    result = parse_mosdepth_regions(regions_file, contig_lengths)

    # contigB should have zero stats
    stats_b = result["contigB"]
    assert stats_b.mean_depth == 0.0
    assert stats_b.median_depth == 0.0
    assert stats_b.breadth_1x == 0.0
    assert stats_b.breadth_10x == 0.0


def test_depth_stats_dataclass():
    """Test DepthStats dataclass creation and attributes."""
    stats = DepthStats(
        mean_depth=45.5,
        median_depth=42.0,
        std_depth=5.2,
        min_depth=10.0,
        max_depth=80.0,
        breadth_1x=0.95,
        breadth_10x=0.85,
    )

    assert stats.mean_depth == 45.5
    assert stats.median_depth == 42.0
    assert stats.std_depth == 5.2
    assert stats.min_depth == 10.0
    assert stats.max_depth == 80.0
    assert stats.breadth_1x == 0.95
    assert stats.breadth_10x == 0.85


def test_parse_mosdepth_varying_window_sizes(tmp_path):
    """Test parsing with varying window sizes (e.g., at contig ends)."""
    regions_file = tmp_path / "depth.regions.bed.gz"

    # Last window is smaller (contig doesn't divide evenly)
    regions_data = textwrap.dedent("""
        contigA\t0\t1000\t20.0
        contigA\t1000\t2000\t30.0
        contigA\t2000\t2500\t40.0
    """).strip()

    with gzip.open(regions_file, "wt") as fh:
        fh.write(regions_data + "\n")

    contig_lengths = {"contigA": 2500}

    result = parse_mosdepth_regions(regions_file, contig_lengths)

    stats = result["contigA"]
    # Weighted mean: (20*1000 + 30*1000 + 40*500) / 2500 = 28.0
    assert abs(stats.mean_depth - 28.0) < 0.01
    assert stats.min_depth == 20.0
    assert stats.max_depth == 40.0
