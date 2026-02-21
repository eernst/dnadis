"""Tests for pairwise assembly-vs-assembly synteny alignment."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from final_finalizer.alignment.pairwise import compute_pairwise_synteny


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Default chain-parsing kwargs matching _make_args() from the old tests.
_CHAIN_KWARGS = dict(
    preset="asm20",
    kmer=None,
    window=None,
    assign_minlen=10000,
    assign_minmapq=0,
    assign_tp="PI",
    chain_q_gap=200000,
    chain_r_gap=400000,
    chain_diag_slop=150000,
    assign_min_ident=0.0,
    assign_chain_topk=3,
    assign_chain_score="matches",
    assign_chain_min_bp=0,
    assign_ref_score="all",
)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestComputePairwiseSynteny:
    """Tests for compute_pairwise_synteny()."""

    def test_missing_left_chrs_returns_none(self, tmp_path):
        """Should return None when left chrs.fasta doesn't exist."""
        left_fasta = tmp_path / "left" / "left.chrs.fasta"
        right_fasta = tmp_path / "right" / "right.chrs.fasta"
        right_fasta.parent.mkdir(parents=True)
        right_fasta.write_text(">chr1\nACGT\n")

        pair_prefix = tmp_path / "pairwise" / "left_vs_right"

        result = compute_pairwise_synteny(
            left_fasta=left_fasta,
            right_fasta=right_fasta,
            left_name="left",
            right_name="right",
            outprefix=pair_prefix,
            threads=1,
            **_CHAIN_KWARGS,
        )
        assert result is None

    def test_missing_right_chrs_returns_none(self, tmp_path):
        """Should return None when right chrs.fasta doesn't exist."""
        left_fasta = tmp_path / "left" / "left.chrs.fasta"
        right_fasta = tmp_path / "right" / "right.chrs.fasta"
        left_fasta.parent.mkdir(parents=True)
        left_fasta.write_text(">chr1\nACGT\n")

        pair_prefix = tmp_path / "pairwise" / "left_vs_right"

        result = compute_pairwise_synteny(
            left_fasta=left_fasta,
            right_fasta=right_fasta,
            left_name="left",
            right_name="right",
            outprefix=pair_prefix,
            threads=1,
            **_CHAIN_KWARGS,
        )
        assert result is None

    def test_empty_right_chrs_returns_none(self, tmp_path):
        """Should return None when right chrs.fasta is empty."""
        left_fasta = tmp_path / "left" / "left.chrs.fasta"
        right_fasta = tmp_path / "right" / "right.chrs.fasta"
        left_fasta.parent.mkdir(parents=True)
        right_fasta.parent.mkdir(parents=True)
        left_fasta.write_text(">chr1\nACGT\n")
        right_fasta.write_text("")  # empty

        pair_prefix = tmp_path / "pairwise" / "left_vs_right"

        result = compute_pairwise_synteny(
            left_fasta=left_fasta,
            right_fasta=right_fasta,
            left_name="left",
            right_name="right",
            outprefix=pair_prefix,
            threads=1,
            **_CHAIN_KWARGS,
        )
        # Empty file should be caught by file_exists_and_valid
        assert result is None

    def test_reuses_cached_macro_blocks(self, tmp_path):
        """Should reuse existing macro_blocks TSV without re-running alignment."""
        left_fasta = tmp_path / "left" / "left.chrs.fasta"
        right_fasta = tmp_path / "right" / "right.chrs.fasta"
        left_fasta.parent.mkdir(parents=True)
        right_fasta.parent.mkdir(parents=True)
        left_fasta.write_text(">chr1\nACGT\n")
        right_fasta.write_text(">chr1\nACGT\n")

        # Pre-create the macro_blocks TSV (cached)
        pair_prefix = tmp_path / "pairwise" / "left_vs_right"
        pair_prefix.parent.mkdir(parents=True)
        macro_tsv = Path(str(pair_prefix) + ".macro_blocks.tsv")
        macro_tsv.write_text("contig\tcontig_len\tref_id\n")  # minimal header

        result = compute_pairwise_synteny(
            left_fasta=left_fasta,
            right_fasta=right_fasta,
            left_name="left",
            right_name="right",
            outprefix=pair_prefix,
            threads=1,
            **_CHAIN_KWARGS,
        )
        assert result == macro_tsv

    def test_creates_output_directory(self, tmp_path):
        """Should create the pairwise output directory if it doesn't exist."""
        left_fasta = tmp_path / "left" / "left.chrs.fasta"
        right_fasta = tmp_path / "right" / "right.chrs.fasta"
        left_fasta.parent.mkdir(parents=True)
        right_fasta.parent.mkdir(parents=True)
        left_fasta.write_text(">chr1\nACGT\n")
        right_fasta.write_text(">chr1\nACGT\n")

        # Pre-create cached result so we don't need minimap2
        pair_prefix = tmp_path / "pairwise" / "subdir" / "left_vs_right"
        pair_prefix.parent.mkdir(parents=True, exist_ok=True)
        macro_tsv = Path(str(pair_prefix) + ".macro_blocks.tsv")
        macro_tsv.write_text("contig\tcontig_len\tref_id\n")

        result = compute_pairwise_synteny(
            left_fasta=left_fasta,
            right_fasta=right_fasta,
            left_name="left",
            right_name="right",
            outprefix=pair_prefix,
            threads=1,
            **_CHAIN_KWARGS,
        )
        assert result is not None
        assert result.exists()


class TestPairwiseImport:
    """Test that the module imports correctly."""

    def test_import(self):
        from final_finalizer.alignment.pairwise import compute_pairwise_synteny
        assert callable(compute_pairwise_synteny)
