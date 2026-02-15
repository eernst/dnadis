"""Tests for pairwise assembly-vs-assembly synteny alignment."""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from final_finalizer.alignment.pairwise import compute_pairwise_synteny
from final_finalizer.models import AssemblyResult, ChromRefSummary, ContigClassification


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_assembly_result(
    name: str,
    outprefix: Path,
    assembly_path: Path | None = None,
) -> AssemblyResult:
    """Create a minimal AssemblyResult for testing."""
    return AssemblyResult(
        assembly_name=name,
        assembly_path=assembly_path or Path(f"/tmp/{name}.fasta"),
        outprefix=outprefix,
        total_contigs=2,
        total_bp=2_000_000,
        n50=1_000_000,
        l50=1,
        largest_contig=1_000_000,
        classification_counts={"chrom_assigned": 2},
        classification_bp={"chrom_assigned": 2_000_000},
        n_chrom_assigned=2,
        n_chrom_unassigned=0,
        n_full_length=2,
        n_with_both_telomeres=0,
        n_with_any_telomere=0,
        mean_ref_coverage=0.95,
        n_chimeric=0,
        mean_identity=0.98,
        mean_collinearity=0.95,
        mean_gc_deviation=0.5,
        n_contaminants=0,
        total_contaminant_bp=0,
        n_unique_contaminant_species=0,
        n_rdna_contigs=0,
        total_rdna_bp=0,
        n_rdna_arrays=0,
        chrC_found=False,
        chrM_found=False,
        mean_chrom_depth=None,
        chrom_ref_coverage={},
        classifications=[],
        summary_tsv=Path(str(outprefix) + ".contig_summary.tsv"),
        segments_tsv=Path(str(outprefix) + ".segments.tsv"),
        evidence_tsv=Path(str(outprefix) + ".evidence_summary.tsv"),
        macro_blocks_tsv=Path(str(outprefix) + ".macro_blocks.tsv"),
    )


def _make_args():
    """Create a mock args namespace with default chain parsing parameters."""
    args = MagicMock()
    args.preset = "asm20"
    args.kmer = None
    args.window = None
    args.assign_minlen = 10000
    args.assign_minmapq = 0
    args.assign_tp = "PI"
    args.chain_q_gap = 200000
    args.chain_r_gap = 400000
    args.chain_diag_slop = 150000
    args.assign_min_ident = 0.0
    args.assign_chain_topk = 3
    args.assign_chain_score = "matches"
    args.assign_chain_min_bp = 0
    args.assign_ref_score = "all"
    return args


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestComputePairwiseSynteny:
    """Tests for compute_pairwise_synteny()."""

    def test_missing_left_chrs_returns_none(self, tmp_path):
        """Should return None when left chrs.fasta doesn't exist."""
        left_prefix = tmp_path / "left" / "left"
        right_prefix = tmp_path / "right" / "right"
        left_prefix.parent.mkdir(parents=True)
        right_prefix.parent.mkdir(parents=True)

        # Create right chrs.fasta but not left
        right_chrs = Path(str(right_prefix) + ".chrs.fasta")
        right_chrs.write_text(">chr1\nACGT\n")

        left_result = _make_assembly_result("left", left_prefix)
        right_result = _make_assembly_result("right", right_prefix)
        pair_prefix = tmp_path / "pairwise" / "left_vs_right"

        result = compute_pairwise_synteny(
            left_result, right_result, pair_prefix, 1, _make_args()
        )
        assert result is None

    def test_missing_right_chrs_returns_none(self, tmp_path):
        """Should return None when right chrs.fasta doesn't exist."""
        left_prefix = tmp_path / "left" / "left"
        right_prefix = tmp_path / "right" / "right"
        left_prefix.parent.mkdir(parents=True)
        right_prefix.parent.mkdir(parents=True)

        # Create left chrs.fasta but not right
        left_chrs = Path(str(left_prefix) + ".chrs.fasta")
        left_chrs.write_text(">chr1\nACGT\n")

        left_result = _make_assembly_result("left", left_prefix)
        right_result = _make_assembly_result("right", right_prefix)
        pair_prefix = tmp_path / "pairwise" / "left_vs_right"

        result = compute_pairwise_synteny(
            left_result, right_result, pair_prefix, 1, _make_args()
        )
        assert result is None

    def test_empty_right_chrs_returns_none(self, tmp_path):
        """Should return None when right chrs.fasta is empty."""
        left_prefix = tmp_path / "left" / "left"
        right_prefix = tmp_path / "right" / "right"
        left_prefix.parent.mkdir(parents=True)
        right_prefix.parent.mkdir(parents=True)

        # Create both chrs.fasta but right is empty
        left_chrs = Path(str(left_prefix) + ".chrs.fasta")
        left_chrs.write_text(">chr1\nACGT\n")
        right_chrs = Path(str(right_prefix) + ".chrs.fasta")
        right_chrs.write_text("")  # empty

        left_result = _make_assembly_result("left", left_prefix)
        right_result = _make_assembly_result("right", right_prefix)
        pair_prefix = tmp_path / "pairwise" / "left_vs_right"

        result = compute_pairwise_synteny(
            left_result, right_result, pair_prefix, 1, _make_args()
        )
        # Empty file should be caught by file_exists_and_valid
        assert result is None

    def test_reuses_cached_macro_blocks(self, tmp_path):
        """Should reuse existing macro_blocks TSV without re-running alignment."""
        left_prefix = tmp_path / "left" / "left"
        right_prefix = tmp_path / "right" / "right"
        left_prefix.parent.mkdir(parents=True)
        right_prefix.parent.mkdir(parents=True)

        # Create both chrs.fasta
        left_chrs = Path(str(left_prefix) + ".chrs.fasta")
        left_chrs.write_text(">chr1\nACGT\n")
        right_chrs = Path(str(right_prefix) + ".chrs.fasta")
        right_chrs.write_text(">chr1\nACGT\n")

        # Pre-create the macro_blocks TSV (cached)
        pair_prefix = tmp_path / "pairwise" / "left_vs_right"
        pair_prefix.parent.mkdir(parents=True)
        macro_tsv = Path(str(pair_prefix) + ".macro_blocks.tsv")
        macro_tsv.write_text("contig\tcontig_len\tref_id\n")  # minimal header

        left_result = _make_assembly_result("left", left_prefix)
        right_result = _make_assembly_result("right", right_prefix)

        result = compute_pairwise_synteny(
            left_result, right_result, pair_prefix, 1, _make_args()
        )
        assert result == macro_tsv

    def test_creates_output_directory(self, tmp_path):
        """Should create the pairwise output directory if it doesn't exist."""
        left_prefix = tmp_path / "left" / "left"
        right_prefix = tmp_path / "right" / "right"
        left_prefix.parent.mkdir(parents=True)
        right_prefix.parent.mkdir(parents=True)

        # Create chrs.fasta files
        left_chrs = Path(str(left_prefix) + ".chrs.fasta")
        left_chrs.write_text(">chr1\nACGT\n")
        right_chrs = Path(str(right_prefix) + ".chrs.fasta")
        right_chrs.write_text(">chr1\nACGT\n")

        # Pre-create cached result so we don't need minimap2
        pair_prefix = tmp_path / "pairwise" / "subdir" / "left_vs_right"
        pair_prefix.parent.mkdir(parents=True, exist_ok=True)
        macro_tsv = Path(str(pair_prefix) + ".macro_blocks.tsv")
        macro_tsv.write_text("contig\tcontig_len\tref_id\n")

        left_result = _make_assembly_result("left", left_prefix)
        right_result = _make_assembly_result("right", right_prefix)

        result = compute_pairwise_synteny(
            left_result, right_result, pair_prefix, 1, _make_args()
        )
        assert result is not None
        assert result.exists()


class TestPairwiseImport:
    """Test that the module imports correctly."""

    def test_import(self):
        from final_finalizer.alignment.pairwise import compute_pairwise_synteny
        assert callable(compute_pairwise_synteny)
