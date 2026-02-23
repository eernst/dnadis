"""Tests for multi-assembly comparison output."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from final_finalizer.models import (
    AssemblyResult,
    ChainEvidenceResult,
    ChromRefSummary,
    ContigClassification,
    DepthStats,
)
from final_finalizer.output.comparison import (
    _compute_n50_l50,
    build_assembly_result,
    write_chromosome_completeness_tsv,
    write_comparison_summary_tsv,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_classification(
    name: str,
    classification: str = "chrom_assigned",
    contig_len: int = 1_000_000,
    assigned_ref_id: str = "chr1",
    ref_coverage: float = 0.95,
    is_full_length: bool = True,
    has_5p_telomere: bool = False,
    has_3p_telomere: bool = False,
    seq_identity: float = 0.98,
    collinearity: float = 0.95,
    gc_deviation: float = 0.5,
    contaminant_sci: str = None,
    contaminant_taxid: int = None,
) -> ContigClassification:
    return ContigClassification(
        original_name=name,
        new_name=name,
        classification=classification,
        reversed=False,
        contaminant_taxid=contaminant_taxid,
        contaminant_sci=contaminant_sci,
        assigned_ref_id=assigned_ref_id if classification == "chrom_assigned" else None,
        ref_gene_proportion=0.8 if classification == "chrom_assigned" else None,
        contig_len=contig_len,
        ref_coverage=ref_coverage if classification == "chrom_assigned" else None,
        is_full_length=is_full_length if classification == "chrom_assigned" else None,
        has_5p_telomere=has_5p_telomere if classification == "chrom_assigned" else None,
        has_3p_telomere=has_3p_telomere if classification == "chrom_assigned" else None,
        seq_identity_vs_ref=seq_identity if classification == "chrom_assigned" else None,
        collinearity_score=collinearity if classification == "chrom_assigned" else None,
        gc_deviation=gc_deviation if classification == "chrom_assigned" else None,
    )


def _make_empty_ev() -> ChainEvidenceResult:
    """Create a minimal ChainEvidenceResult for testing."""
    return ChainEvidenceResult(
        qlens_from_paf={},
        qr_union_bp={},
        qr_matches={},
        qr_alnlen={},
        qr_gene_count={},
        contig_total={},
        contig_refs={},
        qr_score_topk={},
        qr_weight_all={},
        qr_nchains_kept={},
        best_ref={},
        best_score={},
        best_bp={},
        second_ref={},
        second_score={},
        second_bp={},
        chain_segments_rows=[],
        macro_block_rows=[],
        chain_summary_rows=[],
    )


def _make_assembly_result(
    name: str = "asm1",
    classifications: list = None,
    qry_lengths: dict = None,
    ref_lengths_norm: dict = None,
    tmp_path: Path = None,
) -> AssemblyResult:
    """Build an AssemblyResult using build_assembly_result with defaults."""
    if classifications is None:
        classifications = [
            _make_classification("ctg1", "chrom_assigned", 30_000_000, "chr1"),
            _make_classification("ctg2", "chrom_assigned", 25_000_000, "chr2"),
            _make_classification("ctg3", "rDNA", 50_000),
            _make_classification("ctg4", "unclassified", 10_000),
        ]
    if qry_lengths is None:
        qry_lengths = {c.original_name: c.contig_len for c in classifications}
    if ref_lengths_norm is None:
        ref_lengths_norm = {"chr1": 30_000_000, "chr2": 25_000_000}

    base = tmp_path or Path("/tmp")
    outprefix = base / name / name

    ev = _make_empty_ev()
    # Populate ev fields for chrom_assigned contigs
    for c in classifications:
        if c.classification == "chrom_assigned" and c.assigned_ref_id:
            ev.contig_total[c.original_name] = c.contig_len
            ev.contig_refs[c.original_name] = {c.assigned_ref_id}
            ev.best_bp[c.original_name] = c.contig_len
            ev.second_bp[c.original_name] = 0

    return build_assembly_result(
        assembly_name=name,
        assembly_path=base / f"{name}.fa",
        outprefix=outprefix,
        classifications=classifications,
        qry_lengths=qry_lengths,
        ref_lengths_norm=ref_lengths_norm,
        ev=ev,
        contaminants_filtered={},
        chrC_contig=None,
        chrM_contig=None,
        rdna_arrays=[],
        depth_stats={},
        chimera_primary_frac=0.8,
        chimera_secondary_frac=0.2,
        summary_tsv=base / f"{name}.contig_summary.tsv",
        segments_tsv=base / f"{name}.segments.tsv",
        evidence_tsv=base / f"{name}.evidence_summary.tsv",
        macro_blocks_tsv=base / f"{name}.macro_blocks.tsv",
    )


# ---------------------------------------------------------------------------
# N50/L50 computation
# ---------------------------------------------------------------------------
class TestN50L50:
    def test_simple(self):
        """Known distribution: [10, 8, 6, 4, 2] total=30, half=15."""
        n50, l50 = _compute_n50_l50([10, 8, 6, 4, 2])
        assert n50 == 8
        assert l50 == 2

    def test_single_contig(self):
        n50, l50 = _compute_n50_l50([100])
        assert n50 == 100
        assert l50 == 1

    def test_empty(self):
        n50, l50 = _compute_n50_l50([])
        assert n50 == 0
        assert l50 == 0

    def test_equal_lengths(self):
        """Five contigs of equal length."""
        n50, l50 = _compute_n50_l50([100, 100, 100, 100, 100])
        assert n50 == 100
        assert l50 == 3  # 3 * 100 = 300 >= 250

    def test_two_contigs(self):
        """Two contigs: [100, 1] total=101, half=50.5."""
        n50, l50 = _compute_n50_l50([100, 1])
        assert n50 == 100
        assert l50 == 1

    def test_unsorted_input(self):
        """Input does not need to be pre-sorted."""
        n50, l50 = _compute_n50_l50([2, 10, 4, 8, 6])
        assert n50 == 8
        assert l50 == 2


# ---------------------------------------------------------------------------
# build_assembly_result
# ---------------------------------------------------------------------------
class TestBuildAssemblyResult:
    def test_basic_fields(self, tmp_path):
        r = _make_assembly_result("test_asm", tmp_path=tmp_path)
        assert r.assembly_name == "test_asm"
        assert r.total_contigs == 4
        assert r.total_bp == 30_000_000 + 25_000_000 + 50_000 + 10_000

    def test_classification_counts(self, tmp_path):
        r = _make_assembly_result(tmp_path=tmp_path)
        assert r.classification_counts["chrom_assigned"] == 2
        assert r.classification_counts["rDNA"] == 1
        assert r.classification_counts["unclassified"] == 1

    def test_classification_bp(self, tmp_path):
        r = _make_assembly_result(tmp_path=tmp_path)
        assert r.classification_bp["chrom_assigned"] == 30_000_000 + 25_000_000

    def test_n50_computation(self, tmp_path):
        r = _make_assembly_result(tmp_path=tmp_path)
        assert r.n50 > 0
        assert r.l50 > 0
        assert r.largest_contig == 30_000_000

    def test_chromosome_completeness(self, tmp_path):
        r = _make_assembly_result(tmp_path=tmp_path)
        assert r.n_chrom_assigned == 2
        assert r.n_full_length == 2  # both are full-length by default

    def test_quality_metrics(self, tmp_path):
        r = _make_assembly_result(tmp_path=tmp_path)
        assert r.mean_identity == pytest.approx(0.98)
        assert r.mean_collinearity == pytest.approx(0.95)
        assert r.mean_gc_deviation == pytest.approx(0.5)

    def test_no_chrom_assigned(self, tmp_path):
        """Assembly with no chrom_assigned contigs."""
        clfs = [
            _make_classification("ctg1", "unclassified", 10_000),
            _make_classification("ctg2", "rDNA", 50_000),
        ]
        r = _make_assembly_result(
            "no_chrom", classifications=clfs,
            ref_lengths_norm={"chr1": 30_000_000},
            tmp_path=tmp_path,
        )
        assert r.n_chrom_assigned == 0
        assert r.mean_identity is None
        assert r.mean_ref_coverage is None

    def test_chrom_ref_coverage(self, tmp_path):
        """Per-ref-chromosome summaries include all ref chromosomes."""
        r = _make_assembly_result(tmp_path=tmp_path)
        assert "chr1" in r.chrom_ref_coverage
        assert "chr2" in r.chrom_ref_coverage
        crs1 = r.chrom_ref_coverage["chr1"]
        assert crs1.n_contigs == 1
        assert crs1.is_full_length is True

    def test_missing_ref_chromosome(self, tmp_path):
        """Ref chromosomes with no assigned contigs appear with n_contigs=0."""
        clfs = [
            _make_classification("ctg1", "chrom_assigned", 30_000_000, "chr1"),
        ]
        ref_lengths = {"chr1": 30_000_000, "chr2": 25_000_000, "chr3": 20_000_000}
        r = _make_assembly_result(
            classifications=clfs, ref_lengths_norm=ref_lengths, tmp_path=tmp_path,
        )
        assert r.chrom_ref_coverage["chr2"].n_contigs == 0
        assert r.chrom_ref_coverage["chr3"].n_contigs == 0
        assert r.chrom_ref_coverage["chr1"].n_contigs == 1

    def test_organelle_detection(self, tmp_path):
        """chrC_found and chrM_found reflect organelle detection."""
        ev = _make_empty_ev()
        base = tmp_path
        r = build_assembly_result(
            assembly_name="org",
            assembly_path=base / "org.fa",
            outprefix=base / "org" / "org",
            classifications=[_make_classification("ctg1", "unclassified", 10_000)],
            qry_lengths={"ctg1": 10_000},
            ref_lengths_norm={},
            ev=ev,
            contaminants_filtered={},
            chrC_contig="chrC_ctg",
            chrM_contig=None,
            rdna_arrays=[],
            depth_stats={},
            chimera_primary_frac=0.8,
            chimera_secondary_frac=0.2,
            summary_tsv=base / "s.tsv",
            segments_tsv=base / "seg.tsv",
            evidence_tsv=base / "ev.tsv",
            macro_blocks_tsv=base / "mb.tsv",
        )
        assert r.chrC_found is True
        assert r.chrM_found is False

    def test_contaminant_metrics(self, tmp_path):
        clfs = [
            _make_classification(
                "ctg1", "contaminant", 50_000,
                contaminant_sci="E. coli", contaminant_taxid=562,
            ),
            _make_classification(
                "ctg2", "contaminant", 30_000,
                contaminant_sci="E. coli", contaminant_taxid=562,
            ),
            _make_classification(
                "ctg3", "contaminant", 20_000,
                contaminant_sci="S. aureus", contaminant_taxid=1280,
            ),
        ]
        r = _make_assembly_result(
            classifications=clfs, ref_lengths_norm={}, tmp_path=tmp_path,
        )
        assert r.n_contaminants == 3
        assert r.total_contaminant_bp == 100_000
        assert r.n_unique_contaminant_species == 2

    def test_telomere_counts(self, tmp_path):
        clfs = [
            _make_classification(
                "ctg1", "chrom_assigned", 30_000_000, "chr1",
                has_5p_telomere=True, has_3p_telomere=True,
            ),
            _make_classification(
                "ctg2", "chrom_assigned", 25_000_000, "chr2",
                has_5p_telomere=True, has_3p_telomere=False,
            ),
            _make_classification(
                "ctg3", "chrom_assigned", 20_000_000, "chr3",
                has_5p_telomere=False, has_3p_telomere=False,
            ),
        ]
        r = _make_assembly_result(
            classifications=clfs,
            ref_lengths_norm={"chr1": 30_000_000, "chr2": 25_000_000, "chr3": 20_000_000},
            tmp_path=tmp_path,
        )
        assert r.n_with_both_telomeres == 1
        assert r.n_with_any_telomere == 2


# ---------------------------------------------------------------------------
# write_comparison_summary_tsv
# ---------------------------------------------------------------------------
class TestWriteComparisonSummaryTsv:
    def test_writes_correct_columns(self, tmp_path):
        r = _make_assembly_result("asm1", tmp_path=tmp_path)
        out = tmp_path / "comparison.tsv"
        write_comparison_summary_tsv(out, [r])

        lines = out.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header[0] == "assembly"
        assert "n50" in header
        assert "mean_identity" in header
        assert len(lines) == 2  # header + 1 data row

    def test_multiple_assemblies(self, tmp_path):
        r1 = _make_assembly_result("asm1", tmp_path=tmp_path)
        r2 = _make_assembly_result("asm2", tmp_path=tmp_path)
        out = tmp_path / "comparison.tsv"
        write_comparison_summary_tsv(out, [r1, r2])

        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 data rows
        assert lines[1].startswith("asm1\t")
        assert lines[2].startswith("asm2\t")

    def test_column_values(self, tmp_path):
        r = _make_assembly_result("asm1", tmp_path=tmp_path)
        out = tmp_path / "comparison.tsv"
        write_comparison_summary_tsv(out, [r])

        lines = out.read_text().strip().split("\n")
        header = lines[0].split("\t")
        data = lines[1].split("\t")
        col_map = dict(zip(header, data))
        assert col_map["assembly"] == "asm1"
        assert int(col_map["total_contigs"]) == 4
        assert int(col_map["n_chrom_assigned"]) == 2


# ---------------------------------------------------------------------------
# write_chromosome_completeness_tsv
# ---------------------------------------------------------------------------
class TestWriteChromosomeCompletenessTsv:
    def test_complete_matrix(self, tmp_path):
        """Every (ref, assembly) pair has a row, including missing assignments."""
        ref_lengths = {"chr1": 30_000_000, "chr2": 25_000_000}

        # asm1 has chr1 assigned, asm2 has chr2 assigned
        clfs1 = [_make_classification("ctg1", "chrom_assigned", 30_000_000, "chr1")]
        clfs2 = [_make_classification("ctg2", "chrom_assigned", 25_000_000, "chr2")]

        r1 = _make_assembly_result("asm1", clfs1, ref_lengths_norm=ref_lengths, tmp_path=tmp_path)
        r2 = _make_assembly_result("asm2", clfs2, ref_lengths_norm=ref_lengths, tmp_path=tmp_path)

        out = tmp_path / "completeness.tsv"
        write_chromosome_completeness_tsv(out, [r1, r2], ref_lengths)

        lines = out.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header[0] == "ref_id"
        # 2 refs × 2 assemblies = 4 data rows
        assert len(lines) == 5  # header + 4 data rows

    def test_missing_assignment_has_zero(self, tmp_path):
        """Assembly without a ref chromosome gets n_contigs=0."""
        ref_lengths = {"chr1": 30_000_000, "chr2": 25_000_000}
        clfs = [_make_classification("ctg1", "chrom_assigned", 30_000_000, "chr1")]
        r = _make_assembly_result("asm1", clfs, ref_lengths_norm=ref_lengths, tmp_path=tmp_path)

        out = tmp_path / "completeness.tsv"
        write_chromosome_completeness_tsv(out, [r], ref_lengths)

        lines = out.read_text().strip().split("\n")
        header = lines[0].split("\t")
        # Find the chr2 row
        data_rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
        chr2_rows = [d for d in data_rows if d["ref_id"] == "chr2"]
        assert len(chr2_rows) == 1
        assert chr2_rows[0]["n_contigs"] == "0"
        assert chr2_rows[0]["is_full_length"] == "no"

    def test_single_assembly(self, tmp_path):
        ref_lengths = {"chr1": 30_000_000}
        clfs = [_make_classification("ctg1", "chrom_assigned", 30_000_000, "chr1")]
        r = _make_assembly_result("asm1", clfs, ref_lengths_norm=ref_lengths, tmp_path=tmp_path)

        out = tmp_path / "completeness.tsv"
        write_chromosome_completeness_tsv(out, [r], ref_lengths)

        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2  # header + 1 row
        header = lines[0].split("\t")
        data = dict(zip(header, lines[1].split("\t")))
        assert data["ref_id"] == "chr1"
        assert data["assembly"] == "asm1"
        assert int(data["n_contigs"]) == 1
