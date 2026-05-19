"""Unit tests for the phylogeny pipeline (BUSCO grouping, outgroup logic).

These do not exercise MAFFT / trimAl / IQ-TREE — those are run as part of
the end-to-end integration test elsewhere.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnadis.models import CompleasmResult, ContigClassification
from dnadis.phylogeny.busco_extraction import (
    LeafBuscos,
    LeafId,
    build_assembly_leaves,
    build_reference_leaf,
    filter_leaves_by_completeness,
    parse_full_table,
    shared_single_copy_genes,
)
from dnadis.phylogeny.tree import (
    _pair_identity,
    auto_outgroup,
    resolve_outgroup,
)


# --------------------------------------------------------------------------
# parse_full_table
# --------------------------------------------------------------------------
def _write_full_table(path: Path, rows: list) -> Path:
    """Write a compleasm-style full_table_busco_format.tsv."""
    lines = [
        "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength"
    ]
    for r in rows:
        lines.append("\t".join(str(c) for c in r))
    path.write_text("\n".join(lines) + "\n")
    return path


def test_parse_full_table_includes_complete_and_missing(tmp_path):
    p = _write_full_table(
        tmp_path / "full_table.tsv",
        [
            ("g1", "Complete", "chr1", 100, 200, "+", 99.0, 100),
            ("g2", "Missing"),
            ("g3", "Duplicated", "chr1A", 1, 50, "+", 80.0, 49),
            ("g3", "Duplicated", "chr1B", 1, 50, "+", 80.0, 49),
        ],
    )
    hits = parse_full_table(p)
    assert [h.busco_id for h in hits] == ["g1", "g2", "g3", "g3"]
    assert hits[0].contig == "chr1" and hits[0].start == 100
    assert hits[1].contig is None and hits[1].status == "Missing"
    assert hits[2].status == "Duplicated" and hits[2].contig == "chr1A"


def test_parse_full_table_missing_file(tmp_path):
    assert parse_full_table(tmp_path / "nope.tsv") == []


# --------------------------------------------------------------------------
# build_assembly_leaves: subgenome grouping
# --------------------------------------------------------------------------
def _mk_classification(contig: str, ref_id: str, subgenome=None) -> ContigClassification:
    return ContigClassification(
        original_name=contig,
        new_name=contig,
        classification="chrom_assigned",
        reversed=False,
        contaminant_taxid=None,
        contaminant_sci=None,
        assigned_ref_id=ref_id,
        ref_gene_proportion=None,
        contig_len=1_000_000,
        query_subgenome=subgenome,
    )


def _mk_compleasm(full_table: Path, n_total: int = 100) -> CompleasmResult:
    return CompleasmResult(
        lineage="eukaryota_odb12",
        n_total=n_total,
        n_single=n_total,
        n_duplicated=0,
        n_fragmented=0,
        n_interspersed=0,
        n_missing=0,
        pct_single=100.0,
        pct_duplicated=0.0,
        pct_fragmented=0.0,
        pct_interspersed=0.0,
        pct_missing=0.0,
        summary_path=full_table.parent / "summary.txt",
        full_table_path=full_table,
    )


def test_build_assembly_leaves_splits_by_reference_subgenome(tmp_path):
    """Contigs assigned to chr*A and chr*P fall into distinct ref-subgenome leaves.

    Verifies that BUSCOs labelled ``Duplicated`` globally (one hit per
    reference subgenome) still count as single-copy *within each leaf*.
    """
    ft = _write_full_table(
        tmp_path / "full_table.tsv",
        [
            ("gA", "Duplicated", "ctg1A", 1, 100, "+", 90.0, 99),
            ("gA", "Duplicated", "ctg1P", 1, 100, "+", 90.0, 99),
            ("gB", "Duplicated", "ctg1A", 200, 300, "+", 90.0, 99),
            ("gB", "Duplicated", "ctg1P", 200, 300, "+", 90.0, 99),
            ("gC", "Complete",   "ctg1P", 400, 500, "+", 90.0, 99),  # P-only
            ("gD", "Duplicated", "ctg1A", 600, 700, "+", 90.0, 99),  # extra copy on A
            ("gD", "Duplicated", "ctg1A", 800, 900, "+", 90.0, 99),
            ("gD", "Duplicated", "ctg1P", 600, 700, "+", 90.0, 99),
            ("noise", "Complete", "scratch", 1, 50, "+", 50.0, 49),
        ],
    )
    classifications = [
        _mk_classification("ctg1A", "chr1A"),
        _mk_classification("ctg1P", "chr1P"),
        _mk_classification("scratch", "debris"),
    ]
    classifications[2].classification = "debris"

    leaves = build_assembly_leaves("asm1", classifications, _mk_compleasm(ft))
    by_label = {lb.leaf.label: lb for lb in leaves}
    assert set(by_label) == {"asm1_A", "asm1_P"}

    # Globally Duplicated but locally single-copy in each subgenome → both
    # leaves keep gA and gB; gC is P-only; gD is duplicated in A and single in P.
    a = by_label["asm1_A"]
    assert a.single_copy_genes() == {"gA", "gB"}

    p = by_label["asm1_P"]
    assert p.single_copy_genes() == {"gA", "gB", "gC", "gD"}


def test_build_assembly_leaves_combines_ref_and_query_subgenome(tmp_path):
    """Polyploid query on polyploid reference yields ref × query subgenome leaves."""
    ft = _write_full_table(
        tmp_path / "full_table.tsv",
        [
            ("g1", "Complete", "ctg1A_primary", 1, 100, "+", 90.0, 99),
            ("g1", "Duplicated", "ctg1A_alt", 1, 100, "+", 90.0, 99),
            ("g1", "Duplicated", "ctg1P_primary", 1, 100, "+", 90.0, 99),
            ("g1", "Duplicated", "ctg1P_alt", 1, 100, "+", 90.0, 99),
        ],
    )
    classifications = [
        _mk_classification("ctg1A_primary", "chr1A", subgenome=None),
        _mk_classification("ctg1A_alt",     "chr1A", subgenome="B"),
        _mk_classification("ctg1P_primary", "chr1P", subgenome=None),
        _mk_classification("ctg1P_alt",     "chr1P", subgenome="B"),
    ]
    leaves = build_assembly_leaves("asm1", classifications, _mk_compleasm(ft))
    assert {lb.leaf.label for lb in leaves} == {
        "asm1_A", "asm1_A_B", "asm1_P", "asm1_P_B",
    }


def test_build_assembly_leaves_skips_unclassified_contigs(tmp_path):
    ft = _write_full_table(
        tmp_path / "ft.tsv",
        [("g1", "Complete", "junk_contig", 1, 100, "+", 90.0, 99)],
    )
    # No matching ContigClassification → hit is dropped
    leaves = build_assembly_leaves("asm", [], _mk_compleasm(ft))
    assert leaves == []


def test_build_assembly_leaves_returns_empty_when_no_full_table(tmp_path):
    cr = CompleasmResult(
        lineage="x",
        n_total=10, n_single=5, n_duplicated=0, n_fragmented=0,
        n_interspersed=0, n_missing=5,
        pct_single=50.0, pct_duplicated=0.0, pct_fragmented=0.0,
        pct_interspersed=0.0, pct_missing=50.0,
        summary_path=tmp_path / "s.txt",
        full_table_path=None,
    )
    assert build_assembly_leaves("asm", [], cr) == []


# --------------------------------------------------------------------------
# Filtering by completeness
# --------------------------------------------------------------------------
def _mk_leaf(label: str, n_total: int, complete_genes: list) -> LeafBuscos:
    parts = label.split("_")
    leaf = LeafId(
        assembly=parts[0],
        ref_subgenome=parts[1] if len(parts) >= 2 and parts[1] else None,
        query_subgenome=parts[2] if len(parts) >= 3 and parts[2] else None,
    )
    from dnadis.phylogeny.busco_extraction import BuscoHit
    hits = {
        g: [BuscoHit(g, "Complete", "ctg", 1, 100, "+", 99.0, 99)]
        for g in complete_genes
    }
    return LeafBuscos(leaf=leaf, hits_by_busco=hits, translated_protein_path=None, n_lineage_total=n_total)


def test_filter_by_completeness_drops_below_threshold():
    leaves = [
        _mk_leaf("a", 10, ["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8"]),  # 80%
        _mk_leaf("b", 10, ["g1", "g2", "g3", "g4", "g5"]),                    # 50%
        _mk_leaf("c", 10, ["g1", "g2", "g3", "g4", "g5", "g6", "g7"]),        # 70%
    ]
    kept, dropped = filter_leaves_by_completeness(leaves, min_pct=70.0)
    assert {lb.leaf.label for lb in kept} == {"a", "c"}
    assert {lb.leaf.label for lb, _ in dropped} == {"b"}


# --------------------------------------------------------------------------
# Shared single-copy genes (intersection across all leaves)
# --------------------------------------------------------------------------
def test_shared_single_copy_genes_intersects():
    a = _mk_leaf("a", 10, ["g1", "g2", "g3"])
    b = _mk_leaf("b", 10, ["g1", "g3", "g4"])
    c = _mk_leaf("c", 10, ["g1", "g3", "g5"])
    assert shared_single_copy_genes([a, b, c]) == ["g1", "g3"]


def test_shared_single_copy_genes_empty_input_returns_empty():
    assert shared_single_copy_genes([]) == []


# --------------------------------------------------------------------------
# Outgroup resolution
# --------------------------------------------------------------------------
def test_resolve_outgroup_none_returns_none(tmp_path):
    assert resolve_outgroup("none", ["a", "b", "c"], "ref", tmp_path / "sm.faa") is None


def test_resolve_outgroup_reference_returns_label(tmp_path):
    assert resolve_outgroup("reference", ["a", "b", "TAIR10"], "TAIR10", tmp_path / "sm.faa") == "TAIR10"


def test_resolve_outgroup_reference_warns_when_not_in_tree(tmp_path):
    # Reference label is None when --phylo-skip-reference was set; resolve returns None and warns
    assert resolve_outgroup("reference", ["a", "b"], None, tmp_path / "sm.faa") is None


def test_resolve_outgroup_matches_assembly_with_single_subgenome():
    # An assembly with one subgenome "B" should match by bare assembly name
    assert resolve_outgroup("asm1", ["asm1_B", "other"], None, Path("/tmp/x")) == "asm1_B"


def test_resolve_outgroup_unknown_label_returns_none(tmp_path):
    assert resolve_outgroup("nonexistent", ["a", "b"], None, tmp_path / "sm.faa") is None


def test_pair_identity_basic():
    assert _pair_identity("AAAA", "AAAA") == 1.0
    assert _pair_identity("AAAA", "ATAT") == 0.5
    # All-gap columns ignored
    assert _pair_identity("A-A-", "A-A-") == 1.0
    # Mismatch
    assert _pair_identity("AAAA", "TTTT") == 0.0


def test_auto_outgroup_picks_most_divergent(tmp_path):
    fa = tmp_path / "sm.faa"
    fa.write_text(
        ">a\nAAAAAAAA\n"
        ">b\nAAAAAAAA\n"
        ">c\nTTTTTTTT\n"
    )
    # 'c' shares zero identity with a and b → most divergent
    assert auto_outgroup(fa) == "c"


def test_auto_outgroup_two_taxa_returns_none(tmp_path):
    fa = tmp_path / "sm.faa"
    fa.write_text(">a\nAAAA\n>b\nAAAA\n")
    assert auto_outgroup(fa) is None
