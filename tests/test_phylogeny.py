"""Unit tests for the phylogeny pipeline (BUSCO grouping, outgroup logic).

These do not exercise MAFFT / trimAl / IQ-TREE — those are run as part of
the end-to-end integration test elsewhere.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnadis.detection.compleasm import derive_chrs_non_chrs_compleasm
from dnadis.models import CompleasmResult, ContigClassification
from dnadis.phylogeny.busco_extraction import (
    LeafBuscos,
    LeafId,
    build_assembly_leaves,
    build_outgroup_leaves,
    build_outgroup_only_leaves,
    build_reference_leaves,
    filter_leaves_by_completeness,
    infer_assembly_ploidy,
    parse_full_table,
    ref_assignment_quality_by_subgenome,
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
        cobiont_taxid=None,
        cobiont_sci=None,
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


# --------------------------------------------------------------------------
# Derived chrs / non_chrs split from a full-assembly compleasm run
# --------------------------------------------------------------------------
def test_derive_split_counts_correctly(tmp_path):
    """Full-table joined with classifications produces consistent chrs/non_chrs counts."""
    ft = _write_full_table(
        tmp_path / "full_table.tsv",
        [
            # gA: single hit on a chrom_assigned contig → chrs:Complete
            ("gA", "Complete",   "ctg1A", 1, 100, "+", 90.0, 99),
            # gB: hits on both chrs and non_chrs contigs
            ("gB", "Duplicated", "ctg1A", 1, 100, "+", 90.0, 99),
            ("gB", "Duplicated", "ctgX",  1, 100, "+", 90.0, 99),
            # gC: two hits on debris → non_chrs:Duplicated
            ("gC", "Duplicated", "ctgX",  200, 300, "+", 90.0, 99),
            ("gC", "Duplicated", "ctgY",  200, 300, "+", 90.0, 99),
            # gD: fragmented on debris → non_chrs:Fragmented
            ("gD", "Fragmented", "ctgX",  500, 600, "+", 30.0, 99),
        ],
    )
    classifications = [
        _mk_classification("ctg1A", "chr1A"),
        _mk_classification("ctgX",  "debris"),
        _mk_classification("ctgY",  "debris"),
    ]
    classifications[1].classification = "debris"
    classifications[2].classification = "debris"

    full = _mk_compleasm(ft, n_total=10)
    chrs_dir = tmp_path / "chrs"
    non_dir = tmp_path / "non"
    chrs_res, non_res = derive_chrs_non_chrs_compleasm(
        full, classifications, chrs_dir, non_dir,
    )

    # chrs side: gA (1 hit) → S; gB (1 hit) → S; others not present
    assert chrs_res is not None
    assert chrs_res.n_single == 2
    assert chrs_res.n_duplicated == 0
    assert chrs_res.n_fragmented == 0
    assert chrs_res.n_missing == 10 - 2

    # non_chrs side: gB (1 hit) → S; gC (2 hits) → D; gD (1 frag-only) → F
    assert non_res is not None
    assert non_res.n_single == 1
    assert non_res.n_duplicated == 1
    assert non_res.n_fragmented == 1
    assert non_res.n_missing == 10 - 3

    # Sanity: derived summary parses back through parse_compleasm_summary
    assert (chrs_dir / "summary.txt").exists()
    assert (non_dir / "summary.txt").exists()
    # And the derived full_tables are persisted under the lineage subdir
    assert (chrs_dir / full.lineage / "full_table_busco_format.tsv").exists()
    assert (non_dir / full.lineage / "full_table_busco_format.tsv").exists()


# --------------------------------------------------------------------------
# Outgroup helpers: ploidy inference and ref-assignment quality
# --------------------------------------------------------------------------
def test_infer_assembly_ploidy_haploid(tmp_path):
    rows = parse_full_table(_write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Complete", "c", 1, 100, "+", 90.0, 99),
        ("g2", "Complete", "c", 200, 300, "+", 90.0, 99),
        ("g3", "Complete", "c", 400, 500, "+", 90.0, 99),
    ]))
    assert infer_assembly_ploidy(rows) == 1


def test_infer_assembly_ploidy_polyploid(tmp_path):
    # Each gene appears twice (typical allotetraploid) → median = 2
    rows = parse_full_table(_write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Duplicated", "cA", 1, 100, "+", 90.0, 99),
        ("g1", "Duplicated", "cP", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "cA", 200, 300, "+", 90.0, 99),
        ("g2", "Duplicated", "cP", 200, 300, "+", 90.0, 99),
    ]))
    assert infer_assembly_ploidy(rows) == 2


def test_ref_assignment_quality_distant_outgroup():
    """7 assignments scattered across a 30-chrom reference → low quality everywhere."""
    classifications = [
        _mk_classification("c1", "chr1A"),
        _mk_classification("c2", "chr2A"),
        _mk_classification("c3", "chr3A"),
        _mk_classification("c4", "chr1P"),
        _mk_classification("c5", "chr2P"),
        _mk_classification("c6", "chr1T"),
        _mk_classification("c7", "chr2T"),
    ]
    # Pretend the reference has 10 chroms per subgenome
    ref_lengths = {}
    for sg in ("A", "P", "T"):
        for i in range(1, 11):
            ref_lengths[f"chr{i}{sg}"] = 1_000_000
    quality = ref_assignment_quality_by_subgenome(classifications, ref_lengths)
    # A subgenome: 3 chroms assigned out of 10 → 0.3
    assert abs(quality["A"] - 0.3) < 1e-9
    assert abs(quality["P"] - 0.2) < 1e-9
    assert abs(quality["T"] - 0.2) < 1e-9


def test_ref_assignment_quality_well_mapped_one_subgenome():
    """One subgenome densely mapped, others not → keep the dense one."""
    classifications = [_mk_classification(f"c{i}", f"chr{i}A") for i in range(1, 11)]
    ref_lengths = {}
    for sg in ("A", "P"):
        for i in range(1, 11):
            ref_lengths[f"chr{i}{sg}"] = 1_000_000
    quality = ref_assignment_quality_by_subgenome(classifications, ref_lengths)
    assert abs(quality["A"] - 1.0) < 1e-9
    assert quality.get("P", 0.0) == 0.0  # no contigs assigned to P


# --------------------------------------------------------------------------
# build_outgroup_leaves: status decisions
# --------------------------------------------------------------------------
def test_build_outgroup_leaves_haploid_pools_all(tmp_path):
    """Haploid outgroup retains BUSCOs regardless of chromosome assignment."""
    ft = _write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Complete", "ctg_chrom_assigned", 1, 100, "+", 90.0, 99),
        ("g2", "Complete", "ctg_unassigned",     1, 100, "+", 90.0, 99),
        ("g3", "Complete", "ctg_debris",         1, 100, "+", 90.0, 99),
    ])
    classifications = [
        _mk_classification("ctg_chrom_assigned", "chr1A"),
        _mk_classification("ctg_unassigned",     "chr1A"),  # ref_id set ...
        _mk_classification("ctg_debris",         "debris"),
    ]
    classifications[1].classification = "chrom_unassigned"
    classifications[1].assigned_ref_id = None
    classifications[2].classification = "debris"

    full = _mk_compleasm(ft, n_total=10)
    leaves, status = build_outgroup_leaves(
        "outg", classifications, full, ref_lengths_norm={"chr1A": 100_000}, min_ref_assignment=0.5,
    )
    assert status == "haploid_pooled"
    assert len(leaves) == 1
    assert leaves[0].leaf.label == "outg"  # no ref_sg / query_sg suffix
    assert leaves[0].single_copy_genes() == {"g1", "g2", "g3"}


def test_build_outgroup_leaves_polyploid_unusable(tmp_path):
    """Polyploid outgroup with sparse ref-chromosome assignments → unusable."""
    ft = _write_full_table(tmp_path / "ft.tsv", [
        # 4 BUSCOs, each duplicated → ploidy median == 2
        ("g1", "Duplicated", "cA", 1, 100, "+", 90.0, 99),
        ("g1", "Duplicated", "cP", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "cA", 200, 300, "+", 90.0, 99),
        ("g2", "Duplicated", "cP", 200, 300, "+", 90.0, 99),
        ("g3", "Duplicated", "cA", 400, 500, "+", 90.0, 99),
        ("g3", "Duplicated", "cP", 400, 500, "+", 90.0, 99),
        ("g4", "Duplicated", "cA", 600, 700, "+", 90.0, 99),
        ("g4", "Duplicated", "cP", 600, 700, "+", 90.0, 99),
    ])
    # Only 1/10 chroms assigned per subgenome → 0.1 quality, below 0.5 threshold
    classifications = [
        _mk_classification("cA", "chr1A"),
        _mk_classification("cP", "chr1P"),
    ]
    ref_lengths = {}
    for sg in ("A", "P"):
        for i in range(1, 11):
            ref_lengths[f"chr{i}{sg}"] = 1_000_000
    leaves, status = build_outgroup_leaves(
        "outg", classifications, _mk_compleasm(ft, n_total=10),
        ref_lengths_norm=ref_lengths, min_ref_assignment=0.5,
    )
    assert status == "polyploid_unusable"
    assert leaves == []


def test_build_reference_leaves_monoploid(tmp_path):
    """Monoploid reference contigs (no subgenome suffix) → one leaf."""
    ft = _write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Complete", "chr1", 1, 100, "+", 90.0, 99),
        ("g2", "Complete", "chr2", 1, 100, "+", 90.0, 99),
    ])
    leaves = build_reference_leaves("ref", _mk_compleasm(ft, n_total=10))
    assert len(leaves) == 1
    assert leaves[0].leaf.label == "ref"
    assert leaves[0].single_copy_genes() == {"g1", "g2"}


def test_build_reference_leaves_polyploid_splits_by_subgenome(tmp_path):
    """Composed polyploid reference → one leaf per chr subgenome suffix.

    Without this split a polyploid reference produces one pooled leaf
    where every BUSCO hit count equals the reference's ploidy, so
    ``single_copy_genes()`` drops every BUSCO and the reference falls
    below the completeness threshold.
    """
    ft = _write_full_table(tmp_path / "ft.tsv", [
        # BUSCO present once per subgenome (allopolyploid composed reference)
        ("g1", "Duplicated", "chr1A", 1, 100, "+", 90.0, 99),
        ("g1", "Duplicated", "chr1P", 1, 100, "+", 90.0, 99),
        ("g1", "Duplicated", "chr1T", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "chr2A", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "chr2P", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "chr2T", 1, 100, "+", 90.0, 99),
    ])
    leaves = build_reference_leaves("ref", _mk_compleasm(ft, n_total=10))
    by_label = {lb.leaf.label: lb for lb in leaves}
    assert set(by_label) == {"ref_A", "ref_P", "ref_T"}
    for label in ("ref_A", "ref_P", "ref_T"):
        # Each subgenome leaf sees each BUSCO exactly once → single-copy
        assert by_label[label].single_copy_genes() == {"g1", "g2"}


def test_resolve_outgroup_polyploid_reference_returns_clade(tmp_path):
    """--phylo-outgroup=reference on a polyploid ref → comma-joined clade."""
    fa = tmp_path / "sm.faa"
    fa.write_text(">ref_A\nA\n>ref_P\nA\n>ref_T\nA\n>q1\nA\n>q2\nA\n")
    result = resolve_outgroup(
        "reference",
        ["ref_A", "ref_P", "ref_T", "q1", "q2"],
        ["ref_A", "ref_P", "ref_T"],
        fa,
    )
    assert result is not None
    parts = set(result.split(","))
    assert parts == {"ref_A", "ref_P", "ref_T"}


def test_resolve_outgroup_accepts_legacy_string_reference(tmp_path):
    """Legacy single-string reference_labels arg still works."""
    fa = tmp_path / "sm.faa"
    fa.write_text(">ref\nA\n>q1\nA\n")
    assert resolve_outgroup("reference", ["ref", "q1"], "ref", fa) == "ref"


def test_build_outgroup_leaves_polyploid_filtered(tmp_path):
    """Polyploid outgroup with one well-mapped subgenome → keep just that one."""
    ft = _write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Duplicated", "cA1", 1, 100, "+", 90.0, 99),
        ("g1", "Duplicated", "cP",  1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "cA2", 200, 300, "+", 90.0, 99),
        ("g2", "Duplicated", "cP",  200, 300, "+", 90.0, 99),
    ])
    # 2/2 A-chroms assigned (quality 1.0); 1/2 P-chroms assigned (quality 0.5)
    classifications = [
        _mk_classification("cA1", "chr1A"),
        _mk_classification("cA2", "chr2A"),
        _mk_classification("cP",  "chr1P"),
    ]
    ref_lengths = {"chr1A": 1_000_000, "chr2A": 1_000_000,
                   "chr1P": 1_000_000, "chr2P": 1_000_000}
    leaves, status = build_outgroup_leaves(
        "outg", classifications, _mk_compleasm(ft, n_total=10),
        ref_lengths_norm=ref_lengths, min_ref_assignment=0.75,
    )
    assert status == "polyploid_filtered"
    labels = {lb.leaf.label for lb in leaves}
    assert labels == {"outg_A"}  # P dropped (0.5 < 0.75 threshold)


# --------------------------------------------------------------------------
# build_outgroup_only_leaves: phylogeny-only outgroups (name-based subgenomes)
# --------------------------------------------------------------------------
def test_build_outgroup_only_leaves_pools_when_no_subgenome_suffix(tmp_path):
    """A distant outgroup with no subgenome suffixes pools into one leaf."""
    ft = _write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Complete", "ctg000001", 1, 100, "+", 90.0, 99),
        ("g2", "Complete", "ctg000002", 1, 100, "+", 90.0, 99),
        ("g3", "Duplicated", "ctg000003", 1, 100, "+", 90.0, 99),
    ])
    leaves = build_outgroup_only_leaves("la6002", _mk_compleasm(ft, n_total=10))
    assert len(leaves) == 1
    assert leaves[0].leaf.label == "la6002"  # bare name, no suffix
    assert leaves[0].leaf.ref_subgenome is None
    assert leaves[0].leaf.query_subgenome is None
    assert leaves[0].single_copy_genes() == {"g1", "g2", "g3"}


def test_build_outgroup_only_leaves_splits_by_contig_name_suffix(tmp_path):
    """Subgenomes are split from the outgroup's own chr-name suffixes."""
    ft = _write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Duplicated", "chr1A", 1, 100, "+", 90.0, 99),
        ("g1", "Duplicated", "chr1B", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "chr2A", 1, 100, "+", 90.0, 99),
        ("g2", "Duplicated", "chr2B", 1, 100, "+", 90.0, 99),
    ])
    leaves = build_outgroup_only_leaves("outg", _mk_compleasm(ft, n_total=10))
    labels = {lb.leaf.label for lb in leaves}
    assert labels == {"outg_A", "outg_B"}
    by_label = {lb.leaf.label: lb for lb in leaves}
    # Each subgenome leaf is single-copy for both genes
    assert by_label["outg_A"].single_copy_genes() == {"g1", "g2"}
    assert by_label["outg_B"].single_copy_genes() == {"g1", "g2"}
    # Subgenome lands in the query_subgenome slot, ref_subgenome stays None
    assert all(lb.leaf.ref_subgenome is None for lb in leaves)


def test_build_outgroup_only_leaves_no_full_table_returns_empty(tmp_path):
    cr = CompleasmResult(
        lineage="eukaryota_odb12", n_total=10, n_single=10, n_duplicated=0,
        n_fragmented=0, n_interspersed=0, n_missing=0,
        pct_single=100.0, pct_duplicated=0.0, pct_fragmented=0.0,
        pct_interspersed=0.0, pct_missing=0.0,
        summary_path=tmp_path / "summary.txt", full_table_path=None,
    )
    assert build_outgroup_only_leaves("outg", cr) == []


def test_build_outgroup_only_leaves_skips_missing_and_fragmented(tmp_path):
    ft = _write_full_table(tmp_path / "ft.tsv", [
        ("g1", "Complete", "ctg1", 1, 100, "+", 90.0, 99),
        ("g2", "Missing"),
        ("g3", "Fragmented", "ctg1", 1, 100, "+", 90.0, 99),
    ])
    leaves = build_outgroup_only_leaves("outg", _mk_compleasm(ft, n_total=10))
    assert len(leaves) == 1
    assert leaves[0].single_copy_genes() == {"g1"}


# --------------------------------------------------------------------------
# CLI validation for --phylo-outgroup-only
# --------------------------------------------------------------------------
def _touch_fasta(path: Path) -> Path:
    path.write_text(">c1\nACGT\n")
    return path


def _run_cli(argv):
    import dnadis.cli as cli
    old = sys.argv
    sys.argv = ["dnadis.py"] + argv
    try:
        cli.main()
    finally:
        sys.argv = old


def test_phylo_outgroup_only_rejected_in_single_assembly_mode(tmp_path):
    ref = _touch_fasta(tmp_path / "ref.fa")
    qry = _touch_fasta(tmp_path / "q.fa")
    with pytest.raises(SystemExit) as ei:
        _run_cli([
            "-r", str(ref), "-q", str(qry), "-o", str(tmp_path / "out"),
            "--compleasm-lineage", "eukaryota", "--phylo-outgroup-only", "foo",
        ])
    assert "multi-assembly" in str(ei.value)


def test_phylo_outgroup_only_requires_compleasm_lineage(tmp_path):
    ref = _touch_fasta(tmp_path / "ref.fa")
    fofn = tmp_path / "asm.tsv"
    fofn.write_text(f"{_touch_fasta(tmp_path / 'a.fa')}\ta\n{_touch_fasta(tmp_path / 'b.fa')}\tb\n")
    with pytest.raises(SystemExit) as ei:
        _run_cli([
            "-r", str(ref), "--fofn", str(fofn), "-o", str(tmp_path / "out"),
            "--phylo-outgroup-only", "a",
        ])
    assert "compleasm-lineage" in str(ei.value)


def test_phylo_outgroup_only_unknown_name_rejected(tmp_path):
    ref = _touch_fasta(tmp_path / "ref.fa")
    fofn = tmp_path / "asm.tsv"
    fofn.write_text(f"{_touch_fasta(tmp_path / 'a.fa')}\ta\n{_touch_fasta(tmp_path / 'b.fa')}\tb\n")
    with pytest.raises(SystemExit) as ei:
        _run_cli([
            "-r", str(ref), "--fofn", str(fofn), "-o", str(tmp_path / "out"),
            "--compleasm-lineage", "eukaryota", "--phylo-outgroup-only", "nope",
        ])
    assert "not found among assemblies" in str(ei.value)


def test_phylo_outgroup_only_all_assemblies_rejected(tmp_path):
    ref = _touch_fasta(tmp_path / "ref.fa")
    fofn = tmp_path / "asm.tsv"
    fofn.write_text(f"{_touch_fasta(tmp_path / 'a.fa')}\ta\n{_touch_fasta(tmp_path / 'b.fa')}\tb\n")
    with pytest.raises(SystemExit) as ei:
        _run_cli([
            "-r", str(ref), "--fofn", str(fofn), "-o", str(tmp_path / "out"),
            "--compleasm-lineage", "eukaryota",
            "--phylo-outgroup-only", "a", "--phylo-outgroup-only", "b",
        ])
    assert "at least one non-outgroup assembly" in str(ei.value)
