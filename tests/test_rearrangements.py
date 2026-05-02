"""Tests for dnadis.detection.rearrangements.

Focus on the off-target gates and the terminal/whole-arm test that decide
whether marginal protein-mode evidence should produce a translocation call.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnadis.detection.rearrangements import detect_rearrangements


def _macro_row(
    contig: str,
    contig_len: int,
    ref_id: str,
    *,
    qstart: int,
    qend: int,
    union_bp: int,
    n_segments: int,
    gene_count: int,
    strand: str = "+",
    chain_id: int = 1,
    matches: int | None = None,
    aln_len: int | None = None,
    identity: float = 0.99,
    ref_start: int = 0,
    ref_end: int = 0,
) -> tuple:
    """Build a 19-field macro_block tuple matching chain_parsing's output."""
    qspan_bp = qend - qstart
    if matches is None:
        matches = union_bp
    if aln_len is None:
        aln_len = union_bp
    return (
        contig,
        contig_len,
        ref_id,
        ref_id,
        "NA",
        strand,
        chain_id,
        qstart,
        qend,
        qspan_bp,
        union_bp,
        matches,
        aln_len,
        f"{identity:.6f}",
        float(union_bp),
        n_segments,
        gene_count,
        ref_start,
        ref_end,
    )


# ---------------------------------------------------------------------------
# Off-target gate tests
# ---------------------------------------------------------------------------

def test_lm7210_marginal_evidence_no_translocation():
    """The lm7210 chr12→chr1 case: 2 chains, 1+1 genes, ~771 bp matches,
    52 kb genomic span. Should NOT produce any translocation call after
    union_bp/gene/segment gates are applied.
    """
    contig_len = 15_870_478
    rows = [
        # On-target chr12: dominant, many segments and genes
        _macro_row(
            "chr12", contig_len, "chr12",
            qstart=13_485, qend=11_502_071,
            union_bp=1_837_501, n_segments=580, gene_count=577,
        ),
        _macro_row(
            "chr12", contig_len, "chr12",
            qstart=11_765_606, qend=15_814_781, chain_id=2,
            union_bp=774_461, n_segments=153, gene_count=152,
        ),
        # Off-target chr1: two thin chains, 1 gene + 1 gene, ~52 kb span
        _macro_row(
            "chr12", contig_len, "chr1",
            qstart=1_439_213, qend=1_491_777, strand="-",
            union_bp=618, n_segments=2, gene_count=1,
        ),
        _macro_row(
            "chr12", contig_len, "chr1",
            qstart=4_349_123, qend=4_349_276, strand="-", chain_id=2,
            union_bp=153, n_segments=1, gene_count=1,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"chr12": "chr12"},
        ref_lengths={"chr12": 16_000_000, "chr1": 15_000_000},
        query_lengths={"chr12": contig_len},
    )
    # No translocation should fire — gene/segment/union_bp gates reject it.
    assert all(c.partner_ref_id != "chr1" for c in calls)


def test_off_target_below_union_bp_gate_rejected():
    """Off-target span >= 50 kb but aligned content < 10 kb → no call."""
    contig_len = 10_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=0, qend=8_000_000,
            union_bp=2_000_000, n_segments=300, gene_count=300,
        ),
        # Off-target: 60 kb genomic span but only 2 kb aligned content
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=4_000_000, qend=4_060_000,
            union_bp=2_000, n_segments=10, gene_count=5,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 10_000_000, "chrB": 10_000_000},
        query_lengths={"ctg1": contig_len},
    )
    assert all(c.partner_ref_id != "chrB" for c in calls)


def test_off_target_below_gene_floor_rejected():
    """Off-target with sufficient span/segments/union but only 2 genes → no call."""
    contig_len = 10_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=0, qend=8_000_000,
            union_bp=2_000_000, n_segments=300, gene_count=300,
        ),
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=4_000_000, qend=4_080_000,
            union_bp=20_000, n_segments=10, gene_count=2,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 10_000_000, "chrB": 10_000_000},
        query_lengths={"ctg1": contig_len},
    )
    assert all(c.partner_ref_id != "chrB" for c in calls)


def test_off_target_below_segment_floor_rejected():
    """Off-target with sufficient span/genes/union but only 2 segments → no call."""
    contig_len = 10_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=0, qend=8_000_000,
            union_bp=2_000_000, n_segments=300, gene_count=300,
        ),
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=4_000_000, qend=4_080_000,
            union_bp=20_000, n_segments=2, gene_count=5,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 10_000_000, "chrB": 10_000_000},
        query_lengths={"ctg1": contig_len},
    )
    assert all(c.partner_ref_id != "chrB" for c in calls)


# ---------------------------------------------------------------------------
# Terminal/whole-arm tests
# ---------------------------------------------------------------------------

def test_off_target_deep_inside_contig_not_whole_arm():
    """Substantial off-target evidence sitting deep inside the contig
    (well past the 100 kb terminal threshold) must be classified as a
    generic 'translocation', not 'whole_arm_translocation'.
    """
    contig_len = 16_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=0, qend=15_000_000,
            union_bp=4_000_000, n_segments=500, gene_count=400,
        ),
        # Off-target block at 1.4–1.5 Mb — > 100 kb from the 5' terminus.
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=1_400_000, qend=1_460_000,
            union_bp=20_000, n_segments=10, gene_count=5,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 16_000_000, "chrB": 15_000_000},
        query_lengths={"ctg1": contig_len},
    )
    chrB_calls = [c for c in calls if c.partner_ref_id == "chrB"]
    assert chrB_calls, "expected a translocation call for chrB"
    assert all(c.rearrangement_type == "translocation" for c in chrB_calls)


def test_off_target_at_5prime_terminus_is_whole_arm():
    """Off-target evidence reaching within 100 kb of the 5' end → whole_arm."""
    contig_len = 16_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=2_000_000, qend=15_000_000,
            union_bp=4_000_000, n_segments=500, gene_count=400,
        ),
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=10_000, qend=1_500_000,
            union_bp=400_000, n_segments=80, gene_count=50,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 16_000_000, "chrB": 15_000_000},
        query_lengths={"ctg1": contig_len},
    )
    chrB_calls = [c for c in calls if c.partner_ref_id == "chrB"]
    assert chrB_calls
    assert any(c.rearrangement_type == "whole_arm_translocation" for c in chrB_calls)
    wa = next(c for c in chrB_calls if c.rearrangement_type == "whole_arm_translocation")
    assert "5'" in wa.evidence


def test_off_target_at_3prime_terminus_is_whole_arm():
    """Off-target evidence reaching within 100 kb of the 3' end → whole_arm."""
    contig_len = 16_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=0, qend=14_000_000,
            union_bp=4_000_000, n_segments=500, gene_count=400,
        ),
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=14_500_000, qend=15_990_000,
            union_bp=400_000, n_segments=80, gene_count=50,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 16_000_000, "chrB": 15_000_000},
        query_lengths={"ctg1": contig_len},
    )
    chrB_calls = [c for c in calls if c.partner_ref_id == "chrB"]
    assert chrB_calls
    wa = [c for c in chrB_calls if c.rearrangement_type == "whole_arm_translocation"]
    assert wa and "3'" in wa[0].evidence


def test_nucleotide_mode_skips_gene_and_segment_gates():
    """Nucleotide mode never has gene IDs (gene_count=0) and a perfect
    full-length alignment can produce a single segment.  Gene/segment gates
    must be skipped — only span and union_bp gates apply.
    """
    contig_len = 16_000_000
    rows = [
        (
            "ctg1", contig_len, "chrA", "chrA", "NA", "+", 1,
            0, 14_000_000, 14_000_000, 4_000_000,
            4_000_000, 4_000_000, "0.990000", 4_000_000.0,
            1, 0,  # n_segments=1, gene_count=0 — typical nucleotide mode
            0, 14_000_000,
        ),
        (
            "ctg1", contig_len, "chrB", "chrB", "NA", "+", 1,
            14_500_000, 15_990_000, 1_490_000, 600_000,
            600_000, 600_000, "0.990000", 600_000.0,
            1, 0,  # one segment, zero genes
            0, 1_490_000,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 16_000_000, "chrB": 15_000_000},
        query_lengths={"ctg1": contig_len},
        synteny_mode="nucleotide",
    )
    chrB_calls = [c for c in calls if c.partner_ref_id == "chrB"]
    assert chrB_calls
    assert any(c.rearrangement_type == "whole_arm_translocation" for c in chrB_calls)


def test_nucleotide_mode_still_enforces_union_bp_gate():
    """In nucleotide mode the union_bp gate must still reject thin evidence."""
    contig_len = 10_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=0, qend=8_000_000,
            union_bp=2_000_000, n_segments=10, gene_count=0,
        ),
        # 60 kb genomic span, only 2 kb aligned — too thin even in nt mode.
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=4_000_000, qend=4_060_000,
            union_bp=2_000, n_segments=1, gene_count=0,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 10_000_000, "chrB": 10_000_000},
        query_lengths={"ctg1": contig_len},
        synteny_mode="nucleotide",
    )
    assert all(c.partner_ref_id != "chrB" for c in calls)


def test_off_target_spanning_both_ends_not_terminal():
    """Off-target evidence touching both ends should not be 'terminal' (XOR)."""
    contig_len = 16_000_000
    rows = [
        _macro_row(
            "ctg1", contig_len, "chrA",
            qstart=2_000_000, qend=14_000_000,
            union_bp=4_000_000, n_segments=500, gene_count=400,
        ),
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=10_000, qend=500_000,
            union_bp=200_000, n_segments=50, gene_count=30,
        ),
        _macro_row(
            "ctg1", contig_len, "chrB",
            qstart=15_500_000, qend=15_990_000, chain_id=2,
            union_bp=200_000, n_segments=50, gene_count=30,
        ),
    ]
    calls = detect_rearrangements(
        macro_block_rows=rows,
        best_ref={"ctg1": "chrA"},
        ref_lengths={"chrA": 16_000_000, "chrB": 15_000_000},
        query_lengths={"ctg1": contig_len},
    )
    chrB_calls = [c for c in calls if c.partner_ref_id == "chrB"]
    assert chrB_calls
    # Either a generic translocation or fusion, but never whole_arm.
    assert all(c.rearrangement_type != "whole_arm_translocation" for c in chrB_calls)
