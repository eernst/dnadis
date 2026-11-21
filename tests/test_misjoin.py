"""Tests for ``dnadis misjoin`` (junctions, sister evidence, calls, breaking)."""
from __future__ import annotations

import csv
import gzip
import random
import shutil
from pathlib import Path

import pytest

from dnadis.misjoin import breaker
from dnadis.misjoin import identify as idf
from dnadis.misjoin import reads as rd
from dnadis.misjoin import sisters as sis
from dnadis.misjoin.junctions import (
    Block,
    Junction,
    Side,
    build_runs,
    chromosome_number,
    find_junctions,
    find_partners,
    read_macro_blocks,
)

MB_COLS = [
    "contig", "contig_len", "ref_id", "chrom_id", "subgenome", "strand", "chain_id",
    "qstart", "qend", "qspan_bp", "union_bp", "matches", "aln_len", "identity", "score",
    "n_segments", "gene_count_chain", "ref_start", "ref_end",
]


def blk(contig, ref, strand, q0, q1, r0, r1, union=None, length=10_000_000, sg="A"):
    return Block(contig, length, ref, sg, strand, q0, q1, r0, r1,
                 union if union is not None else q1 - q0)


def write_macro_blocks(path: Path, blocks):
    with path.open("w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MB_COLS)
        for b in blocks:
            w.writerow([b.contig, b.contig_len, b.ref_id, b.ref_id[:-1], b.subgenome, b.strand, 1,
                        b.qstart, b.qend, b.qend - b.qstart, b.union_bp, 0, 0, 0.9, 0, 1, 0,
                        b.ref_start, b.ref_end])


# ---------------------------------------------------------------------------
# Junctions
# ---------------------------------------------------------------------------
def test_switch_gives_one_junction_with_reference_positions():
    blocks = [
        blk("c1", "chr10A", "-", 0, 1_000_000, 6_000_000, 7_000_000),
        blk("c1", "chr10A", "+", 1_100_000, 1_400_000, 4_700_000, 5_000_000),
        blk("c1", "chr15A", "-", 1_600_000, 2_000_000, 8_400_000, 8_800_000),
        blk("c1", "chr15A", "+", 2_100_000, 3_000_000, 9_000_000, 9_900_000),
    ]
    js = find_junctions({"c1": blocks})
    assert len(js) == 1
    j = js[0]
    assert (j.start, j.end) == (1_400_000, 1_600_000)
    assert j.left.ref_id == "chr10A" and j.right.ref_id == "chr15A"
    # + block ending at the junction: its ref_end; - block starting there: its ref_end
    assert j.left.ref_pos == 5_000_000 and j.left.strand == "+"
    assert j.right.ref_pos == 8_800_000 and j.right.strand == "-"
    assert j.left.aligned_bp == 1_300_000
    assert (j.left.qstart, j.right.qend) == (0, 3_000_000)
    assert j.junction_id == "c1:j1"


def test_small_foreign_blocks_do_not_split_a_run():
    blocks = [
        blk("c1", "chr1A", "+", 0, 500_000, 0, 500_000),
        blk("c1", "chr13A", "-", 520_000, 530_000, 1_260_000, 1_270_000),   # repeat copy
        blk("c1", "chr7T", "+", 540_000, 570_000, 3_000_000, 3_030_000),     # < min_side_bp
        blk("c1", "chr1A", "+", 600_000, 1_200_000, 600_000, 1_200_000),
    ]
    runs = build_runs(blocks, min_block_bp=10_000, min_side_bp=50_000)
    assert [r.ref_id for r in runs] == ["chr1A"]
    assert find_junctions({"c1": blocks}) == []


def test_short_side_below_min_side_bp_is_not_a_junction():
    blocks = [
        blk("c1", "chr5A", "+", 0, 2_600_000, 18_000_000, 20_600_000),
        blk("c1", "chr11A", "+", 2_700_000, 2_740_000, 700_000, 740_000),
    ]
    assert find_junctions({"c1": blocks}, min_side_bp=50_000) == []
    assert len(find_junctions({"c1": blocks}, min_side_bp=30_000)) == 1


def test_insertion_gives_two_junctions():
    blocks = [
        blk("c1", "chr16A", "-", 0, 800_000, 5_900_000, 6_700_000),
        blk("c1", "chr14A", "+", 900_000, 980_000, 7_600_000, 7_680_000),
        blk("c1", "chr16A", "+", 1_100_000, 2_000_000, 5_000_000, 5_900_000),
    ]
    js = find_junctions({"c1": blocks})
    assert [(j.left.ref_id, j.right.ref_id) for j in js] == [("chr16A", "chr14A"), ("chr14A", "chr16A")]
    assert [j.index for j in js] == [1, 2]


def test_reciprocal_partners():
    a = [blk("a", "chr15A", "-", 0, 2_400_000, 2_000_000, 4_400_000),
         blk("a", "chr14A", "-", 2_500_000, 8_000_000, 5_000_000, 12_570_000)]
    b = [blk("b", "chr15A", "+", 0, 600_000, 1_380_000, 1_980_000),
         blk("b", "chr14A", "+", 700_000, 1_000_000, 13_030_000, 13_330_000)]
    c = [blk("c", "chr15A", "+", 0, 600_000, 9_000_000, 9_600_000),
         blk("c", "chr14A", "+", 700_000, 1_000_000, 1_000_000, 1_300_000)]
    js = find_junctions({"a": a, "b": b, "c": c})
    p = find_partners(js, max_ref_dist=1_000_000)
    assert p == {"a:j1": "b:j1", "b:j1": "a:j1"}


def test_organelle_blocks_are_ignored(tmp_path):
    path = tmp_path / "x.macro_blocks.tsv"
    write_macro_blocks(path, [
        blk("c1", "chrC", "+", 0, 76_000, 100_000, 176_000),
        blk("c1", "chr8P", "+", 80_000, 150_000, 20_000_000, 20_070_000, sg="P"),
    ])
    blocks = read_macro_blocks(path)
    assert [b.ref_id for b in blocks["c1"]] == ["chr8P"]
    assert find_junctions(blocks) == []


def test_unplaced_scaffold_blocks_are_transparent(tmp_path):
    from dnadis.misjoin.junctions import read_reference_lengths
    ref = tmp_path / "reference"
    ref.mkdir()
    (ref / "reference.ref_lengths.tsv").write_text(
        "ref_id\tchrom_id\tsubgenome\tref_len\tkind\n"
        "chr1A\tchr1\tA\t30000000\tnuclear\n"
        "chr2A\tchr2\tA\t30000000\tnuclear\n"
        "scaf_47\tscaf_47\tNA\t2000000\tother\n"
        "scaf_59\tscaf_59\tNA\t2000000\tother\n")
    lengths, nuclear = read_reference_lengths(tmp_path)
    assert nuclear == {"chr1A", "chr2A"} and lengths["scaf_47"] == 2_000_000
    path = tmp_path / "x.macro_blocks.tsv"
    write_macro_blocks(path, [
        blk("c1", "chr1A", "+", 0, 1_000_000, 0, 1_000_000),
        blk("c1", "scaf_47", "+", 1_000_000, 1_400_000, 0, 400_000),     # knob array
        blk("c1", "scaf_59", "-", 1_400_000, 1_800_000, 0, 400_000),
        blk("c1", "chr1A", "+", 1_800_000, 3_000_000, 1_100_000, 2_300_000),
        blk("c1", "scaf_59", "+", 3_000_000, 3_300_000, 500_000, 800_000),
        blk("c1", "chr2A", "+", 3_300_000, 4_500_000, 5_000_000, 6_200_000),
    ])
    js = find_junctions(read_macro_blocks(path, nuclear))
    assert [(j.left.ref_id, j.right.ref_id) for j in js] == [("chr1A", "chr2A")]
    assert (js[0].start, js[0].end) == (3_000_000, 3_300_000)
    # without the kind column every non-organelle block counts
    assert len(find_junctions(read_macro_blocks(path))) == 5


def test_chromosome_number():
    assert chromosome_number("chr05A") == "5"
    assert chromosome_number("Chr12_P") == "12"
    assert chromosome_number("chrUn") is None


# ---------------------------------------------------------------------------
# Sister evidence
# ---------------------------------------------------------------------------
def _junction():
    blocks = [
        blk("q", "chr10A", "+", 0, 2_600_000, 2_000_000, 4_600_000, length=5_600_000),
        blk("q", "chr15A", "+", 3_000_000, 5_600_000, 8_000_000, 10_600_000, length=5_600_000),
    ]
    (j,) = find_junctions({"q": blocks})
    return j


def hit(r0, r1, sister, s0, s1, slen, strand="+", ident=0.99, primary=True):
    return sis.Hit(r0, r1, sister, s0, s1, slen, strand, ident, primary)


P = sis.SisterParams()


def test_sister_carrying_the_join():
    j = _junction()
    # one sister contig collinear through the junction (50 kb insertion in the sister)
    hits = [
        hit(2_300_000, 2_950_000, "s1", 1_300_000, 1_950_000, 9_000_000),
        hit(2_950_000, 3_300_000, "s1", 2_000_000, 2_350_000, 9_000_000),
    ]
    ev = sis.evaluate_unit(j, hits, "sisA", {"s1"}, P)
    assert ev.status == "carries"


def test_same_genotype_carrier_does_not_count():
    j = _junction()
    hits = [hit(2_300_000, 3_300_000, "s1", 1_300_000, 2_300_000, 9_000_000, ident=1.0)]
    ev = sis.evaluate_unit(j, hits, "sisA", {"s1"}, P)
    assert ev.status == "carries_same_genotype"


def test_sister_contradicting_the_join_refines_the_break():
    j = _junction()
    # left side continues on sister contig s1 past the point where the
    # junction contig switches; the right side is on s2.
    hits = [
        hit(2_300_000, 2_941_000, "s1", 4_000_000, 4_641_000, 20_000_000),
        hit(2_959_000, 3_300_000, "s2", 7_000_000, 7_341_000, 30_000_000),
    ]
    ev = sis.evaluate_unit(j, hits, "sisA", {"s1", "s2"}, P)
    assert ev.status == "contradicts"
    assert ev.left_continues and ev.right_continues
    assert (ev.x_left, ev.x_right) == (2_941_000, 2_959_000)


def test_sister_contradicting_on_reverse_strand():
    j = _junction()
    # s1 aligned in reverse: contig 2.3-2.94 Mb = s1 5.0-4.36 Mb; s1 continues below 4.36 Mb
    hits = [
        hit(2_300_000, 2_940_000, "s1", 4_360_000, 5_000_000, 5_100_000, strand="-"),
        hit(2_960_000, 3_300_000, "s2", 0, 340_000, 20_000_000),
    ]
    ev = sis.evaluate_unit(j, hits, "sisA", {"s1", "s2"}, P)
    assert ev.status == "contradicts"
    assert ev.left_continues            # 4.36 Mb of s1 before the aligned part
    assert not ev.right_continues       # s2 starts at the junction


def test_sister_contigs_ending_at_the_junction():
    j = _junction()
    hits = [
        hit(2_300_000, 2_940_000, "s1", 0, 640_000, 650_000),
        hit(2_960_000, 3_300_000, "s2", 10_000, 350_000, 2_000_000),
    ]
    ev = sis.evaluate_unit(j, hits, "sisA", {"s1", "s2"}, P)
    assert ev.status == "ends"


def test_hits_outside_the_unit_or_below_identity_are_ignored():
    j = _junction()
    hits = [
        hit(2_300_000, 3_300_000, "other", 0, 1_000_000, 9_000_000),
        hit(2_300_000, 3_300_000, "s1", 0, 1_000_000, 9_000_000, ident=0.8),
    ]
    assert sis.evaluate_unit(j, hits, "sisA", {"s1"}, P).status == "absent"


def test_shared_segment_gives_two_breaks():
    j = _junction()
    # both sides' sister chains run through a shared 11 kb segment
    units = []
    for k in range(3):
        units.append(sis.UnitEvidence(f"u{k}", "contradicts", "a", "b",
                                      x_left=2_959_400 + k, x_right=2_948_200 - k,
                                      left_continues=True, right_continues=True))
    breaks, core = idf.locate_breaks(j, units, shared_min_bp=1_000)
    # the furthest left chain end and the earliest right chain start
    assert core == ("shared", 2_948_198, 2_959_402)
    assert [(b.pos, b.kind) for b in breaks] == [
        (2_948_198, "shared_left_edge"), (2_959_402, "shared_right_edge")]


def test_gap_gives_one_break_between_the_chains():
    j = _junction()
    units = [sis.UnitEvidence("u", "contradicts", "a", "b", x_left=2_941_000,
                              x_right=2_959_000, left_continues=True)]
    breaks, core = idf.locate_breaks(j, units, shared_min_bp=1_000)
    assert [(b.pos, b.lo, b.hi, b.kind) for b in breaks] == [(2_950_000, 2_941_000, 2_959_000, "single")]
    assert core == ("gap", 2_941_000, 2_959_000)


def test_no_sister_chain_falls_back_to_reference_interval():
    j = _junction()
    breaks, core = idf.locate_breaks(j, [], shared_min_bp=1_000)
    assert [(b.pos, b.lo, b.hi) for b in breaks] == [(2_800_000, 2_600_000, 3_000_000)]
    assert core is None


def _units(statuses):
    return [sis.UnitEvidence(f"u{i}", s, "a", "b", x_left=2_941_000, x_right=2_959_000,
                             left_continues=True, right_continues=True)
            for i, s in enumerate(statuses)]


def test_decide_carrier_keeps_the_join():
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 4 + ["carries"]), idf.IdentifyParams(), [])
    assert (jc.call, jc.action) == ("sister_supported", "keep")


def _units_sg(spec):
    """Units from (status, same_genotype) pairs."""
    us = _units([s for s, _ in spec])
    for u, (_, sg) in zip(us, spec):
        u.same_genotype = sg
    return us


def test_decide_genotype_specific_when_same_genotype_assemblies_agree():
    p = idf.IdentifyParams()
    units = _units_sg([("carries_same_genotype", True)] * 3 + [("contradicts", False)] * 8)
    jc = idf.decide("asm", _junction(), units, p, [])
    assert (jc.call, jc.action) == ("genotype_specific", "review")
    # one carrier of the same genotype could be built from the target's contigs
    units = _units_sg([("carries_same_genotype", True)] + [("contradicts", False)] * 8)
    assert idf.decide("asm", _junction(), units, p, []).call == "misjoin"


def test_decide_same_genotype_contradiction_restores_the_misjoin():
    units = _units_sg([("carries_same_genotype", True)] * 4 + [("contradicts", True)]
                      + [("contradicts", False)] * 7)
    jc = idf.decide("asm", _junction(), units, idf.IdentifyParams(), [])
    assert (jc.call, jc.action, jc.confidence) == ("misjoin", "break", "high")


def test_decide_units_that_end_lower_confidence():
    p = idf.IdentifyParams()
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 9 + ["ends"]), p, [])
    assert jc.confidence == "high"
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 8 + ["ends"] * 2), p, [])
    assert (jc.call, jc.confidence) == ("misjoin", "medium")


def test_decide_sparse_sides_are_reviewed():
    blocks = [  # two runs of scattered hits, each ~4% aligned
        blk("q", "chr8A", "+", 0, 100_000, 1_000, 101_000, length=12_000_000),
        blk("q", "chr8A", "+", 2_900_000, 3_000_000, 2_000, 102_000, length=12_000_000),
        blk("q", "chr9A", "+", 3_100_000, 3_160_000, 5_000, 65_000, length=12_000_000),
        blk("q", "chr9A", "+", 4_500_000, 4_560_000, 9_000, 69_000, length=12_000_000),
    ]
    (j,) = find_junctions({"q": blocks})
    units = [sis.UnitEvidence(f"u{i}", "contradicts", x_left=3_050_000, x_right=3_050_000,
                              left_continues=True, right_continues=True) for i in range(4)]
    jc = idf.decide("asm", j, units, idf.IdentifyParams(), [])
    assert (jc.call, jc.action) == ("misjoin", "review") and "sparse_sides" in jc.flags
    # one dense side is enough
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 4), idf.IdentifyParams(), [])
    assert jc.action == "break" and "sparse_sides" not in jc.flags


def _call(contig, index, left, right, lq, rq, call="misjoin", flags=()):
    side = lambda ref, q: Side(ref, "A", 0, "+", q[1] - q[0], q[0], q[1], 0, 1)  # noqa: E731
    j = Junction(contig, 5_000_000, index, side(left, lq), side(right, rq), lq[1], rq[0])
    action = "break" if call in ("misjoin", "repeat_misjoin") else "keep"
    return idf.JunctionCall("asm", j, call, action, "high", [idf.Break(lq[1], lq[1], lq[1], "single")],
                            flags=list(flags))


def test_short_run_between_one_chromosome_is_an_insertion():
    a = _call("c", 1, "chr1A", "chr2A", (0, 2_000_000), (2_010_000, 2_066_000))
    b = _call("c", 2, "chr2A", "chr1A", (2_010_000, 2_066_000), (2_070_000, 5_000_000))
    idf.mark_short_inserts([a, b], 300_000)
    assert [(x.call, x.action) for x in (a, b)] == [("insertion", "review")] * 2
    assert all("short_insert" in x.flags for x in (a, b))
    # a long run, or one between two chromosomes, is left alone
    a = _call("c", 1, "chr1A", "chr2A", (0, 2_000_000), (2_010_000, 2_500_000))
    b = _call("c", 2, "chr2A", "chr1A", (2_010_000, 2_500_000), (2_510_000, 5_000_000))
    c = _call("c", 3, "chr1A", "chr3A", (2_510_000, 4_000_000), (4_010_000, 4_050_000))
    d = _call("c", 4, "chr3A", "chr4A", (4_010_000, 4_050_000), (4_060_000, 5_000_000))
    idf.mark_short_inserts([a, b, c, d], 300_000)
    assert [x.call for x in (a, b, c, d)] == ["misjoin"] * 4
    # a scaffold join keeps its break: cutting both gaps lifts the contig out
    a = _call("c", 1, "chr1A", "chr2A", (0, 2_000_000), (2_010_000, 2_066_000), flags=["scaffold_gap"])
    b = _call("c", 2, "chr2A", "chr1A", (2_010_000, 2_066_000), (2_070_000, 5_000_000), flags=["scaffold_gap"])
    idf.mark_short_inserts([a, b], 300_000)
    assert [x.action for x in (a, b)] == ["break", "break"]


def _sister(tmp_path, name, contigs, lengths):
    s = sis.SisterAssembly(name, tmp_path / f"{name}.fa")
    for cname, ref, sg, group, length in contigs:
        s.contigs[cname] = sis.SisterContig(cname, sg, group, length, ref)
    s.ref_lengths = lengths
    return s


def test_sister_units_need_the_junction_chromosomes(tmp_path):
    lengths = {"chr10A": 20_000_000, "chr15A": 20_000_000}
    full = _sister(tmp_path, "full", [("a", "chr10A", "A", "", 19_000_000),
                                      ("b", "chr15A", "A", "", 12_000_000)], lengths)
    empty = _sister(tmp_path, "empty", [("a", "chr10A", "A", "", 300_000),
                                        ("b", "chr15A", "A", "", 2_000_000)], lengths)
    other_ref = _sister(tmp_path, "other", [("a", "Chr10", "A", "", 300_000)], {"Chr10": 20_000_000})
    derived = _sister(tmp_path, "scaf", [("a", "chr10A", "A", "", 19_000_000),
                                         ("b", "chr15A", "A", "", 19_000_000)], lengths)
    asms = {s.name: s for s in (full, empty, other_ref, derived)}
    units = idf._sister_units("asm", _junction(), asms, {}, sis.SisterParams(),
                              {"asm": {"scaf"}, "scaf": {"asm"}})
    # untestable sides (another reference's names) do not exclude a unit
    assert [u.unit for u in units] == ["full", "other"]
    assert full.covers("full", [("chr10A", "A"), ("chr15P", "P")], 0.5)   # P not held


def test_query_files_are_keyed_by_region_set(tmp_path):
    s = _sister(tmp_path, "sis", [("a", "chr1A", "A", "", 1_000)], {})
    seqs = {("t1", "c"): "ACGT" * 100, ("t2", "c"): "TTGA" * 100}
    r1 = sis.Region("t1|c|0", "t1", "c", 0, 100, ("A",))
    r2 = sis.Region("t2|c|0", "t2", "c", 0, 100, ("A",))
    (_, q1, paf1), = idf._query_files({"sis": s}, [r1], seqs, tmp_path, "regions")
    (_, q2, paf2), = idf._query_files({"sis": s}, [r2], seqs, tmp_path, "regions")
    # another run's targets: another query and another PAF for the same sister
    assert q1 != q2 and paf1 != paf2
    assert idf._query_files({"sis": s}, [r1], seqs, tmp_path, "regions")[0][2] == paf1


def test_decide_same_genotype_carrier_does_not_keep():
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 3 + ["carries_same_genotype"]), idf.IdentifyParams(), [])
    assert (jc.call, jc.action, jc.confidence) == ("misjoin", "break", "high")
    assert [b.pos for b in jc.breaks] == [2_950_000]


def test_decide_repeat_between_the_sides_breaks_both_edges():
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 2), idf.IdentifyParams(), [(2_941_000, 2_959_000, 30)])
    assert jc.call == "repeat_misjoin"
    assert [(b.pos, b.kind) for b in jc.breaks] == [
        (2_941_000, "shared_left_edge"), (2_959_000, "shared_right_edge")]
    rows = idf.call_rows(jc)
    assert [r["break_id"] for r in rows] == ["q:j1a", "q:j1b"]
    assert rows[0]["core_copies"] == "30"


def test_decide_repeat_is_cut_out_but_not_extended():
    # sister chains leave 1.7 kb between the sides, inside a repeat copy:
    # the repeat labels that segment, it does not move the breaks
    units = [sis.UnitEvidence("u", "contradicts", "a", "b", x_left=2_496_212, x_right=2_497_956,
                              left_continues=True, right_continues=True)] * 3
    jc = idf.decide("asm", _junction(), units, idf.IdentifyParams(),
                    [(2_496_200, 2_507_300, 28), (2_700_000, 2_710_000, 5)])
    assert jc.call == "repeat_misjoin"
    assert jc.core == ("repeat", 2_496_212, 2_497_956) and jc.core_copies == 28
    assert [b.pos for b in jc.breaks] == [2_496_212, 2_497_956]


def test_decide_partner_shared_segment_sets_both_edges():
    units = [sis.UnitEvidence("u", "contradicts", "a", "b", x_left=2_496_212, x_right=2_497_956,
                              left_continues=True, right_continues=True)] * 3
    jc = idf.decide("asm", _junction(), units, idf.IdentifyParams(),
                    [(2_496_200, 2_505_800, 28)], shared=(2_496_210, 2_507_386))
    assert jc.core == ("partner_shared", 2_496_210, 2_507_386) and jc.core_copies == 28
    assert [b.pos for b in jc.breaks] == [2_496_210, 2_507_386]


def test_decide_scaffold_gap_takes_the_break_and_reads_do_not_count():
    units = _units(["contradicts"] * 3)
    jc = idf.decide("asm", _junction(), units, idf.IdentifyParams(),
                    gaps=[(100, 200), (2_950_100, 2_950_200)])
    assert [(b.pos, b.lo, b.hi, b.kind) for b in jc.breaks] == [
        (2_950_100, 2_950_100, 2_950_200, "scaffold_gap")]
    assert "scaffold_gap" in jc.flags
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=0, spanning=0)], rd.ReadParams())
    assert jc.action == "break"
    jc2 = idf.decide("asm", _junction(), units, idf.IdentifyParams(), gaps=[(2_950_100, 2_950_200)])
    idf.apply_reads(jc2, [rd.ReadEvidence(split_reads=20, split_pos=2_950_100, spanning=0)],
                    rd.ReadParams())
    assert "reads_split" not in jc2.flags and jc2.breaks[0].pos == 2_950_100


def test_decide_spanning_sisters_make_an_insertion():
    units = _units(["spans", "spans", "contradicts"])
    jc = idf.decide("asm", _junction(), units, idf.IdentifyParams(), [])
    assert (jc.call, jc.action) == ("insertion", "review")


def test_decide_without_sisters_is_reference_only():
    jc = idf.decide("asm", _junction(), [], idf.IdentifyParams(), [])
    assert (jc.call, jc.action) == ("reference_only", "review")
    assert idf.call_rows(jc)[0]["evidence"] == "reference only"


def test_decide_single_contradictor_is_unresolved():
    jc = idf.decide("asm", _junction(), _units(["contradicts", "ends", "absent"]), idf.IdentifyParams(), [])
    assert (jc.call, jc.action) == ("unresolved", "review")


def test_copy_profile_and_repeat_runs():
    # probe contig 0-30 kb aligned to three sisters: once over its whole
    # length, plus 30 copies of 10,000-21,200 in two of them
    p = sis.SisterParams()
    own = hit(0, 30_000, "c0", 5_000_000, 5_030_000, 9_000_000)
    copies = [hit(10_000, 21_200, f"c{k}", k * 100_000, k * 100_000 + 11_200, 9_000_000,
                  primary=False) for k in range(1, 31)]
    hits = {"s1": [own] + copies, "s2": [own] + copies, "s3": [own]}
    prof = sis.copy_profile(hits, 0, 30_000, p)
    assert prof[0] == (0, 1) and prof[150] == (15_000, 31)
    assert sis.repeat_runs(prof, 3, 1_000, p.profile_bin) == [(10_000, 21_200, 31)]
    hits["s2"] = [own]
    assert sis.repeat_runs(sis.copy_profile(hits, 0, 30_000, p), 3, 1_000, 100) == []
    # low-identity copies do not count
    hits = {"s1": [own] + [hit(10_000, 21_200, "x", 0, 11_200, 9e6, ident=0.9)] * 5}
    assert sis.repeat_runs(sis.copy_profile(hits, 0, 30_000, p), 3, 1_000, 100) == []


def test_sister_spanning_an_inserted_run():
    # contig: chr16A, 80 kb of chr14A, chr16A; the sister's chr16A runs across
    blocks = [
        blk("q", "chr16A", "+", 0, 800_000, 5_000_000, 5_800_000, length=2_000_000),
        blk("q", "chr14A", "+", 820_000, 900_000, 7_600_000, 7_680_000, length=2_000_000),
        blk("q", "chr16A", "+", 920_000, 2_000_000, 5_800_000, 6_880_000, length=2_000_000),
    ]
    j1, _ = find_junctions({"q": blocks})
    hits = [
        hit(500_000, 810_000, "s16", 1_000_000, 1_310_000, 9_000_000),
        hit(905_000, 1_200_000, "s16", 1_312_000, 1_607_000, 9_000_000),
        hit(825_000, 895_000, "s14", 3_000_000, 3_070_000, 9_000_000),
    ]
    ev = sis.evaluate_unit(j1, hits, "u", {"s16", "s14"}, P)
    assert ev.status == "spans"


def test_nested_runs_have_an_empty_flank():
    j = _junction()
    j.right.qstart, j.right.qend = 2_000_000, 2_500_000     # right run inside the left one
    (l0, l1), (r0, r1) = sis.junction_flanks(j, 300_000)
    assert r1 <= r0


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------
def test_apply_reads_overrides_action():
    p = rd.ReadParams()
    jc = idf.decide("asm", _junction(), [], idf.IdentifyParams(), [])
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=12, split_pos=2_951_755, spanning=0)], p)
    assert (jc.action, jc.confidence) == ("break", "high")
    assert jc.breaks[0].pos == 2_951_755

    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 3), idf.IdentifyParams(), [])
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=0, spanning=20)], p)
    assert jc.action == "keep" and "reads_span" in jc.flags

    # reads that both split and span: a sister-called break goes to review,
    # the position stays, and a keep stays a keep
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 3), idf.IdentifyParams(), [])
    assert jc.action == "break"
    pos = jc.breaks[0].pos
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=3, split_pos=pos + 6_000, spanning=2)], p)
    assert (jc.action, jc.confidence) == ("review", "low") and "reads_mixed" in jc.flags
    assert jc.breaks[0].pos == pos

    jc = idf.decide("asm", _junction(), _units(["carries"]), idf.IdentifyParams(), [])
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=15, split_pos=1, spanning=19)], p)
    assert jc.action == "keep" and "reads_mixed" in jc.flags


def _array_contig(seed=7):
    """30 kb unique, a 40,160 bp array of a 251 bp unit (1% substitutions), 30 kb unique."""
    rng = random.Random(seed)
    unit = "".join(rng.choice("ACGT") for _ in range(251))
    arr = "".join(c if rng.random() > 0.01 else rng.choice("ACGT") for c in unit * 160)
    left = "".join(rng.choice("ACGT") for _ in range(30_000))
    right = "".join(rng.choice("ACGT") for _ in range(30_000))
    return left + arr + right


def test_tandem_run_finds_array():
    seq = _array_contig()
    p = rd.ReadParams()
    lo, hi = rd.tandem_run(seq, 50_000, p)
    assert 29_000 <= lo <= 31_000 and 69_000 <= hi <= 71_000
    assert rd.tandem_run(seq, 31_500, p) == (lo, hi)      # near the edge
    assert rd.tandem_run(seq, 10_000, p) is None
    assert rd.tandem_run(seq, 85_000, p) is None


def test_mark_tandem_flags_break():
    seq = _array_contig()
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 3), idf.IdentifyParams(), [])
    jc.breaks[0].pos = 50_000
    idf.mark_tandem(jc, seq, rd.ReadParams())
    assert jc.breaks[0].tandem is not None and "tandem_repeat" in jc.flags
    row = idf.call_rows(jc)[0]
    assert row["tandem_start"] and "tandem array" in row["evidence"]


def test_telomere_arrays_merge_runs_per_strand():
    from dnadis.detection.telomere import telomere_arrays
    fwd = "TTTAGGG" * 50 + "TTCAGGG" + "TTTAGGG" * 50     # one variant repeat
    rev = "CCCTAAA" * 20
    seq = "ACGT" * 250 + fwd + "ACGT" * 500 + rev
    arrays = telomere_arrays(seq, offset=10)
    assert [(a[0], a[1], a[2], a[3]) for a in arrays] == [
        (1_010, 1_010 + len(fwd), 700, "+"),
        (1_010 + len(fwd) + 2_000, 1_010 + len(fwd) + 2_000 + 140, 140, "-"),
    ]


def _telomere_contig(arrays):
    """A 5.6 Mb contig without telomere repeats, with (position, strand, copies) arrays."""
    seq = list("ACGT" * 1_400_000)
    for pos, strand, n in arrays:
        unit = "TTTAGGG" if strand == "+" else "CCCTAAA"
        seq[pos:pos + 7 * n] = unit * n
    return "".join(seq)


@pytest.mark.parametrize("arrays, side, pos", [
    ([(2_945_000, "-", 300)], "right", 2_945_000),            # right side starts with a telomere
    ([(2_955_000 - 2_100, "+", 300)], "left", 2_955_000),      # left side ends with one
    ([(2_946_000, "+", 300), (2_948_200, "-", 300)], "both", 2_948_100),
])
def test_mark_telomere_at_break(arrays, side, pos):
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 2), idf.IdentifyParams(), [])
    assert (jc.call, jc.confidence) == ("misjoin", "medium")
    idf.mark_telomere(jc, _telomere_contig(arrays), idf.IdentifyParams())
    assert jc.telomere[0] == side and jc.breaks[0].pos == pos
    assert jc.confidence == "high" and "interstitial_telomere" in jc.flags
    row = idf.call_rows(jc)[0]
    assert row["telomere_side"] == side and row["telomere_bp"] in ("2100", "4200")
    assert "telomere array" in row["evidence"]


@pytest.mark.parametrize("arrays", [
    [(2_965_000, "+", 300)],     # right of the break, its chromosome end pointing away
    [(2_930_000, "-", 300)],     # left of the break, pointing away
    [(2_950_000, "-", 100)],     # 700 bp: too short
    [(2_900_000, "-", 300)],     # beyond the window
])
def test_mark_telomere_ignores_arrays_not_at_a_chromosome_end(arrays):
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 2), idf.IdentifyParams(), [])
    idf.mark_telomere(jc, _telomere_contig(arrays), idf.IdentifyParams())
    assert jc.telomere is None and jc.confidence == "medium" and jc.breaks[0].pos == 2_950_000


def test_mark_telomere_keeps_a_sister_supported_join():
    jc = idf.decide("asm", _junction(), _units(["carries"] * 2), idf.IdentifyParams(), [])
    idf.mark_telomere(jc, _telomere_contig([(2_945_000, "-", 300)]), idf.IdentifyParams())
    assert (jc.action, jc.telomere[0]) == ("keep", "right")
    assert "interstitial_telomere" in jc.flags


def test_read_evidence_inside_tandem_array(tmp_path):
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "ctg", "LN": 100_000}]}
    path = tmp_path / "r.bam"
    reads = []
    # reads that lie inside the array (30-70 kb) and cover the break at 50 kb
    for k in range(4):
        reads.append((f"inside{k}", 32_000 + k, [(0, 36_000)]))
    # reads clipped inside the array, as reads from another array would be
    for k in range(5):
        reads.append((f"clip{k}", 55_000 + k, [(4, 8_000), (0, 25_000)]))
    # reads clipped at the array edge, and one read across the whole array
    for k in range(3):
        reads.append((f"edge{k}", 70_000, [(4, 5_000), (0, 10_000)]))
    reads.append(("across", 25_000, [(0, 50_000)]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, cigar in sorted(reads, key=lambda x: x[1]):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 4))
            a.reference_id = 0
            a.reference_start = start
            a.cigartuples = cigar
            a.mapping_quality = 60
            out.write(a)
    pysam.index(str(path))
    bam = rd.open_alignments(path)
    p = rd.ReadParams()
    plain = rd.read_evidence(bam, "ctg", 49_000, 62_000, 50_000, 50_000, p)
    assert (plain.split_reads, plain.split_pos, plain.spanning) == (5, 55_002, 5)
    ev = rd.read_evidence(bam, "ctg", 49_000, 62_000, 50_000, 50_000, p, repeat=(30_000, 70_000))
    assert ev.spanning == 1          # only "across" covers the array plus flank
    ev = rd.read_evidence(bam, "ctg", 49_000, 72_000, 50_000, 50_000, p, repeat=(30_000, 70_000))
    assert (ev.split_reads, ev.split_pos) == (3, 70_000)


def test_read_evidence_from_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "ctg", "LN": 100_000}]}
    path = tmp_path / "r.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        reads = []
        # five reads clipped at 50,000 (aligned left of it), three spanning reads
        for k in range(5):
            reads.append((f"split{k}", 40_000 + k * 10, [(0, 10_000 - k * 10), (4, 3_000)]))
        for k in range(3):
            reads.append((f"span{k}", 45_000 + k, [(0, 10_000)]))
        for name, start, cigar in sorted(reads, key=lambda x: x[1]):
            a = pysam.AlignedSegment()
            a.query_name = name
            qlen = sum(n for op, n in cigar if op in (0, 4))
            a.query_sequence = "A" * qlen
            a.reference_id = 0
            a.reference_start = start
            a.cigartuples = cigar
            a.mapping_quality = 60
            out.write(a)
    pysam.index(str(path))
    bam = rd.open_alignments(path)
    ev = rd.read_evidence(bam, "ctg", 49_000, 51_000, 50_000, 50_000, rd.ReadParams())
    assert (ev.split_reads, ev.split_pos, ev.spanning) == (5, 50_000, 3)
    assert rd.read_evidence(bam, "other", 0, 1, 0, 1, rd.ReadParams()) is None



def _write_bam(path, reads, length=200_000):
    """Indexed BAM of (name, start, cigar) primary alignments at MAPQ 60."""
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "ctg", "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, cigar in sorted(reads, key=lambda x: x[1]):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 4))
            a.reference_id = 0
            a.reference_start = start
            a.cigartuples = cigar
            a.mapping_quality = 60
            out.write(a)
    pysam.index(str(path))
    return rd.open_alignments(path)


def test_tiled_reads_need_every_point_of_the_interval(tmp_path):
    # reads tile 50-150 kb, three deep, except for a stretch at 100-110 kb
    # that only one read runs through
    reads = [(f"r{k}{s}", s + k * 100, [(0, 20_000)])
             for s in range(40_000, 150_000, 10_000) for k in range(3)
             if not 90_000 <= s <= 100_000]
    reads.append(("one", 95_000, [(0, 20_000)]))
    bam = _write_bam(tmp_path / "t.bam", reads)
    p = rd.ReadParams()
    assert rd.tiled_reads(bam, "ctg", 60_000, 80_000, p) >= 3
    assert rd.tiled_reads(bam, "ctg", 60_000, 140_000, p) == 1
    assert rd.tiled_reads(bam, "ctg", 190_000, 195_000, p) == 0
    ev = rd.read_evidence(bam, "ctg", 60_000, 140_000, 60_000, 140_000, p)
    assert ev.spanning == 1
    # the shared-segment rule: one read across all of it
    ev = rd.read_evidence(bam, "ctg", 60_000, 80_000, 60_000, 80_000, p, whole=True)
    assert ev.spanning == 0


def test_split_cluster_with_reads_running_through_is_not_a_break(tmp_path):
    # five reads end at 50 kb (clipped); 12 run through it
    reads = [(f"split{k}", 40_000 + k * 10, [(0, 10_000 - k * 10), (4, 3_000)]) for k in range(5)]
    reads += [(f"through{k}", 44_000 + k * 100, [(0, 12_000)]) for k in range(12)]
    bam = _write_bam(tmp_path / "s.bam", reads)
    p = rd.ReadParams()
    ev = rd.read_evidence(bam, "ctg", 49_000, 51_000, 49_000, 51_000, p)
    assert (ev.split_reads, ev.split_pos, ev.split_through) == (5, 50_000, 12)
    assert not ev.split_supported(p)
    assert rd.ReadEvidence(split_reads=5, split_through=4).split_supported(p)
    # sisters call a misjoin; the reads neither confirm nor refute it
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 3), idf.IdentifyParams(), [])
    jc.breaks[0].lo, jc.breaks[0].hi = 2_900_000, 3_000_000
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=80, split_pos=2_920_000,
                                         split_through=120, spanning=0)], p)
    assert jc.action == "break" and jc.breaks[0].pos == 2_950_000
    assert not {"reads_split", "reads_mixed", "reads_span"} & set(jc.flags)


def test_reads_do_not_break_a_sister_supported_join():
    jc = idf.decide("asm", _junction(), _units(["carries"] * 3), idf.IdentifyParams(), [])
    idf.apply_reads(jc, [rd.ReadEvidence(split_reads=12, split_pos=2_950_000, spanning=0)],
                    rd.ReadParams())
    assert (jc.action, jc.confidence) == ("review", "low") and "reads_split" in jc.flags


def test_mark_tandem_finds_an_array_at_the_interval_end():
    seq = _array_contig()            # array at 30,000-70,160
    jc = idf.decide("asm", _junction(), _units(["contradicts"] * 3), idf.IdentifyParams(), [])
    b = jc.breaks[0]
    b.pos, b.lo, b.hi = 10_000, 5_000, 31_000
    idf.mark_tandem(jc, seq, rd.ReadParams())
    assert b.tandem is not None and b.tandem[0] <= 31_000 and "tandem_repeat" in jc.flags
    b.tandem, b.hi = None, 20_000
    jc.flags.remove("tandem_repeat")
    idf.mark_tandem(jc, seq, rd.ReadParams())
    assert b.tandem is None and "tandem_repeat" not in jc.flags


# ---------------------------------------------------------------------------
# Breaking
# ---------------------------------------------------------------------------
def test_plan_pieces():
    assert breaker.plan_pieces("c", 100, []) == [breaker.Piece("c", "c", 0, 100)]
    assert breaker.plan_pieces("c", 100, [60, 30, 30, 0, 100, 150]) == [
        breaker.Piece("c_b1", "c", 0, 30),
        breaker.Piece("c_b2", "c", 30, 60),
        breaker.Piece("c_b3", "c", 60, 100),
    ]


def test_plan_pieces_trims_gap_ns():
    seq = "ACGTA" + "N" * 10 + "CCGGT"
    assert breaker.plan_pieces("c", len(seq), [5], seq) == [
        breaker.Piece("c_b1", "c", 0, 5), breaker.Piece("c_b2", "c", 15, 20)]
    assert breaker.plan_pieces("c", len(seq), [5, 15], seq) == [
        breaker.Piece("c_b1", "c", 0, 5), breaker.Piece("c_b2", "c", 15, 20)]


def _write_tsv(path: Path, rows):
    with path.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=idf.COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in idf.COLUMNS})


def test_break_fasta_and_agp(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">a desc\nACGTACGTAC\nGTACGTACGT\n>b\nTTTTTGGGGG\n>c\nCCCCC\n")
    tsv = tmp_path / "m.tsv"
    _write_tsv(tsv, [
        {"assembly": "x", "contig": "a", "break_pos": "5", "action": "break"},
        {"assembly": "x", "contig": "a", "break_pos": "12", "action": "break"},
        {"assembly": "x", "contig": "b", "break_pos": "5", "action": "review"},
    ])
    breaks = breaker.read_breaks(tsv)
    assert breaks == {"a": [5, 12]}
    assert breaker.read_breaks(tsv, None) == {"a": [5, 12], "b": [5]}
    out = tmp_path / "out.fasta.gz"
    agp = tmp_path / "out.agp"
    pieces = breaker.break_fasta(fa, breaks, out, agp, wrap=4)
    assert [p.name for p in pieces] == ["a_b1", "a_b2", "a_b3"]
    with gzip.open(out, "rt") as fh:
        text = fh.read()
    recs = {}
    for block in text.split(">")[1:]:
        lines = block.strip().split("\n")
        recs[lines[0]] = "".join(lines[1:])
    assert recs == {"a_b1": "ACGTA", "a_b2": "CGTACGT", "a_b3": "ACGTACGT",
                    "b": "TTTTTGGGGG", "c": "CCCCC"}
    rows = [line.split("\t") for line in agp.read_text().splitlines() if not line.startswith("#")]
    assert rows[1] == ["a_b2", "1", "7", "1", "W", "a", "6", "12", "+"]
    assert rows[3] == ["b", "1", "10", "1", "W", "b", "1", "10", "+"]
    # each piece's sequence is the input at its AGP coordinates
    src = "ACGTACGTACGTACGTACGT"
    for r in rows[:3]:
        assert recs[r[0]] == src[int(r[6]) - 1:int(r[7])]


def test_break_rejects_unknown_contig_and_multi_assembly_tsv(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">a\nACGT\n")
    with pytest.raises(ValueError, match="not in"):
        breaker.break_fasta(fa, {"zz": [2]}, tmp_path / "o.fa", tmp_path / "o.agp")
    tsv = tmp_path / "m.tsv"
    _write_tsv(tsv, [{"assembly": "x", "contig": "a", "break_pos": "2", "action": "break"},
                     {"assembly": "y", "contig": "a", "break_pos": "2", "action": "break"}])
    with pytest.raises(ValueError, match="several assemblies"):
        breaker.read_breaks(tsv)
    assert breaker.read_breaks(tsv, assembly="y") == {"a": [2]}


# ---------------------------------------------------------------------------
# End to end on a synthetic results folder
# ---------------------------------------------------------------------------
def _rand(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(seq, rate, rng):
    s = list(seq)
    for i in range(0, len(s), int(1 / rate)):
        s[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[i]]
    return "".join(s)


def _write_assembly(root: Path, name: str, contigs, summary_rows, blocks=()):
    d = root / name
    d.mkdir(parents=True)
    with gzip.open(d / f"{name}.chrs.fasta.gz", "wt") as fh:
        for c, s in contigs.items():
            fh.write(f">{c}\n{s}\n")
    cols = ["contig", "original_name", "classification", "assigned_subgenome", "query_subgenome",
            "length", "reversed"]
    with (d / f"{name}.contig_summary.tsv").open("w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in summary_rows:
            w.writerow([r.get(k, "") for k in cols])
    write_macro_blocks(d / f"{name}.macro_blocks.tsv", blocks)


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not installed")
def test_identify_end_to_end(tmp_path):
    rng = random.Random(7)
    chr1 = _rand(400_000, rng)
    chr2 = _rand(400_000, rng)
    # target contig: chr1[0:200k] + chr2[200k:400k], a join no sister has
    target = chr1[:200_000] + chr2[200_000:]
    root = tmp_path / "out"
    blocks = [
        blk("ptgX", "chr1A", "+", 0, 190_000, 0, 190_000, length=400_000),
        blk("ptgX", "chr2A", "+", 210_000, 400_000, 210_000, 400_000, length=400_000),
    ]
    _write_assembly(root, "tgt", {"chr1A_f1": target},
                    [{"contig": "chr1A_f1", "original_name": "ptgX", "classification": "chrom_assigned",
                      "assigned_subgenome": "A", "length": 400_000, "reversed": "no"}], blocks)
    for k, name in enumerate(["sisA", "sisB"]):
        s1, s2 = _mutate(chr1, 0.01 + k * 0.002, rng), _mutate(chr2, 0.01 + k * 0.002, rng)
        _write_assembly(root, name, {"chr1A": s1, "chr2A": s2}, [
            {"contig": "chr1A", "original_name": "a", "classification": "chrom_assigned",
             "assigned_subgenome": "A", "length": len(s1)},
            {"contig": "chr2A", "original_name": "b", "classification": "chrom_assigned",
             "assigned_subgenome": "A", "length": len(s2)},
        ])
    # a sister built from the target's own contig carries the join, at 100% identity
    _write_assembly(root, "derived", {"chr1A": target}, [
        {"contig": "chr1A", "original_name": "ptgX", "classification": "chrom_assigned",
         "assigned_subgenome": "A", "length": len(target)}])

    from dnadis.misjoin.cli import main
    assert main(["identify", str(root), "-a", "tgt", "-t", "2",
                 "--min-extend", "20000", "--flank", "150000"]) == 0
    rows = list(csv.DictReader((root / "tgt" / "tgt.misjoins.tsv").open(), delimiter="\t"))
    assert len(rows) == 1
    r = rows[0]
    assert (r["call"], r["action"]) == ("misjoin", "break")
    assert r["contradictors"] == "sisA,sisB"
    assert r["same_genotype_carriers"] == "derived"
    assert abs(int(r["break_pos"]) - 200_000) <= 1_000

    # break the input with that row
    fa = tmp_path / "in.fa"
    fa.write_text(f">ptgX\n{target}\n")
    assert main(["break", str(fa), str(root / "tgt" / "tgt.misjoins.tsv"),
                 "-o", str(tmp_path / "broken.fa")]) == 0
    agp = (tmp_path / "broken.agp").read_text().splitlines()
    assert [ln.split("\t")[0] for ln in agp if not ln.startswith("#")] == ["ptgX_b1", "ptgX_b2"]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not installed")
def test_partner_segments_find_the_swapped_repeat(tmp_path):
    rng = random.Random(3)
    rep = _rand(8_000, rng)
    a1, b1, a2, b2 = (_rand(60_000, rng) for _ in range(4))
    # contig x: chrA - R - chrB; contig y: chrB' - R - chrA' (the reciprocal product)
    x = a1 + rep + b1
    y = b2 + rep + a2
    seqs = {("asm", "x"): x, ("asm", "y"): y}
    probes = {
        "px": sis.Region("px", "asm", "x", 30_000, 98_000),
        "py": sis.Region("py", "asm", "y", 30_000, 98_000),
    }
    got = idf._partner_segments(
        junctions={}, partners={"asm": {"x:j1": "y:j1", "y:j1": "x:j1"}},
        probe_of={("asm", "x:j1"): "px", ("asm", "y:j1"): "py"},
        interval_of={("asm", "x:j1"): (59_000, 61_000), ("asm", "y:j1"): (67_000, 69_000)},
        probes=probes, seqs=seqs, work_dir=tmp_path, params=idf.IdentifyParams(),
    )
    assert got == {("asm", "x:j1"): (60_000, 68_000), ("asm", "y:j1"): (60_000, 68_000)}
