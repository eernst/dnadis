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
