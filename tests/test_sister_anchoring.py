"""Tests for sister-assembly anchoring (multi-assembly second assignment pass)."""
from __future__ import annotations

import gzip
import os
import time
from pathlib import Path

from dnadis.alignment.sister_anchoring import (
    SisterCall,
    SisterVote,
    apply_sister_calls,
    build_call,
    decide,
    load_anchor_segments,
    load_ref_candidates,
    parse_sister_votes,
    segment_label,
    select_anchors,
    select_candidates,
    write_inputs,
)

MACRO_HEADER = [
    "contig", "contig_len", "ref_id", "chrom_id", "subgenome", "strand", "chain_id",
    "qstart", "qend", "qspan_bp", "union_bp", "matches", "aln_len", "identity",
    "score", "n_segments", "gene_count_chain", "ref_start", "ref_end",
]


def _write_macro_blocks(path: Path, blocks):
    """blocks: (contig, ref_id, strand, qstart, qend)"""
    with path.open("w") as out:
        out.write("\t".join(MACRO_HEADER) + "\n")
        for c, ref, strand, s, e in blocks:
            row = [c, "30000000", ref, ref[:-1], ref[-1], strand, "1", str(s), str(e),
                   str(e - s), str(e - s), "0", "0", "0.9", "0", "1", "0", "0", "0"]
            out.write("\t".join(row) + "\n")


def _paf_line(q, qlen, qs, qe, strand, t, ts, te, blk, de, primary=True):
    return "\t".join([
        q, str(qlen), str(qs), str(qe), strand, t, "30000000", str(ts), str(te),
        str(blk), str(blk), "60", f"tp:A:{'P' if primary else 'S'}", f"de:f:{de}",
    ]) + "\n"


def _row(name, cls, length, aligned, ref="NA"):
    return {"original_name": name, "classification": cls, "length": str(length),
            "total_aligned_bp": str(aligned), "assigned_ref_id": ref}


# --- candidate / anchor selection ---------------------------------------------

def test_select_candidates_and_anchors():
    rows = [
        _row("weak_frag", "chrom_fragment", 500_000, 50_000, "chr9T"),
        _row("short", "unclassified", 50_000, 0),
        _row("cobiont", "cobiont", 900_000, 0),
        _row("well_anchored", "chrom_assigned", 20_000_000, 8_000_000, "chr13A"),
        _row("weak_chrom", "chrom_assigned", 5_000_000, 500_000, "chr10A"),
        _row("good_frag", "chrom_fragment", 400_000, 200_000, "chr2A"),
    ]
    cands = {r["original_name"] for r in select_candidates(rows, 100_000, 0.25)}
    assert cands == {"weak_frag", "weak_chrom"}
    assert set(select_anchors(rows, 0.25)) == {"well_anchored", "good_frag"}


# --- anchor segment labels ----------------------------------------------------

def test_anchor_segments_keep_major_references(tmp_path):
    mb = tmp_path / "x.macro_blocks.tsv"
    _write_macro_blocks(mb, [
        # translocation carrier: 14A then 4A, plus a small spurious 6T block
        ("der", "chr14A", "+", 0, 8_000_000),
        ("der", "chr6T", "+", 8_100_000, 8_150_000),
        ("der", "chr4A", "-", 8_400_000, 15_000_000),
        ("other", "chr1A", "+", 0, 1_000_000),
    ])
    segs = load_anchor_segments(mb, {"der"})
    assert set(segs) == {"der"}
    assert {b[2] for b in segs["der"]} == {"chr14A", "chr4A"}
    # Overlap decides; the spurious block is ignored
    assert segment_label(segs["der"], 8_110_000, 8_140_000) == ("chr14A", "+")
    assert segment_label(segs["der"], 10_000_000, 10_100_000) == ("chr4A", "-")
    # No overlap: nearest major segment
    assert segment_label(segs["der"], 15_500_000, 15_600_000) == ("chr4A", "-")
    assert segment_label([], 0, 1) is None


# --- per-sister votes ---------------------------------------------------------

def test_parse_sister_votes_filters_and_labels(tmp_path):
    segs = {
        "anc13": [(0, 9_000_000, "chr13A", "+")],
        "anc2": [(0, 9_000_000, "chr2A", "-")],
    }
    paf = tmp_path / "sis.paf.gz"
    with gzip.open(paf, "wt") as fh:
        # placed: 600 kb of a 1 Mb contig on chr13A, gap-compressed identity 0.998
        fh.write(_paf_line("asmA|ctg1", 1_000_000, 0, 600_000, "+", "anc13", 0, 600_000, 700_000, 0.002))
        # low identity: ignored
        fh.write(_paf_line("asmA|ctg1", 1_000_000, 600_000, 900_000, "+", "anc2", 0, 300_000, 300_000, 0.05))
        # short block: ignored
        fh.write(_paf_line("asmA|ctg1", 1_000_000, 900_000, 905_000, "+", "anc2", 0, 5_000, 5_000, 0.0))
        # secondary: ignored
        fh.write(_paf_line("asmA|ctg1", 1_000_000, 0, 600_000, "+", "anc2", 0, 600_000, 600_000, 0.0, primary=False))
        # chimeric: two chromosomes each >= min_cov
        fh.write(_paf_line("asmA|ctg2", 1_000_000, 0, 400_000, "+", "anc13", 0, 400_000, 400_000, 0.001))
        fh.write(_paf_line("asmA|ctg2", 1_000_000, 500_000, 900_000, "-", "anc2", 0, 400_000, 400_000, 0.001))
        # the sister's own contig: ignored
        fh.write(_paf_line("sis|ctg9", 1_000_000, 0, 600_000, "+", "anc13", 0, 600_000, 600_000, 0.0))
        # reverse strand to a reverse-oriented anchor -> forward relative to reference
        fh.write(_paf_line("asmA|ctg3", 500_000, 0, 400_000, "-", "anc2", 0, 400_000, 400_000, 0.001))
    v = parse_sister_votes(paf, "sis", segs, min_aln=10_000, min_identity=0.98, min_cov=0.2)
    assert set(v) == {"asmA|ctg1", "asmA|ctg2", "asmA|ctg3"}
    assert v["asmA|ctg1"].label == "chr13A"
    assert abs(v["asmA|ctg1"].coverage - 0.6) < 1e-9
    assert v["asmA|ctg1"].reversed is False
    assert v["asmA|ctg2"].label == "CHIMERIC"
    assert v["asmA|ctg3"].label == "chr2A" and v["asmA|ctg3"].reversed is False


# --- decisions ----------------------------------------------------------------

def _votes(*labels_ids):
    return {f"s{i}": SisterVote(lab, 0.8, ident, False) for i, (lab, ident) in enumerate(labels_ids)}


def test_decide_outcomes():
    assert decide({}, "", set(), 0.75)[0] == "no_vote"
    # No reference evidence: unrestricted
    d = decide(_votes(("chr13A", 0.99), ("chr13A", 0.995)), "", set(), 0.75)
    assert d[:3] == ("assigned", "chr13A", 2)
    # Reference supports chr13A among its candidates: reassignment allowed
    d = decide(_votes(("chr13A", 0.99), ("chr13A", 0.99)), "chr9T", {"chr9T", "chr13A"}, 0.75)
    assert d[0] == "reassigned"
    # Same chromosome as the first pass
    votes = _votes(("chr13A", 0.99), ("chr13A", 0.99))
    assert decide(votes, "chr13A", {"chr13A"}, 0.75)[0] == "confirmed"
    # Chimeric majority
    assert decide(_votes(("CHIMERIC", 0.99)), "", set(), 0.75)[0] == "no_call_chimeric"
    # Too little agreement
    d = decide(_votes(("chr1A", 0.99), ("chr2A", 0.98)), "", set(), 0.75)
    assert d[0] == "no_call_disagree"
    # Majority agrees but the most similar sister does not
    d = decide(_votes(("chr1A", 0.98), ("chr1A", 0.98), ("chr1A", 0.98), ("chr2A", 0.999)), "", set(), 0.75)
    assert d[0] == "no_call_closest"


def test_decide_requires_two_agreeing_sisters():
    # One sister cannot outvote the reference, nor place a contig alone.
    one = _votes(("chr13A", 0.999))
    assert decide(one, "chr9T", {"chr9T", "chr13A"}, 0.75)[0] == "no_call_few_votes"
    assert decide(one, "", set(), 0.75)[0] == "no_call_few_votes"
    assert decide(one, "", set(), 0.75, min_votes=1)[0] == "assigned"


def test_decide_restricts_to_reference_supported_chromosomes():
    # Reference evidence only for chr1A; sisters (e.g. sharing a translocation)
    # place the contig on chr2A: keep the reference result.
    votes = _votes(("chr2A", 0.99), ("chr2A", 0.99), ("chr2A", 0.99))
    assert decide(votes, "chr1A", {"chr1A"}, 0.75)[0] == "conflict_keep_ref"
    # A homeolog of a supported chromosome may change the subgenome.
    votes = _votes(("chr18A", 0.99), ("chr18A", 0.99))
    assert decide(votes, "", {"chr18P"}, 0.75)[:2] == ("assigned", "chr18A")


def test_build_call_summarizes_agreeing_sisters():
    row = _row("ctg", "chrom_fragment", 1_000_000, 50_000, "chr9T")
    votes = {
        "a": SisterVote("chr13A", 0.9, 0.999, True),
        "b": SisterVote("chr13A", 0.5, 0.990, True),
        "c": SisterVote("chr13A", 0.7, 0.995, None),
    }
    call = build_call("asm", row, votes, {"chr9T", "chr13A"}, 0.75)
    assert call.decision == "reassigned" and call.applies
    assert call.ref_id == "chr13A" and call.votes_str == "3/3"
    assert call.closest_sister == "a"
    assert abs(call.identity - 0.995) < 1e-9 and abs(call.coverage - 0.7) < 1e-9
    assert call.reversed is True


def test_apply_sister_calls_only_applies_changes():
    def call(decision, ref):
        return SisterCall("asm", "x", 1, "debris", "", set(), decision, ref)
    best_ref = {"a": "", "b": "chr9T", "c": "chr1A", "d": ""}
    calls = {
        "a": call("assigned", "chr13A"),
        "b": call("reassigned", "chr13A"),
        "c": call("conflict_keep_ref", "chr2A"),
        "d": call("no_vote", ""),
    }
    applied = apply_sister_calls(best_ref, calls, {"chr13A": "chr13A_norm"})
    assert applied == {"a", "b"}
    assert best_ref == {"a": "chr13A_norm", "b": "chr13A_norm", "c": "chr1A", "d": ""}


# --- inputs -------------------------------------------------------------------

def test_load_ref_candidates_threshold(tmp_path):
    ev = tmp_path / "x.evidence_summary.tsv"
    ev.write_text(
        "contig\tref_id\tunion_bp\n"
        "c1\tchr13A\t30000\n"
        "c1\tchr9T\t41000\n"
        "c1\tchr6T\t9000\n"
        "c1\tchrC\t50000\n"
    )
    assert load_ref_candidates(ev, 10_000) == {"c1": {"chr13A", "chr9T"}}


def test_write_inputs_prefixes_and_keeps_unchanged_files(tmp_path):
    fa_a = tmp_path / "a.fa"
    fa_a.write_text(">c1 desc\nACGT\n>c2\nGGGG\n")
    fa_b = tmp_path / "b.fa"
    fa_b.write_text(">x1\nTTTT\n")
    work = tmp_path / "work"
    cand, anchors = write_inputs(
        {"a": fa_a, "b": fa_b}, {"a": {"c1"}}, {"a": {"c2"}, "b": {"x1"}}, work,
    )
    assert cand.read_text() == ">a|c1\nACGT\n"
    assert anchors["a"].read_text() == ">c2\nGGGG\n"
    assert set(anchors) == {"a", "b"}
    # Unchanged inputs keep their mtime so cached alignments stay valid
    old = cand.stat().st_mtime_ns
    os.utime(cand, ns=(old - 10_000_000_000, old - 10_000_000_000))
    old = cand.stat().st_mtime_ns
    time.sleep(0.01)
    write_inputs({"a": fa_a, "b": fa_b}, {"a": {"c1"}}, {"a": {"c2"}, "b": {"x1"}}, work)
    assert cand.stat().st_mtime_ns == old
