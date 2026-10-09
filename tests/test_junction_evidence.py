"""Tests for fusion-junction depth/alignment evidence."""
from __future__ import annotations

import gzip

from dnadis.analysis.junction_evidence import (
    FusionThresholds, annotate_fusion_evidence, assess_alignment_continuity,
    classify_fusion_evidence, find_fusion_junctions, score_depth_uniformity,
)
from dnadis.analysis.read_depth import parse_clipping_sam, read_depth_windows
from dnadis.models import FusionJunction, RearrangementCall
from dnadis.output.tsv_output import write_fusion_junctions_tsv, write_rearrangements_tsv

MB = 1_000_000


def _row(contig, clen, ref, qs, qe, rs=0, re_=0, strand="+"):
    return (contig, clen, ref, ref, "", strand, 1, qs, qe, qe - qs, qe - qs,
            qe - qs, qe - qs, 0.98, 100.0, 3, 0, rs, re_ or (qe - qs))


def _rows(gap=0):
    return [_row("ctg", 18 * MB, "chr1", 0, 12 * MB),
            _row("ctg", 18 * MB, "chr2", 12 * MB + gap, 18 * MB)]


def _junction(**kw):
    base = dict(contig="ctg", left_ref="chr1", right_ref="chr2",
                junction_start=12 * MB, junction_end=12 * MB,
                junction_point=12 * MB, contig_len=18 * MB,
                left_span=12 * MB, right_span=6 * MB)
    base.update(kw)
    return FusionJunction(**base)


def _windows(depths_by_start, win=1000):
    return lambda c, s, e: [(ws, ws + win, d) for ws, d in depths_by_start
                            if ws + win > s and ws < e]


def _flat(depth, jp=12 * MB, junction_depth=None, span=8000, win=1000):
    out = []
    for ws in range(jp - span, jp + span, win):
        in_j = jp - win <= ws < jp + win
        out.append((ws, junction_depth if (in_j and junction_depth is not None) else depth))
    return out


def test_junction_location_abutting_and_gapped():
    j = find_fusion_junctions(_rows(), {"ctg": 18 * MB})[0]
    assert (j.left_ref, j.right_ref) == ("chr1", "chr2")
    assert j.junction_point == 12 * MB and j.aln_gap_bp == 0
    j = find_fusion_junctions(_rows(gap=200_000), {"ctg": 18 * MB})[0]
    assert j.junction_start == 12 * MB and j.junction_end == 12 * MB + 200_000
    assert j.junction_point == 12 * MB + 100_000


def test_junction_overlap_and_minor_arm_ignored():
    rows = [_row("ctg", 18 * MB, "chr1", 0, 12 * MB),
            _row("ctg", 18 * MB, "chr2", 11 * MB, 18 * MB),
            _row("ctg", 18 * MB, "chr3", 5 * MB, 5 * MB + 100_000)]
    js = find_fusion_junctions(rows, {"ctg": 18 * MB})
    assert len(js) == 1 and js[0].aln_gap_bp < 0
    assert js[0].junction_point == 11 * MB + 500_000


def test_depth_uniform_drop_spike_insufficient():
    thr = FusionThresholds()
    j = _junction()
    score_depth_uniformity(j, _windows(_flat(30)), thr)
    assert j.depth_verdict == "uniform" and abs(j.depth_junction_ratio - 1) < 1e-9
    j = _junction()
    score_depth_uniformity(j, _windows(_flat(30, junction_depth=12)), thr)
    assert j.depth_verdict == "drop"
    j = _junction()
    score_depth_uniformity(j, _windows(_flat(30, junction_depth=70)), thr)
    assert j.depth_verdict == "spike"
    j = _junction()
    score_depth_uniformity(j, _windows(_flat(4)), thr)
    assert j.depth_verdict == "insufficient_depth"
    j = _junction()
    score_depth_uniformity(j, lambda c, s, e: [], thr)
    assert j.depth_verdict == "insufficient_depth"


def test_depth_asymmetric_flanks():
    w = [(ws, 20.0 if ws < 12 * MB else 40.0) for ws, _ in _flat(0)]
    j = _junction()
    score_depth_uniformity(j, _windows(w), FusionThresholds())
    assert j.depth_verdict == "asymmetric" and j.depth_flank_ratio == 0.5


def test_alignment_continuity():
    thr = FusionThresholds()
    for gap, exp in ((0, "abutting"), (150_000, "gapped"), (-5000, "overlapping")):
        j = _junction(aln_gap_bp=gap)
        assess_alignment_continuity(j, thr)
        assert j.aln_verdict == exp


def test_classification_matrix():
    j = _junction(depth_verdict="uniform", aln_verdict="abutting")
    classify_fusion_evidence(j)
    assert (j.fusion_call, j.fusion_confidence) == ("biological_fusion", "medium")
    assert j.evidence_codes == "Aln,Depth"

    j = _junction(depth_verdict="uniform", aln_verdict="abutting",
                  clip_spanning=20, clip_clipped=0)
    classify_fusion_evidence(j)
    assert (j.fusion_call, j.fusion_confidence) == ("biological_fusion", "high")

    j = _junction(depth_verdict="drop", aln_verdict="gapped",
                  clip_spanning=2, clip_clipped=18)
    classify_fusion_evidence(j)
    assert (j.fusion_call, j.fusion_confidence) == ("assembly_artifact", "high")
    assert set(j.evidence_codes.split(",")) <= {"Aln", "Depth", "Split"}

    j = _junction(depth_verdict="drop", aln_verdict="abutting")
    classify_fusion_evidence(j)
    assert (j.fusion_call, j.fusion_confidence) == ("assembly_artifact", "low")

    j = _junction(depth_verdict="insufficient_depth", aln_verdict="abutting")
    classify_fusion_evidence(j)
    assert j.fusion_call == "inconclusive" and j.evidence_codes == ""

    # overlapping arms: no defined junction, so depth cannot decide
    for dv in ("uniform", "drop"):
        j = _junction(depth_verdict=dv, aln_verdict="overlapping")
        classify_fusion_evidence(j)
        assert j.fusion_call == "inconclusive" and j.evidence_codes == ""


def test_read_depth_windows(tmp_path):
    bed = tmp_path / "d.regions.bed.gz"
    with gzip.open(bed, "wt") as fh:
        for c in ("a", "b"):
            for i in range(10):
                fh.write(f"{c}\t{i*1000}\t{(i+1)*1000}\t{i+1}.5\n")
    w = read_depth_windows(bed, "b", 2500, 5000)
    assert [x[0] for x in w] == [2000, 3000, 4000]
    assert read_depth_windows(bed, "zzz", 0, 100) == []
    assert read_depth_windows(tmp_path / "missing.gz", "a", 0, 10) == []
    assert read_depth_windows(None, "a", 0, 10) == []


def test_parse_clipping_sam():
    def rec(pos, cigar, extra=""):
        return f"r\t0\tctg\t{pos}\t60\t{cigar}\t*\t0\t0\tA\tA{extra}"
    lines = ["@HD\tVN:1.6",
             rec(1000, "5000M"),                       # spans 1000..6000
             rec(1000, "2000M3000S"),                  # clip at 3000
             rec(2000, "1000M", "\tSA:Z:x,1,+,1M,60,0;")]  # split read ending at 3000
    assert parse_clipping_sam(lines, 3000, 3000) == (1, 2)


def test_annotate_and_tsv_roundtrip(tmp_path):
    bed = tmp_path / "depth.regions.bed.gz"
    with gzip.open(bed, "wt") as fh:
        for ws in range(0, 18 * MB, 1000):
            if 12 * MB - 8000 <= ws < 12 * MB + 8000:
                fh.write(f"ctg\t{ws}\t{ws+1000}\t30\n")
    call = RearrangementCall("ctg", "ctg", "chr1", "fusion", "chr2", 0, 18 * MB,
                             0, 0, 18 * MB, "+", "medium", "e")
    other = RearrangementCall("ctg", "ctg", "chr1", "inversion", None, 0, 1, 0, 0,
                              1, "-", "low", "e")
    js = annotate_fusion_evidence([call, other], _rows(), {"ctg": 18 * MB}, bed,
                                  None, FusionThresholds())
    assert len(js) == 1 and call.fusion_call == "biological_fusion"
    assert call.confidence == "high" and call.junction_start == 12 * MB
    assert other.fusion_call == ""

    out = tmp_path / "r.tsv"
    write_rearrangements_tsv(out, [call, other])
    lines = out.read_text().splitlines()
    hdr = lines[0].split("\t")
    assert hdr[-1] == "evidence_codes" and "fusion_call" in hdr
    assert all(len(l.split("\t")) == len(hdr) for l in lines[1:])
    assert lines[2].split("\t")[hdr.index("fusion_call")] == ""

    out2 = tmp_path / "j.tsv"
    write_fusion_junctions_tsv(out2, js)
    rows = out2.read_text().splitlines()
    assert len(rows) == 2 and "biological_fusion" in rows[1]


def test_annotate_without_depth_is_alignment_only():
    call = RearrangementCall("ctg", "ctg", "chr1", "fusion", "chr2", 0, 1, 0, 0,
                             1, "+", "medium", "e")
    annotate_fusion_evidence([call], _rows(), {"ctg": 18 * MB}, None, None,
                             FusionThresholds())
    assert call.fusion_call == "inconclusive" and call.confidence == "medium"


def test_config_template_has_fusion_section():
    import argparse
    from dnadis.utils.config import dump_config_template, CONFIG_SCHEMA
    ns = argparse.Namespace(**{d: 1 for sec in CONFIG_SCHEMA.values() for d in sec.values()})
    out = dump_config_template(ns)
    assert "[fusion_evidence]" in out and "fusion_split_reads" in out


def test_only_junctions_with_a_fusion_call_are_scored():
    rows = _rows() + [_row("other", 1 * MB, "chr3", 0, 500_000),
                      _row("other", 1 * MB, "chr4", 500_000, 1 * MB)]
    call = RearrangementCall("ctg", "ctg", "chr1", "fusion", "chr2", 0, 1, 0, 0,
                             1, "+", "medium", "e")
    js = annotate_fusion_evidence([call], rows, {"ctg": 18 * MB, "other": MB},
                                  None, None, FusionThresholds())
    assert [j.contig for j in js] == ["ctg"]
