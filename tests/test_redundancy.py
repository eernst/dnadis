"""Tests for own-assembly redundancy evidence."""
from __future__ import annotations

import gzip

import pytest

from dnadis.detection.redundancy import best_hit, parse_redundancy_paf, region_depth_median


def _paf(q, qlen, qs, qe, t, ts, te, blen, de, tp="P"):
    return "\t".join(map(str, [
        q, qlen, qs, qe, "+", t, 10_000_000, ts, te, blen - 10, blen, 60,
        f"tp:A:{tp}", f"de:f:{de}",
    ])) + "\n"


def test_collinear_copy_vs_scattered_repeats(tmp_path):
    paf = tmp_path / "x.paf"
    lines = [
        # allelic copy: two blocks covering 95% of a 100 kb contig, adjacent on chrA
        _paf("hap", 100_000, 0, 60_000, "chrA", 1_000_000, 1_060_500, 60_500, 0.01),
        _paf("hap", 100_000, 60_000, 95_000, "chrA", 1_061_000, 1_096_000, 35_000, 0.02),
        # repeat-only: 8 kb hits scattered over 5 Mb of chrB
        *[_paf("rep", 100_000, i * 10_000, i * 10_000 + 8_000, "chrB",
               i * 600_000, i * 600_000 + 8_000, 8_000, 0.05) for i in range(8)],
        # secondary alignments are ignored
        _paf("hap", 100_000, 0, 100_000, "chrC", 0, 100_000, 100_000, 0.0, tp="S"),
    ]
    paf.write_text("".join(lines))
    hits = parse_redundancy_paf(paf)

    hap = best_hit(hits["hap"])
    assert hap.target == "chrA"
    assert hap.coverage == pytest.approx(0.95)
    assert hap.identity == pytest.approx(1 - (0.01 * 60_500 + 0.02 * 35_000) / 95_500)
    assert (hap.target_start, hap.target_end) == (1_000_000, 1_096_000)
    assert "chrC" not in hits["hap"]

    rep = best_hit(hits["rep"])
    assert rep.coverage < 0.2  # 0.64 in total, but not within one window
    assert rep.n_alignments <= 2


def test_region_depth_median(tmp_path):
    bed = tmp_path / "d.regions.bed.gz"
    with gzip.open(bed, "wt") as fh:
        for i, d in enumerate([10, 20, 30, 40, 50]):
            fh.write(f"chrA\t{i * 1000}\t{(i + 1) * 1000}\t{d}\n")
        fh.write("chrB\t0\t1000\t99\n")
    med = region_depth_median(bed, {"chrA": [(1000, 4000)]})
    assert med == {("chrA", 1000, 4000): 30}
