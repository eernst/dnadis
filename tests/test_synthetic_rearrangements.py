"""Tests for the synthetic rearrangement helper library plus an
end-to-end integration test that runs the dnadis pipeline against
permuted Arabidopsis genomes and reports detection sensitivity.

The unit tests use small in-memory genomes and verify each operation
mechanically.  The integration test (gated on ``DNADIS_TAIR10_FASTA``)
generates N permuted assemblies, runs dnadis on each, parses the
rearrangements.tsv, matches detected calls against ground truth, and
prints a sensitivity report.  Per the testing policy this last test
never hard-fails on a sensitivity miss — it logs and warns.
"""
from __future__ import annotations

import csv
import os
import random
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import List, Tuple

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tests.synthetic_rearrangements import (  # noqa: E402
    Rearrangement,
    apply_fission,
    apply_fusion,
    apply_inversion,
    apply_reciprocal_translocation,
    apply_whole_arm_translocation,
    generate_random_assembly,
    load_genome,
    revcomp,
    write_genome,
)


# ---------------------------------------------------------------------------
# Unit tests — operations on tiny synthetic genomes
# ---------------------------------------------------------------------------
def _toy_genome():
    return {
        "chr1": "AAAA" "CCCC" "GGGG" "TTTT",   # 16 bp
        "chr2": "ACGT" "ACGT" "ACGT" "ACGT",   # 16 bp
        "chr3": "AAAA" "TTTT" "AAAA" "TTTT",   # 16 bp
    }


def test_revcomp_basic():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAACCC") == "GGGTTT"
    # Case is preserved per-base but the sequence is reversed.
    assert revcomp("NnNn") == "nNnN"


def test_load_write_genome_roundtrip(tmp_path):
    g = _toy_genome()
    path = tmp_path / "toy.fa"
    write_genome(g, path)
    loaded = load_genome(path)
    assert loaded == g


def test_apply_inversion_reverses_region():
    g = _toy_genome()
    r = apply_inversion(g, "chr1", 4, 12)  # invert "CCCCGGGG"
    assert g["chr1"] == "AAAA" + revcomp("CCCCGGGG") + "TTTT"
    assert g["chr1"] == "AAAA" + "CCCCGGGG" + "TTTT"  # palindromic test region
    assert r.rearrangement_type == "inversion"
    assert r.primary_ref == "chr1"
    assert r.partner_ref is None
    assert r.ref_breakpoint == 4
    assert r.span_bp == 8


def test_apply_inversion_non_palindromic():
    g = {"chr1": "A" * 4 + "ACGT" + "T" * 4}
    apply_inversion(g, "chr1", 4, 8)
    assert g["chr1"][4:8] == revcomp("ACGT") == "ACGT"  # palindrome of ACGT
    # Use a clearly non-palindromic stretch
    g = {"chr1": "A" * 4 + "AAAC" + "T" * 4}
    apply_inversion(g, "chr1", 4, 8)
    assert g["chr1"][4:8] == revcomp("AAAC") == "GTTT"


def test_apply_inversion_out_of_range_raises():
    g = _toy_genome()
    with pytest.raises(ValueError):
        apply_inversion(g, "chr1", -1, 5)
    with pytest.raises(ValueError):
        apply_inversion(g, "chr1", 5, 100)
    with pytest.raises(ValueError):
        apply_inversion(g, "chr1", 5, 5)


def test_apply_reciprocal_translocation_swaps_arms():
    g = {
        "chr1": "A" * 10 + "X" * 10,   # 20 bp
        "chr2": "B" * 10 + "Y" * 10,   # 20 bp
    }
    rs = apply_reciprocal_translocation(g, "chr1", 10, "chr2", 10)
    # chr1 = left arm of chr1 + right arm of chr2 (the swap)
    assert g["chr1"] == "A" * 10 + "Y" * 10
    assert g["chr2"] == "B" * 10 + "X" * 10
    assert len(rs) == 1
    r = rs[0]
    assert r.rearrangement_type == "reciprocal_translocation"
    assert {r.primary_ref, r.partner_ref} == {"chr1", "chr2"}
    assert r.ref_breakpoint == 10
    assert r.partner_breakpoint == 10


def test_apply_whole_arm_translocation_moves_arm():
    g = {"chr1": "A" * 5 + "X" * 5, "chr2": "B" * 5}  # chr1=10, chr2=5
    r = apply_whole_arm_translocation(g, "chr1", 5, "chr2")
    assert g["chr1"] == "A" * 5
    assert g["chr2"] == "B" * 5 + "X" * 5  # acceptor gains donor's arm
    assert r.rearrangement_type == "whole_arm_translocation"
    assert r.primary_ref == "chr1"
    assert r.partner_ref == "chr2"
    assert r.ref_breakpoint == 5
    assert r.partner_breakpoint == 5  # length of chr2 BEFORE append
    assert r.span_bp == 5


def test_apply_fusion_concatenates_and_drops_originals():
    g = {"chr1": "A" * 5, "chr2": "B" * 5, "chr3": "C" * 5}
    r = apply_fusion(g, "chr1", "chr2")
    assert "chr1" not in g
    assert "chr2" not in g
    fused = "fus_chr1_chr2"
    assert g[fused] == "A" * 5 + "B" * 5
    assert g["chr3"] == "C" * 5
    assert r.rearrangement_type == "fusion"
    assert r.primary_ref == "chr1"
    assert r.partner_ref == "chr2"
    assert r.ref_breakpoint == 5  # junction within fused contig
    assert fused in r.query_contigs


def test_apply_fusion_preserves_chromosome_order():
    g = {"chr1": "A", "chr2": "B", "chr3": "C", "chr4": "D"}
    apply_fusion(g, "chr2", "chr3")
    # chr1 first, then the fused product where chr2 was, then chr4
    assert list(g.keys()) == ["chr1", "fus_chr2_chr3", "chr4"]


def test_apply_fission_splits_chromosome():
    g = {"chr1": "AAAAACCCCC", "chr2": "BBBBB"}
    r = apply_fission(g, "chr1", 5)
    assert "chr1" not in g
    assert g["fis_chr1_a"] == "AAAAA"
    assert g["fis_chr1_b"] == "CCCCC"
    assert g["chr2"] == "BBBBB"
    assert r.rearrangement_type == "fission"
    assert r.ref_breakpoint == 5
    assert r.span_bp == 10


def test_generate_random_assembly_writes_valid_fasta(tmp_path):
    # Build a tiny reference; bypass min_chrom_length / min_segment_length
    # by passing small values.
    ref = tmp_path / "ref.fa"
    write_genome({
        "chr1": "A" * 20_000,
        "chr2": "C" * 20_000,
        "chr3": "G" * 20_000,
        "chrM": "T" * 5_000,   # organelle, should be excluded
    }, ref)
    out = tmp_path / "asm.fa"
    truth = generate_random_assembly(
        ref, out,
        seed=42,
        n_rearrangements=3,
        min_chrom_length=10_000,
        min_segment_length=2_000,
    )
    # Reproducibility: same seed → same operation sequence
    out2 = tmp_path / "asm2.fa"
    truth2 = generate_random_assembly(
        ref, out2, seed=42, n_rearrangements=3,
        min_chrom_length=10_000, min_segment_length=2_000,
    )
    assert load_genome(out) == load_genome(out2)
    assert len(truth) == len(truth2)
    for a, b in zip(truth, truth2):
        assert a.rearrangement_type == b.rearrangement_type
        assert a.primary_ref == b.primary_ref
        assert a.ref_breakpoint == b.ref_breakpoint

    # Different seed → different output
    out3 = tmp_path / "asm3.fa"
    truth3 = generate_random_assembly(
        ref, out3, seed=123, n_rearrangements=3,
        min_chrom_length=10_000, min_segment_length=2_000,
    )
    # Either operation sequence or contig contents should differ.
    assert (load_genome(out) != load_genome(out3)
            or [(r.rearrangement_type, r.primary_ref) for r in truth]
               != [(r.rearrangement_type, r.primary_ref) for r in truth3])


def test_generate_random_assembly_skips_organelles(tmp_path):
    ref = tmp_path / "ref.fa"
    write_genome({
        "chr1": "A" * 20_000,
        "chr2": "C" * 20_000,
        "chrM": "T" * 30_000,  # long enough to pass min_chrom_length…
        "chrC": "G" * 30_000,  # …but flagged as organelle
    }, ref)
    out = tmp_path / "asm.fa"
    truth = generate_random_assembly(
        ref, out, seed=7, n_rearrangements=5,
        min_chrom_length=10_000, min_segment_length=2_000,
    )
    for r in truth:
        assert r.primary_ref not in ("chrM", "chrC")
        assert r.partner_ref not in ("chrM", "chrC")


# ---------------------------------------------------------------------------
# Integration test — synthetic assemblies → dnadis → validation
# ---------------------------------------------------------------------------
TAIR10_FASTA = (
    Path(os.environ["DNADIS_TAIR10_FASTA"])
    if os.environ.get("DNADIS_TAIR10_FASTA") else None
)
DNADIS_PY = Path(__file__).resolve().parents[1] / "dnadis.py"


def _parse_rearrangements_tsv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def _norm_ref(name: str) -> str:
    """Mirror dnadis's normalize_ref_id: lowercase 'chr' prefix.  Used
    so capitalized reference IDs (e.g. TAIR10's Chr1) compare equal to
    the lowercase form (chr1) that dnadis writes to rearrangements.tsv.
    """
    from dnadis.utils.reference_utils import normalize_ref_id
    return normalize_ref_id(name) if name else ""


# dnadis emits multiple labels for the same inter-chromosomal event
# (a whole-arm translocation that moves an entire arm to the end of
# another chromosome can be described as a fusion, a whole_arm
# translocation, or both).  The test reports two parallel sensitivity
# metrics:
#   * STRICT — exact rearrangement_type match; tracks how often dnadis
#     applies the same label the ground truth uses.  Improvements to
#     the detector's labelling show up here.
#   * BROAD  — any inter-chromosomal label satisfies any inter-chromosomal
#     truth; tracks whether the underlying event was detected at all.
#     Inversion and fission stay strict in BROAD too — their signatures
#     are distinct and a label swap would indicate a real miss.
INTER_CHROM_TYPES = frozenset({
    "fusion", "translocation",
    "reciprocal_translocation", "whole_arm_translocation",
})


def _chrom_and_interval_ok(detected: dict, truth: Rearrangement,
                           bp_tolerance: int) -> bool:
    """Shared check used by both strict and broad matchers: chromosomes
    overlap and (for non-fusion events with a meaningful ref interval)
    the intervals overlap within tolerance.
    """
    det_refs = {
        _norm_ref(detected.get("assigned_ref_id", "")),
        _norm_ref(detected.get("partner_ref_id", "") or ""),
    } - {""}
    truth_refs = {_norm_ref(truth.primary_ref)}
    if truth.partner_ref:
        truth_refs.add(_norm_ref(truth.partner_ref))
    if not (det_refs & truth_refs):
        return False
    try:
        det_start = int(detected.get("ref_start", 0) or 0)
        det_end = int(detected.get("ref_end", 0) or 0)
    except ValueError:
        return False
    # Fusion calls carry no meaningful ref interval; chromosome
    # involvement is the only signature.
    if (truth.rearrangement_type == "fusion"
            or detected.get("rearrangement_type") == "fusion"
            or (det_start == 0 and det_end == 0)):
        return True
    truth_lo = truth.ref_breakpoint
    truth_hi = truth.ref_breakpoint + truth.span_bp
    overlap_lo = max(det_start, truth_lo) - bp_tolerance
    overlap_hi = min(det_end, truth_hi) + bp_tolerance
    return overlap_hi > overlap_lo


def _matches_truth_strict(detected: dict, truth: Rearrangement,
                          bp_tolerance: int = 500_000) -> bool:
    """Detection's rearrangement_type matches the truth's exactly, plus
    chromosomes overlap and (for non-fusion) reference intervals overlap.
    """
    if detected.get("rearrangement_type") != truth.rearrangement_type:
        return False
    return _chrom_and_interval_ok(detected, truth, bp_tolerance)


def _matches_truth_broad(detected: dict, truth: Rearrangement,
                         bp_tolerance: int = 500_000) -> bool:
    """As _matches_truth_strict, but treat all inter-chromosomal labels
    as a single equivalence class.  Inversion and fission stay strict.
    """
    dt = detected.get("rearrangement_type", "")
    tt = truth.rearrangement_type
    if not (dt == tt or (tt in INTER_CHROM_TYPES and dt in INTER_CHROM_TYPES)):
        return False
    return _chrom_and_interval_ok(detected, truth, bp_tolerance)


# Back-compat shim for any external caller.
_matches_truth = _matches_truth_broad


@pytest.mark.integration
def test_synthetic_rearrangement_detection(tmp_path):
    """Generate N permuted Arabidopsis assemblies, run dnadis on each,
    compare detected rearrangements to ground truth.  Per-type sensitivity
    is logged and warned-on but never hard-fails (see TODO testing item)."""
    if TAIR10_FASTA is None or not TAIR10_FASTA.exists():
        pytest.skip("Set DNADIS_TAIR10_FASTA to run synthetic rearrangement tests")

    n_asms = int(os.environ.get("DNADIS_SYNTHETIC_N", "10"))
    seed = int(os.environ.get("DNADIS_SYNTHETIC_SEED", str(int(time.time()))))
    threads = int(os.environ.get("DNADIS_SYNTHETIC_THREADS",
                                  str(min(8, os.cpu_count() or 4))))
    print(f"\n[synthetic-rearr] seed={seed} n_asms={n_asms} threads={threads}")

    rng = random.Random(seed)
    all_truth: List[Tuple[int, Rearrangement]] = []
    all_detected: List[Tuple[int, dict]] = []

    for i in range(n_asms):
        asm_seed = rng.randint(0, 2 ** 30)
        asm_dir = tmp_path / f"asm_{i:02d}"
        asm_dir.mkdir()
        asm_fa = asm_dir / "permuted.fasta"
        n_rearr = rng.randint(1, 3)
        truth = generate_random_assembly(
            TAIR10_FASTA, asm_fa,
            seed=asm_seed, n_rearrangements=n_rearr,
        )
        for r in truth:
            all_truth.append((i, r))
        print(f"  asm_{i:02d}: seed={asm_seed} n_rearr={n_rearr}", flush=True)
        for r in truth:
            print(f"    {r.rearrangement_type} {r.primary_ref}"
                  f"{':' + r.partner_ref if r.partner_ref else ''}"
                  f" @ {r.ref_breakpoint} ({r.span_bp/1e6:.1f} Mb)")

        out_dir = asm_dir / "out"
        cmd = [
            sys.executable, str(DNADIS_PY),
            "-r", str(TAIR10_FASTA),
            "-q", str(asm_fa),
            "-o", str(out_dir),
            "--skip-organelles", "--skip-rdna",
            "--skip-telomeres", "--skip-rdna-consensus",
            "--skip-compleasm", "--skip-plot",
            "-t", str(threads),
        ]
        # Stream stdout/stderr to log files instead of capture_output=True.
        # subprocess.run with capture_output buffers both pipes in memory
        # and only reads them after the child exits — on verbose pipelines
        # the OS pipe buffer (~64 KB) fills up mid-run and the dnadis
        # child deadlocks waiting for the parent to drain the pipe.
        stdout_log = asm_dir / "dnadis.stdout.log"
        stderr_log = asm_dir / "dnadis.stderr.log"
        try:
            with stdout_log.open("w") as so, stderr_log.open("w") as se:
                subprocess.run(cmd, check=True, stdout=so, stderr=se,
                               timeout=1800)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            print(f"  asm_{i:02d}: dnadis failed: {exc} "
                  f"(see {stderr_log})", flush=True)
            continue

        rearr_tsv = out_dir / "permuted" / "permuted.rearrangements.tsv"
        detected = _parse_rearrangements_tsv(rearr_tsv)
        for d in detected:
            all_detected.append((i, d))
        print(f"  asm_{i:02d}: dnadis detected {len(detected)} rearrangements",
              flush=True)
        for d in detected:
            print(
                f"      {d.get('rearrangement_type','?'):<24s}"
                f" {d.get('assigned_ref_id','?')}"
                f"{':' + d['partner_ref_id'] if d.get('partner_ref_id') else ''}"
                f" ref=[{d.get('ref_start','?')}-{d.get('ref_end','?')}]"
                f" span={d.get('span_bp','?')}"
                f" {d.get('strand','')}"
                f" conf={d.get('confidence','?')}"
            )

    # Compute strict and broad matches separately.  Use many-to-many
    # semantics under both: a single ground-truth event can match
    # multiple sub-events (e.g. dnadis splits a fission into two pieces).
    from collections import Counter, defaultdict

    matched_truth_strict: set = set()
    matched_det_strict: set = set()
    matched_truth_broad: set = set()
    matched_det_broad: set = set()
    # truth -> set of distinct detection labels that matched it broadly
    truth_emitted_labels: dict = defaultdict(set)

    for ti, (truth_asm, truth_r) in enumerate(all_truth):
        for di, (det_asm, det_d) in enumerate(all_detected):
            if det_asm != truth_asm:
                continue
            if _matches_truth_strict(det_d, truth_r):
                matched_truth_strict.add(ti)
                matched_det_strict.add(di)
            if _matches_truth_broad(det_d, truth_r):
                matched_truth_broad.add(ti)
                matched_det_broad.add(di)
                truth_emitted_labels[ti].add(
                    det_d.get("rearrangement_type", "?"))

    per_type_total = Counter(r.rearrangement_type for _, r in all_truth)
    per_type_strict_hit = Counter(
        r.rearrangement_type for ti, (_, r) in enumerate(all_truth)
        if ti in matched_truth_strict
    )
    per_type_broad_hit = Counter(
        r.rearrangement_type for ti, (_, r) in enumerate(all_truth)
        if ti in matched_truth_broad
    )

    def _fmt_row(label, hit, tot):
        pct = (100.0 * hit / tot) if tot else 0
        return f"  {label:<28s} {hit:3d}/{tot:<3d}  ({pct:5.1f}%)"

    print(f"\n[synthetic-rearr] sensitivity report (seed={seed})")
    print("\n  --- BROAD (any inter-chromosomal label satisfies any "
          "inter-chromosomal truth) ---")
    print("  Tracks whether the underlying event was detected at all.")
    for op in sorted(per_type_total):
        print(_fmt_row(op, per_type_broad_hit.get(op, 0), per_type_total[op]))
    overall_broad = sum(per_type_broad_hit.values())
    overall_tot = sum(per_type_total.values())
    print(_fmt_row("OVERALL", overall_broad, overall_tot))

    print("\n  --- STRICT (rearrangement_type label must match exactly) ---")
    print("  Tracks whether the detector emits the correct label.  Lower")
    print("  than BROAD means dnadis recognised the event but labelled it")
    print("  with a different (compatible) inter-chromosomal type.")
    for op in sorted(per_type_total):
        print(_fmt_row(op, per_type_strict_hit.get(op, 0), per_type_total[op]))
    overall_strict = sum(per_type_strict_hit.values())
    print(_fmt_row("OVERALL", overall_strict, overall_tot))

    # Label-confusion view: for each truth type, what labels did the
    # detector actually emit on matching detections?
    truth_type_to_label_counts: dict = defaultdict(Counter)
    for ti, (_, r) in enumerate(all_truth):
        for lbl in truth_emitted_labels.get(ti, ()):
            truth_type_to_label_counts[r.rearrangement_type][lbl] += 1
    print("\n  --- Label confusion (truth_type → emitted labels on matching detections) ---")
    for op in sorted(per_type_total):
        labels = truth_type_to_label_counts.get(op, Counter())
        if not labels:
            print(f"  {op:<28s} (no matching detections)")
            continue
        label_str = ", ".join(f"{k}={v}" for k, v in labels.most_common())
        print(f"  {op:<28s} {label_str}")

    # Extras: detections that don't match any truth, even broadly.
    n_extras = len(all_detected) - len(matched_det_broad)
    print(f"\n  --- Detection stats ---")
    print(f"  total detections: {len(all_detected)}")
    print(f"  detections matching ≥1 truth (broad): "
          f"{len(matched_det_broad)}")
    print(f"  detections matching no truth (extras): {n_extras}")
    if n_extras:
        print("    (extras may indicate real false positives or detector")
        print("     behaviour on chromosome boundaries.  Listed below:)")
        for di, (asm, d) in enumerate(all_detected):
            if di in matched_det_broad:
                continue
            print(f"      asm_{asm:02d}  {d.get('rearrangement_type','?')}"
                  f" {d.get('assigned_ref_id','?')}"
                  f"{':' + d['partner_ref_id'] if d.get('partner_ref_id') else ''}"
                  f" ref=[{d.get('ref_start','?')}-{d.get('ref_end','?')}]")

    broad_pct = (100.0 * overall_broad / overall_tot) if overall_tot else 0
    strict_pct = (100.0 * overall_strict / overall_tot) if overall_tot else 0
    if broad_pct < 70.0:
        warnings.warn(
            f"Synthetic-rearrangement BROAD sensitivity {broad_pct:.1f}% < 70% "
            f"(seed={seed}).  Detection coverage gap; investigate the per-type "
            f"breakdown above."
        )
    if strict_pct < 70.0:
        warnings.warn(
            f"Synthetic-rearrangement STRICT sensitivity {strict_pct:.1f}% < 70% "
            f"(seed={seed}).  Detection works but type labels frequently "
            f"disagree with truth.  Review the label-confusion table above."
        )
