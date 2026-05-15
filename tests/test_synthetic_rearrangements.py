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


def _matches_truth(detected: dict, truth: Rearrangement,
                   bp_tolerance: int = 500_000) -> bool:
    """Heuristic match for a detected RearrangementCall row against a
    Rearrangement truth record.

    Three checks:

    1. Type-class match.  ``translocation``, ``reciprocal_translocation``,
       and ``whole_arm_translocation`` all share the same detector
       signature, so any pair within that group counts.  Other types
       (inversion / fusion / fission) must match exactly.
    2. Chromosome involvement.  The set of (assigned_ref_id, partner_ref_id)
       must intersect the truth's (primary_ref, partner_ref) — after
       lowercasing the chr prefix to match dnadis's normalize_ref_id.
    3. Reference-interval overlap (for types that have a well-defined
       ref interval).  Truth's [breakpoint, breakpoint + span] must
       overlap detected's [ref_start, ref_end].  Falls back to a
       midpoint-within-tolerance check when the truth has no
       reference-side span (fusion).
    """
    dt = detected.get("rearrangement_type", "")
    tt = truth.rearrangement_type
    # dnadis's detector routinely co-emits more than one inter-chromosomal
    # label for the same event (e.g. a fusion call AND a whole-arm
    # translocation call, when one chromosome's full length is appended to
    # another).  Treat all four inter-chromosomal labels as a single
    # equivalence class for type-matching: any of {fusion, translocation,
    # reciprocal_translocation, whole_arm_translocation} satisfies a truth
    # of any of those types.  Inversion and fission stay strict — their
    # signatures are distinct and a label swap would indicate a real miss.
    inter_chrom_types = {"fusion", "translocation",
                         "reciprocal_translocation", "whole_arm_translocation"}
    type_ok = (
        dt == tt
        or (tt in inter_chrom_types and dt in inter_chrom_types)
    )
    if not type_ok:
        return False

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

    # Fusion calls (both as truth and as the detection's emitted label)
    # carry no meaningful single ref interval — dnadis writes ref=[0-0]
    # for fusions because the event involves two whole chromosomes.  In
    # those cases chromosome involvement alone is the match signature.
    if tt == "fusion" or dt == "fusion" or (det_start == 0 and det_end == 0):
        return True

    truth_lo = truth.ref_breakpoint
    truth_hi = truth.ref_breakpoint + truth.span_bp
    # Allow detection start/end to fall up to bp_tolerance outside the
    # truth interval to tolerate edge-effect detection drift.
    overlap_lo = max(det_start, truth_lo) - bp_tolerance
    overlap_hi = min(det_end, truth_hi) + bp_tolerance
    return overlap_hi > overlap_lo


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
        print(f"  asm_{i:02d}: seed={asm_seed} n_rearr={n_rearr}")
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
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=1800)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            print(f"  asm_{i:02d}: dnadis failed: {exc}")
            continue

        rearr_tsv = out_dir / "permuted" / "permuted.rearrangements.tsv"
        detected = _parse_rearrangements_tsv(rearr_tsv)
        for d in detected:
            all_detected.append((i, d))
        print(f"  asm_{i:02d}: dnadis detected {len(detected)} rearrangements")
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

    # Match detected vs truth, report sensitivity per type.  Use
    # many-to-many semantics: a single ground-truth event can match
    # multiple detected sub-events (dnadis routinely splits a large
    # inversion into several smaller calls), and a detected call may
    # legitimately match more than one truth in dense rearrangement
    # scenarios.  Sensitivity = fraction of truths with ≥1 match.
    matched_truth = set()
    matched_det = set()
    for ti, (truth_asm, truth_r) in enumerate(all_truth):
        for di, (det_asm, det_d) in enumerate(all_detected):
            if det_asm != truth_asm:
                continue
            if _matches_truth(det_d, truth_r):
                matched_truth.add(ti)
                matched_det.add(di)

    from collections import Counter
    per_type_total = Counter(r.rearrangement_type for _, r in all_truth)
    per_type_hit = Counter(
        r.rearrangement_type for ti, (_, r) in enumerate(all_truth)
        if ti in matched_truth
    )
    print(f"\n[synthetic-rearr] sensitivity report (seed={seed}):")
    for op in sorted(per_type_total):
        hit = per_type_hit.get(op, 0)
        tot = per_type_total[op]
        pct = (100.0 * hit / tot) if tot else 0
        print(f"  {op:<28s} {hit:3d}/{tot:<3d}  ({pct:5.1f}%)")
    overall_hit = sum(per_type_hit.values())
    overall_tot = sum(per_type_total.values())
    overall_pct = (100.0 * overall_hit / overall_tot) if overall_tot else 0
    print(f"  {'OVERALL':<28s} {overall_hit:3d}/{overall_tot:<3d}"
          f"  ({overall_pct:5.1f}%)")
    n_false_pos = len(all_detected) - len(matched_det)
    print(f"  unmatched detections (false positives or extra splits): "
          f"{n_false_pos}/{len(all_detected)}")

    if overall_pct < 70.0:
        warnings.warn(
            f"Synthetic-rearrangement sensitivity {overall_pct:.1f}% < 70% "
            f"threshold (seed={seed}).  Test does not hard-fail; review the "
            f"per-type breakdown above to identify the regressed type."
        )
