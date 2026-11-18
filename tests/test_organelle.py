"""Organelle (chrC/chrM) calls from BLAST evidence.

Each case is a small synthetic BLAST table: one row per HSP with query and
reference (subject) coordinates.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pytest

from dnadis.detection.organelle import (
    OrganelleThresholds,
    call_organelles,
    parse_organelle_blast,
)
from dnadis.utils.sequence_utils import is_hifiasm_circular

CHRC = 140_000
CHRM = 470_000
REF = {"chrC": CHRC, "chrM": CHRM}


def _hsp(contig, organelle, qs, qe, ss, se, pident=99.9):
    """One BLAST outfmt-6 row (1-based inclusive coordinates)."""
    length = qe - qs + 1
    return "\t".join(map(str, [
        contig, organelle, pident, length, 0, 0, qs, qe, ss, se, 0.0, 1000, 0,
    ]))


def _call(tmp_path: Path, rows, lengths, circular=None, **kw):
    blast = tmp_path / "organelle_blast.txt"
    blast.write_text("\n".join(rows) + "\n")
    return call_organelles(
        parse_organelle_blast(blast), lengths, REF,
        OrganelleThresholds(**kw), circular_contigs=circular,
    )


def _tiled(contig, organelle, contig_len, ref_start, ref_end, pident=99.9):
    """HSPs that align the whole contig to ref[ref_start:ref_end], repeating
    the reference window as often as needed (a concatemer when the contig is
    longer than the window)."""
    rows, q, unit = [], 1, ref_end - ref_start
    while q <= contig_len:
        n = min(unit, contig_len - q + 1)
        rows.append(_hsp(contig, organelle, q, q + n - 1,
                         ref_start + 1, ref_start + n, pident))
        q += n
    return rows


class TestCompleteCalls:
    def test_oversized_mitogenome_is_complete(self, tmp_path):
        """A 1.85x contig covering all of chrM (repeat-graph walk) is complete."""
        L = int(1.85 * CHRM)
        calls = _call(tmp_path, _tiled("ptg1l", "chrM", L, 0, CHRM), {"ptg1l": L})
        assert calls.chrM_contig == "ptg1l"
        hit = calls.hits["ptg1l"]
        assert hit.status == "complete" and hit.ref_coverage == pytest.approx(1.0)
        assert hit.length_ratio == pytest.approx(1.85, abs=1e-3)

    def test_mitogenome_above_cap_is_not_complete(self, tmp_path):
        L = int(3.2 * CHRM)
        calls = _call(tmp_path, _tiled("ptg1l", "chrM", L, 0, CHRM), {"ptg1l": L})
        assert calls.chrM_contig is None
        # Still the only contig covering chrM: a one-contig partial set
        assert calls.partial["chrM"] == ["ptg1l"]

    def test_collapsed_mito_repeats_shorter_than_reference(self, tmp_path):
        """0.66x of the reference length but covering all of it (the reference
        carries repeat copies the assembly collapses) is complete."""
        L = int(0.66 * CHRM)
        rows = [_hsp("ptg2c", "chrM", 1, L, 1, L),
                _hsp("ptg2c", "chrM", 1, CHRM - L, L + 1, CHRM)]
        calls = _call(tmp_path, rows, {"ptg2c": L})
        assert calls.chrM_contig == "ptg2c"

    def test_query_coverage_alone_is_not_enough(self, tmp_path):
        """Contig fully aligned, about reference length, but covering 68% of
        chrM: not complete (the old query-coverage test passed it)."""
        L = CHRM
        cov = int(0.68 * CHRM)
        rows = [_hsp("ptg3l", "chrM", 1, cov, 1, cov),
                _hsp("ptg3l", "chrM", cov + 1, L, 1, L - cov)]
        calls = _call(tmp_path, rows, {"ptg3l": L})
        assert calls.chrM_contig is None
        assert calls.hits["ptg3l"].ref_coverage == pytest.approx(0.68, abs=1e-3)
        assert "ptg3l" in calls.debris

    def test_chrC_in_length_window_but_partial_reference(self, tmp_path):
        """chrC-length contig covering 77.5% of chrC is not complete."""
        L = int(0.98 * CHRC)
        cov = int(0.775 * CHRC)
        rows = [_hsp("ptg4l", "chrC", 1, cov, 1, cov),
                _hsp("ptg4l", "chrC", cov + 1, L, 1, L - cov)]
        calls = _call(tmp_path, rows, {"ptg4l": L})
        assert calls.chrC_contig is None

    def test_chrC_missing_inverted_repeat_is_below_length_floor(self, tmp_path):
        """A contig lacking one inverted-repeat copy covers the whole reference
        (its IR aligns to both copies) but is 0.85x long: not complete."""
        L = int(0.85 * CHRC)
        rows = [_hsp("ptg5l", "chrC", 1, L, 1, L),
                _hsp("ptg5l", "chrC", L - (CHRC - L) + 1, L, L + 1, CHRC)]
        calls = _call(tmp_path, rows, {"ptg5l": L})
        assert calls.hits["ptg5l"].ref_coverage == pytest.approx(1.0)
        assert calls.chrC_contig is None


class TestConcatemers:
    def test_chrC_concatemer_above_cap_and_one_unit_preferred(self, tmp_path):
        lengths = {"cat7l": int(7.5 * CHRC), "cat14l": int(1.4 * CHRC), "unit1l": CHRC + 300}
        rows = []
        for c, L in lengths.items():
            rows += _tiled(c, "chrC", L, 0, CHRC)
        calls = _call(tmp_path, rows, lengths)
        assert calls.chrC_contig == "unit1l"
        # Candidates that lose the ranking stay organelle debris
        assert {"cat7l", "cat14l"} <= calls.debris
        assert calls.hits["cat14l"].status == "debris"

    def test_only_concatemer_becomes_partial(self, tmp_path):
        L = int(7.5 * CHRC)
        calls = _call(tmp_path, _tiled("cat7l", "chrC", L, 0, CHRC), {"cat7l": L})
        assert calls.chrC_contig is None
        assert calls.partial["chrC"] == ["cat7l"]

    def test_circular_preferred_over_closer_length(self, tmp_path):
        lengths = {"ptg10l": CHRC, "ptg11c": int(1.3 * CHRC)}
        rows = []
        for c, L in lengths.items():
            rows += _tiled(c, "chrC", L, 0, CHRC)
        calls = _call(tmp_path, rows, lengths)
        assert calls.chrC_contig == "ptg11c"
        assert calls.debris == {"ptg10l"}


class TestPartialSets:
    def test_split_mitogenome_forms_partial_set(self, tmp_path):
        """Three subgenomic contigs (55%, 50%, 30% of chrM, overlapping) and one
        redundant piece: the fewest that reach 90% form the set."""
        spans = {
            "ptg20c": (0, int(0.55 * CHRM)),
            "ptg21l": (int(0.45 * CHRM), int(0.95 * CHRM)),
            "ptg22l": (int(0.70 * CHRM), CHRM),
            "ptg23l": (int(0.10 * CHRM), int(0.30 * CHRM)),  # inside ptg20c
        }
        rows, lengths = [], {}
        for c, (s, e) in spans.items():
            lengths[c] = e - s
            rows.append(_hsp(c, "chrM", 1, e - s, s + 1, e))
        calls = _call(tmp_path, rows, lengths, circular={"ptg20c"})
        assert calls.chrM_contig is None
        assert calls.partial["chrM"] == ["ptg20c", "ptg21l"]
        assert calls.set_ref_coverage["chrM"] == pytest.approx(0.95, abs=1e-3)
        assert {"ptg22l", "ptg23l"} <= calls.debris
        assert calls.hits["ptg21l"].status == "partial"

    def test_set_below_threshold_is_debris(self, tmp_path):
        rows = [_hsp("a", "chrM", 1, 200_000, 1, 200_000),
                _hsp("b", "chrM", 1, 150_000, 200_001, 350_000)]
        calls = _call(tmp_path, rows, {"a": 200_000, "b": 150_000})
        assert "chrM" not in calls.partial
        assert calls.debris == {"a", "b"}
        assert calls.organelle_contigs() == {"a", "b"}

    def test_mostly_nuclear_contig_cannot_join_set(self, tmp_path):
        """A contig with only 60% organelle sequence (e.g. a nuclear insertion
        with flanks) adds reference coverage but is not a set member."""
        rows = [_hsp("m1", "chrM", 1, 300_000, 1, 300_000),
                _hsp("numt", "chrM", 50_001, 170_000, 300_001, 420_000)]
        lengths = {"m1": 300_000, "numt": 200_000}
        calls = _call(tmp_path, rows, lengths)
        assert "chrM" not in calls.partial
        assert calls.hits["numt"].coverage == pytest.approx(0.6)
        assert "numt" in calls.debris

    def test_set_threshold_flag(self, tmp_path):
        # Together 400 kb of 470 kb (0.85): below the 0.90 default
        rows = [_hsp("a", "chrM", 1, 250_000, 1, 250_000),
                _hsp("b", "chrM", 1, 150_000, 250_001, 400_000)]
        lengths = {"a": 250_000, "b": 150_000}
        assert "chrM" not in _call(tmp_path, rows, lengths).partial
        loose = _call(tmp_path, rows, lengths, set_min_ref_cov=0.85)
        assert loose.partial["chrM"] == ["a", "b"]


class TestOrganelleAssignment:
    def test_plastid_fragment_with_mtpt_hit_is_chrC(self, tmp_path):
        """A chloroplast fragment also matches plastid-derived stretches of the
        mitogenome (85% at 97.7%); it is chrC debris, not chrM."""
        L = 20_000
        rows = [_hsp("frag", "chrC", 1, L, 30_001, 30_000 + L, 99.7),
                _hsp("frag", "chrM", 1, int(0.85 * L), 100_001, 100_000 + int(0.85 * L), 97.7)]
        calls = _call(tmp_path, rows, {"frag": L})
        assert calls.hits["frag"].organelle_type == "chrC"
        assert calls.debris == {"frag"}

    def test_low_query_coverage_is_not_called(self, tmp_path):
        rows = [_hsp("nuc", "chrM", 1, 30_000, 1, 30_000)]
        calls = _call(tmp_path, rows, {"nuc": 100_000})
        assert calls.organelle_contigs() == set()
        assert calls.hits == {}

    def test_blast_coordinates_are_inclusive(self, tmp_path):
        rows = [_hsp("x", "chrC", 1, CHRC, 1, CHRC)]
        calls = _call(tmp_path, rows, {"x": CHRC})
        assert calls.hits["x"].coverage == 1.0
        assert calls.hits["x"].ref_coverage == 1.0


@pytest.mark.parametrize("name, circular", [
    ("ptg000010c", True), ("ptg000010l", False),
    ("h1tg000049c", True), ("h2tg000019c", True), ("h1tg000049l", False),
    ("utg000003c", True), ("ptg000010c_RagTag", True), ("ptg000010ctg", False),
    ("LPlike_H1_h2tg000039c", False),
])
def test_hifiasm_circular_names(name, circular):
    assert is_hifiasm_circular(name) is circular


def test_legacy_length_tolerance_config_warns(tmp_path, caplog):
    from dnadis.utils.config import load_config

    cfg = tmp_path / "c.toml"
    cfg.write_text("[thresholds_organelle]\nchrM_len_tolerance = 0.2\norganelle_min_ref_cov = 0.9\n")
    with caplog.at_level(logging.WARNING):
        flat = load_config(cfg)
    assert flat == {"organelle_min_ref_cov": 0.9}
    assert "chrM_len_tolerance" in caplog.text


def test_classifier_names_partial_set_and_routes_fasta(tmp_path):
    """Partial members keep chrM_1/chrM_2 names and go to organelles FASTA."""
    from dnadis.classification.classifier import assign_contig_names, classify_all_contigs
    from dnadis.models import ChainEvidenceResult, OrganelleHit
    from dnadis.output.fasta_output import _CATEGORY_BY_CLASS

    fa = tmp_path / "q.fa"
    seqs = {"m1": "ACGT" * 50, "m2": "ACGT" * 40, "c1": "ACGT" * 30, "d1": "ACGT" * 10}
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    lengths = {k: len(v) for k, v in seqs.items()}
    hits = {
        "m1": OrganelleHit("chrM", 1.0, 0.999, 0.6, False, 0.55, "partial"),
        "m2": OrganelleHit("chrM", 1.0, 0.999, 0.5, False, 0.45, "partial"),
        "c1": OrganelleHit("chrC", 1.0, 0.999, 1.0, True, 0.97, "complete"),
        "d1": OrganelleHit("chrC", 0.55, 0.99, 0.1, False, 0.05, "debris"),
    }
    ev = ChainEvidenceResult(
        qlens_from_paf={}, qr_union_bp={}, qr_matches={}, qr_alnlen={},
        qr_gene_count={}, contig_total={}, contig_refs={}, qr_score_topk={},
        qr_weight_all={}, qr_nchains_kept={}, best_ref={}, best_score={},
        best_bp={}, second_ref={}, second_score={}, second_bp={},
        chain_segments_rows=[], macro_block_rows=[], chain_summary_rows=[],
        qr_ref_span_bp={}, qr_best_chain_ident={},
    )
    clfs = classify_all_contigs(
        query_fasta=fa, query_lengths=lengths, best_ref={}, chr_like_minlen=10**9,
        ev=ev, ref_gene_counts={}, chrC_contig="c1", chrM_contig=None,
        organelle_debris={"d1"}, rdna_contigs=set(), cobionts={},
        chromosome_debris=set(), other_debris=set(), add_subgenome_suffix=None,
        organelle_partial={"chrM": ["m1", "m2"]}, organelle_hits=hits,
    )
    by = {c.original_name: c for c in clfs}
    assert (by["m1"].classification, by["m1"].new_name) == ("organelle_partial", "chrM_1")
    assert (by["m2"].new_name, by["m2"].assigned_ref_id) == ("chrM_2", "chrM")
    assert by["m1"].organelle_ref_cov == 0.55
    assert (by["c1"].classification, by["c1"].new_name) == ("organelle_complete", "chrC")
    assert by["c1"].classification_confidence == "medium"  # ref cov < 0.98
    assert by["d1"].classification == "organelle_debris"
    assert by["d1"].classification_confidence == "low"
    assign_contig_names(clfs, lengths, None, overwrite=True)
    assert {c.new_name for c in clfs if c.original_name in ("m1", "m2", "c1")} == {"chrM_1", "chrM_2", "chrC"}
    assert _CATEGORY_BY_CLASS["organelle_partial"] == "organelles"
