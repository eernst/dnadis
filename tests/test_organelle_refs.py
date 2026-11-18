"""Organelle references from reference FASTAs: names, species tags, calls."""
from __future__ import annotations

import logging
import random

import pytest

from dnadis.detection.organelle import (
    call_organelles,
    parse_organelle_blast,
    prepare_organelle_references,
)
from dnadis.utils.reference_utils import (
    find_organelle_sequences,
    organelle_from_description,
    parse_organelle_id,
    ref_seq_kind,
)
from dnadis.utils.sequence_utils import read_fasta_sequences

TAGS = {"A", "P", "T"}


@pytest.mark.parametrize("name, tags, bare, expected", [
    ("chrC", TAGS, False, ("chrC", None)),
    ("ChrM", TAGS, False, ("chrM", None)),
    ("Pt", set(), False, ("chrC", None)),
    ("Mt", set(), False, ("chrM", None)),
    ("chrMt", set(), False, ("chrM", None)),
    ("chloroplast", set(), False, ("chrC", None)),
    ("mitochondrion", set(), False, ("chrM", None)),
    ("chrC_A", TAGS, False, ("chrC", "A")),
    ("chrCP", TAGS, False, ("chrC", "P")),
    ("Pt.T", TAGS, False, ("chrC", "T")),
    ("Mt-A", TAGS, False, ("chrM", "A")),
    ("chrC_X", TAGS, False, ("chrC", None)),   # suffix that is not a nuclear tag
    ("C", set(), True, ("chrC", None)),
    ("M", set(), True, ("chrM", None)),
    ("C", set(), False, None),
    ("chr1A", TAGS, True, None),
    ("scaffold_12", TAGS, True, None),
    ("CpG", TAGS, False, None),
])
def test_parse_organelle_id(name, tags, bare, expected):
    assert parse_organelle_id(name, tags, bare) == expected


def test_organelle_from_description():
    assert organelle_from_description("Sorghum bicolor chloroplast, complete genome") == "chrC"
    assert organelle_from_description("Zea mays mitochondrion, complete genome") == "chrM"
    assert organelle_from_description("Zea mays plastid complete sequence") == "chrC"
    # Not a complete organelle genome, or both organelles named
    assert organelle_from_description("chloroplast-derived insertion on chr5") is None
    assert organelle_from_description("chloroplast and mitochondrion, complete genome") is None


def test_tagged_organelle_kind():
    assert ref_seq_kind("chrC_P") == "organelle"
    assert ref_seq_kind("PtA") == "organelle"
    assert ref_seq_kind("chr1A") == "nuclear"


def test_find_bare_names_only_with_numbered_or_roman_chromosomes():
    numbered = [("1", "", 10**7), ("2", "", 10**7), ("C", "", 140_000), ("M", "", 400_000)]
    assert find_organelle_sequences(numbered) == {("chrC", None): "C", ("chrM", None): "M"}
    roman = [("I", "", 10**7), ("II", "", 10**7), ("C", "", 140_000)]
    assert find_organelle_sequences(roman) == {("chrC", None): "C"}
    named = [("chr1", "", 10**7), ("C", "", 140_000)]
    assert find_organelle_sequences(named) == {}


def test_find_multi_species_and_untagged(caplog):
    headers = [("chr1A", "", 10**7), ("chr1P", "", 10**7), ("chr1T", "", 10**7),
               ("chrC_A", "", 164_000), ("chrCP", "", 165_000), ("chrM", "", 262_000)]
    assert find_organelle_sequences(headers) == {
        ("chrC", "A"): "chrC_A", ("chrC", "P"): "chrCP", ("chrM", None): "chrM"}


def test_find_clash_uses_longest(caplog):
    headers = [("chr1", "", 10**7), ("chrC", "", 140_000),
               ("NC_000001.1", "Zea mays chloroplast, complete genome", 140_384)]
    with caplog.at_level(logging.WARNING):
        found = find_organelle_sequences(headers)
    assert found == {("chrC", None): "NC_000001.1"}
    assert "using the longest" in caplog.text


def _rand(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_prepare_organelle_references(tmp_path):
    rng = random.Random(0)
    ref = tmp_path / "ref.fa"
    seqs = {"chr1A": _rand(5000, rng), "chr1P": _rand(5000, rng),
            "chrC_A": _rand(1640, rng), "PtP": _rand(1650, rng), "Mt": _rand(2620, rng)}
    ref.write_text("".join(f">{k} desc\n{v}\n" for k, v in seqs.items()))
    c, m = prepare_organelle_references(ref, None, None, tmp_path / "w")
    assert read_fasta_sequences(c) == {"chrC_A": seqs["chrC_A"], "chrC_P": seqs["PtP"]}
    assert read_fasta_sequences(m) == {"chrM": seqs["Mt"]}

    # A user file replaces that organelle; tagged names keep their tag
    user = tmp_path / "user_cp.fa"
    u = {"NC_1": _rand(1600, rng), "chrC_T": _rand(1700, rng)}
    user.write_text("".join(f">{k}\n{v}\n" for k, v in u.items()))
    c2, m2 = prepare_organelle_references(ref, str(user), None, tmp_path / "w2")
    # T is not a nuclear tag here: both untagged, the longer one is used
    assert read_fasta_sequences(c2) == {"chrC": u["chrC_T"]}
    assert read_fasta_sequences(m2) == {"chrM": seqs["Mt"]}


def test_prepare_single_untagged_organelle_in_multi_species_reference(tmp_path):
    rng = random.Random(1)
    ref = tmp_path / "ref.fa"
    seqs = {"chr1A": _rand(5000, rng), "chr1T": _rand(5000, rng), "chrC": _rand(1640, rng)}
    ref.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    c, m = prepare_organelle_references(ref, None, None, tmp_path / "w")
    assert read_fasta_sequences(c) == {"chrC": seqs["chrC"]} and m is None


def _hsp(contig, ref, qs, qe, ss, se, pident=99.9):
    return "\t".join(map(str, [contig, ref, pident, qe - qs + 1, 0, 0, qs, qe, ss, se, 0.0, 1000, 0]))


def test_calls_label_species_and_margin(tmp_path):
    """Each contig takes its best species reference; the call records it and
    the identity margin over the other species."""
    L = 140_000
    rows = [_hsp("own", "chrC_A", 1, L, 1, L, 99.3), _hsp("own", "chrC_P", 1, L, 1, L, 98.8),
            _hsp("own", "chrC_T", 1, L, 1, L, 96.6)]
    blast = tmp_path / "b.txt"
    blast.write_text("\n".join(rows) + "\n")
    calls = call_organelles(parse_organelle_blast(blast), {"own": L},
                            {"chrC_A": L, "chrC_P": L + 1000, "chrC_T": L + 2000})
    hit = calls.hits["own"]
    assert calls.chrC_contig == "own"
    assert hit.ref_name == "chrC_A" and hit.organelle_type == "chrC"
    assert hit.ref_margin == pytest.approx(0.005, abs=1e-6)
    assert calls.ref_name == {"chrC": "chrC_A"}


def test_single_reference_has_no_margin(tmp_path):
    L = 140_000
    blast = tmp_path / "b.txt"
    blast.write_text(_hsp("x", "chrC", 1, L, 1, L) + "\n")
    calls = call_organelles(parse_organelle_blast(blast), {"x": L}, {"chrC": L})
    assert calls.hits["x"].ref_name == "chrC" and calls.hits["x"].ref_margin is None


def test_partial_set_within_one_species_reference(tmp_path):
    """Set members must share a reference: coordinates differ between species."""
    R = 400_000
    rows = [_hsp("a", "chrM_A", 1, 200_000, 1, 200_000),
            _hsp("b", "chrM_P", 1, 200_000, 200_001, 400_000),   # other species
            _hsp("c", "chrM_A", 1, 190_000, 200_001, 390_000)]
    blast = tmp_path / "b.txt"
    blast.write_text("\n".join(rows) + "\n")
    lengths = {"a": 200_000, "b": 200_000, "c": 190_000}
    calls = call_organelles(parse_organelle_blast(blast), lengths, {"chrM_A": R, "chrM_P": R})
    assert calls.partial["chrM"] == ["a", "c"]
    assert calls.ref_name["chrM"] == "chrM_A"
    assert "b" in calls.debris


def test_species_chosen_by_identity_not_aligned_fraction(tmp_path):
    """Slightly more of the contig aligns to another species' reference (e.g.
    an extra repeat hit), but the identity is higher to its own."""
    L = 140_000
    rows = [_hsp("own", "chrC_A", 1, 133_000, 1, 133_000, 98.9),
            _hsp("own", "chrC_T", 1, 139_000, 1, 139_000, 96.7)]
    blast = tmp_path / "b.txt"
    blast.write_text("\n".join(rows) + "\n")
    calls = call_organelles(parse_organelle_blast(blast), {"own": L},
                            {"chrC_A": L, "chrC_T": L})
    hit = calls.hits["own"]
    assert hit.ref_name == "chrC_A"
    assert hit.ref_margin == pytest.approx(0.022, abs=1e-6)
