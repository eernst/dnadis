"""Tests for placing chromosome fragments into input scaffold gaps."""
from __future__ import annotations

import random

from dnadis.models import ContigClassification
from dnadis.output import scaffolding
from dnadis.output.scaffolding import (
    BRACKET_SLACK,
    OFF_TARGET_MAX_BP,
    _off_target_bp,
    _gap_bracket,
    _index_ref_blocks,
    _place_fragments_in_gaps,
    _RefBlock,
    _scaffold_from_ragtag_agp,
    orientations_from_agp,
    scaffold_chromosomes,
)
from dnadis.output.tsv_output import chimera_fractions
from dnadis.utils.sequence_utils import reverse_complement

REF = "chr1A"


def _row(contig, strand, qstart, qend, ref_start, ref_end, bp=100_000, ref=REF):
    """A 19-field macro block tuple as chain parsing writes it."""
    return (contig, 0, ref, ref[:-1], ref[-1], strand, 1, qstart, qend, qend - qstart,
            bp, 0, 0, 0.9, 0.0, 1, 0, ref_start, ref_end)


def _blk(strand, qstart, qend, ref_start, ref_end, bp=100_000):
    return _RefBlock(qstart, qend, strand, ref_start, ref_end, bp)


def _seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _clf(name, new_name, cls, **kw):
    return ContigClassification(
        original_name=name, new_name=new_name, classification=cls,
        reversed=False, cobiont_taxid=None, cobiont_sci=None,
        assigned_ref_id=REF, ref_gene_proportion=None, contig_len=0,
        query_subgenome_grp=1, seq_identity_vs_ref=0.9, **kw,
    )


# --- Gap brackets ------------------------------------------------------------

def test_bracket_from_nearest_flanks():
    blocks = [
        _blk("+", 0, 400, 1_000_000, 1_400_000),
        _blk("+", 500, 900, 1_500_000, 1_900_000),
        _blk("+", 1_100, 1_500, 2_300_000, 2_700_000),
    ]
    br = _gap_bracket("a", (1_000, 1_010), blocks)
    assert (br.left_ref, br.right_ref) == (1_900_000, 2_300_000)
    assert (br.left_strand, br.right_strand) == ("+", "+")
    # On the minus strand, the gap-side ends are the blocks' reference starts.
    minus = [_blk("-", 0, 900, 2_300_000, 3_000_000), _blk("-", 1_100, 1_500, 1_000_000, 1_900_000)]
    br = _gap_bracket("a", (1_000, 1_010), minus)
    assert (br.left_ref, br.right_ref) == (2_300_000, 1_900_000)


def test_closed_gaps_have_no_bracket():
    # A block aligns across the gap.
    assert _gap_bracket("a", (1_000, 1_010), [
        _blk("+", 0, 400, 0, 400_000), _blk("+", 500, 1_500, 500_000, 1_500_000),
    ]) is None
    # No flank on one side.
    assert _gap_bracket("a", (1_000, 1_010), [_blk("+", 0, 900, 0, 900_000)]) is None
    # Same-strand flanks running backwards on the reference.
    assert _gap_bracket("a", (1_000, 1_010), [
        _blk("+", 0, 900, 1_000_000, 2_000_000),
        _blk("+", 1_100, 1_500, 1_000_000 - BRACKET_SLACK - 1, 1_500_000),
    ]) is None


def test_short_blocks_are_not_flanks():
    blocks = [
        _blk("+", 0, 800, 1_000_000, 1_800_000),
        _blk("+", 850, 950, 5_000_000, 5_010_000, bp=10_000),  # out of place
        _blk("+", 1_100, 1_500, 2_300_000, 2_700_000),
    ]
    br = _gap_bracket("a", (1_000, 1_010), blocks)
    assert br.left_ref == 1_800_000


def test_fit_between_same_strand_flanks():
    br = _gap_bracket("a", (1_000, 1_010), [
        _blk("+", 0, 900, 1_000_000, 1_900_000), _blk("+", 1_100, 1_500, 2_300_000, 2_700_000),
    ])
    assert br.fit(1_950_000, 2_250_000, "+")[0] == "+"
    assert br.fit(1_950_000, 2_250_000, "-")[0] == "-"
    # Within the slack at either side, and not beyond.
    assert br.fit(1_900_000 - BRACKET_SLACK, 2_300_000 + BRACKET_SLACK, "+") is not None
    assert br.fit(1_000_000, 1_500_000, "+") is None


def test_fit_with_flanks_on_opposite_strands():
    # Left flank on the minus strand (an inversion next to the gap), right
    # flank on the plus strand at reference 2,300,000.
    br = _gap_bracket("a", (1_000, 1_010), [
        _blk("-", 0, 900, 1_000_000, 1_500_000), _blk("+", 1_100, 1_500, 2_300_000, 2_700_000),
    ])
    # Ends just before the right flank: continues it, oriented with it.
    assert br.fit(2_100_000, 2_299_000, "+")[0] == "+"
    assert br.fit(2_100_000, 2_299_000, "-")[0] == "-"
    # Starts just below the left flank's reference start: continues the
    # minus-strand left flank.
    assert br.fit(800_000, 990_000, "+")[0] == "-"
    # Continues neither flank.
    assert br.fit(1_700_000, 1_800_000, "+") is None


# --- Placement ---------------------------------------------------------------

def _place(anchor_rows, frag_rows, gaps, lengths, chimeric=None):
    blocks = _index_ref_blocks(anchor_rows + frag_rows)
    frags = sorted({r[0] for r in frag_rows} | {f for f in lengths if f != "a"})
    return _place_fragments_in_gaps(["a"], frags, REF, gaps, blocks, lengths, chimeric=chimeric)


def test_fragment_placed_into_its_gap():
    rows = [_row("a", "+", 0, 900, 1_000_000, 1_900_000), _row("a", "+", 1_100, 2_000, 2_300_000, 3_000_000)]
    comps, placed, unplaced = _place(
        rows, [_row("f", "-", 0, 300, 1_950_000, 2_250_000)],
        {"a": [(1_000, 1_010)]}, {"a": 2_000, "f": 300},
    )
    assert comps == {"a": [
        ("W", "a", 1, 1_000, "+"), ("N", "", 0, 0, ""), ("W", "f", 1, 300, "-"),
        ("N", "", 0, 0, ""), ("W", "a", 1_011, 2_000, "+"),
    ]}
    assert set(placed) == {"f"} and unplaced == {}


def test_fragment_fitting_two_gaps_or_none_is_not_placed():
    # Two gaps with nothing aligned between them share a bracket.
    rows = [_row("a", "+", 0, 900, 1_000_000, 1_900_000), _row("a", "+", 2_100, 3_000, 2_300_000, 3_000_000)]
    comps, placed, unplaced = _place(
        rows,
        [_row("f", "+", 0, 300, 1_950_000, 2_250_000), _row("g", "+", 0, 300, 5_000_000, 5_300_000)],
        {"a": [(1_000, 1_010), (2_000, 2_010)]},
        {"a": 3_000, "f": 300, "g": 300, "h": 300},
    )
    assert comps == {} and placed == {}
    assert unplaced == {"f": "gap_ambiguous", "g": "no_gap_bracket", "h": "no_gap_bracket"}


def test_fragments_sharing_a_gap_follow_the_reference():
    frags = [
        _row("late", "+", 0, 100, 2_150_000, 2_250_000),
        _row("early", "+", 0, 100, 1_950_000, 2_050_000),
    ]
    lengths = {"a": 2_000, "early": 100, "late": 100}
    # Anchor on the plus strand: early, then late.
    rows = [_row("a", "+", 0, 900, 1_000_000, 1_900_000), _row("a", "+", 1_100, 2_000, 2_300_000, 3_000_000)]
    comps, _, _ = _place(rows, frags, {"a": [(1_000, 1_010)]}, lengths)
    assert [p[1] for p in comps["a"] if p[0] == "W"] == ["a", "early", "late", "a"]
    # Anchor on the minus strand: late, then early, both reverse-complemented.
    rows = [_row("a", "-", 0, 900, 2_300_000, 3_000_000), _row("a", "-", 1_100, 2_000, 1_000_000, 1_900_000)]
    comps, _, _ = _place(rows, frags, {"a": [(1_000, 1_010)]}, lengths)
    assert [(p[1], p[4]) for p in comps["a"] if p[0] == "W"] == [
        ("a", "+"), ("late", "-"), ("early", "-"), ("a", "+"),
    ]


def test_overlapping_fragments_in_one_gap_are_not_placed():
    rows = [_row("a", "+", 0, 900, 1_000_000, 1_900_000), _row("a", "+", 1_100, 2_000, 2_300_000, 3_000_000)]
    comps, placed, unplaced = _place(
        rows,
        [_row("f", "+", 0, 300, 1_950_000, 2_250_000), _row("g", "+", 0, 300, 2_000_000, 2_200_000)],
        {"a": [(1_000, 1_010)]}, {"a": 2_000, "f": 300, "g": 300},
    )
    assert comps == {} and placed == {}
    assert unplaced == {"f": "gap_ambiguous", "g": "gap_ambiguous"}


def test_fragments_beside_opposite_strand_flanks_keep_flank_order():
    # Left flank +, right flank - (inverted relative to the reference): one
    # fragment continues each, and they are laid out left-flank first.
    rows = [_row("a", "+", 0, 900, 1_000_000, 1_900_000), _row("a", "-", 1_100, 2_000, 300_000, 900_000)]
    comps, placed, _ = _place(
        rows,
        [_row("r", "+", 0, 100, 920_000, 1_000_000), _row("l", "-", 0, 100, 1_910_000, 2_100_000)],
        {"a": [(1_000, 1_010)]}, {"a": 2_000, "l": 100, "r": 100},
    )
    assert [(p[1], p[4]) for p in comps["a"] if p[0] == "W"] == [
        ("a", "+"), ("l", "-"), ("r", "-"), ("a", "+"),
    ]


# --- Fragments with sequence from other chromosomes --------------------------

# Anchor "a" with one gap; reference 1.95-2.25 Mb lies between its flanks.
_ANCHOR = [
    _row("a", "+", 0, 900_000, 1_000_000, 1_900_000),
    _row("a", "+", 1_100_000, 2_000_000, 2_300_000, 3_000_000),
]
_ANCHOR_GAPS = {"a": [(1_000_000, 1_000_100)]}


def _fitting_fragment(length=400_000):
    """A fragment whose first 300 kb fits the anchor's gap."""
    return [_row("f", "+", 0, 300_000, 1_950_000, 2_250_000)], {"a": 2_000_000, "f": length}


def test_off_target_bp_counts_sequence_on_other_references_only():
    rows, _ = _fitting_fragment()
    rows += [
        _row("f", "+", 300_000, 330_000, 0, 30_000, ref="chr5A"),
        _row("f", "-", 340_000, 380_000, 0, 40_000, ref="chr9B"),
        # Aligned to another reference where it also aligns to its own (a
        # homeolog): not off-target.
        _row("f", "+", 0, 200_000, 1_950_000, 2_150_000, ref="chr1B"),
        # Overlaps the chr5A block: counted once.
        _row("f", "+", 320_000, 335_000, 500_000, 515_000, ref="chr2A"),
    ]
    by_ref = {ref: blocks for (_, ref), blocks in _index_ref_blocks(rows).items()}
    assert _off_target_bp(by_ref, REF) == 35_000 + 40_000
    assert _off_target_bp({REF: by_ref[REF]}, REF) == 0


def test_fragment_with_other_chromosome_pieces_is_not_placed():
    # Pieces between the fragment's internal gaps map to two other
    # references, together above the limit.
    rows, lengths = _fitting_fragment()
    rows += [
        _row("f", "+", 300_000, 330_000, 0, 30_000, ref="chr5A"),
        _row("f", "-", 350_000, 380_000, 0, 30_000, ref="chr9B"),
    ]
    comps, placed, unplaced = _place(_ANCHOR, rows, _ANCHOR_GAPS, lengths)
    assert comps == {} and placed == {}
    assert unplaced == {"f": "chimeric_fragment"}


def test_off_target_end_block_at_the_limit_is_not_placed():
    # A gapless fragment whose end, which would face the gap flank, maps to
    # another reference.
    rows, lengths = _fitting_fragment()
    end = [_row("f", "+", 300_000, 300_000 + OFF_TARGET_MAX_BP, 0, OFF_TARGET_MAX_BP, ref="chr5A")]
    _, placed, unplaced = _place(_ANCHOR, rows + end, _ANCHOR_GAPS, lengths)
    assert placed == {} and unplaced == {"f": "chimeric_fragment"}
    # Just below the limit, it is placed.
    end = [_row("f", "+", 300_000, 299_999 + OFF_TARGET_MAX_BP, 0, OFF_TARGET_MAX_BP, ref="chr5A")]
    _, placed, unplaced = _place(_ANCHOR, rows + end, _ANCHOR_GAPS, lengths)
    assert set(placed) == {"f"} and unplaced == {}


def test_chimeric_flag_keeps_fragment_out_of_gaps():
    rows, lengths = _fitting_fragment()
    comps, placed, unplaced = _place(_ANCHOR, rows, _ANCHOR_GAPS, lengths, chimeric={"f"})
    assert comps == {} and placed == {}
    assert unplaced == {"f": "chimeric_fragment"}


def test_chimera_fractions():
    # Best reference 70% of aligned bp, second 25%: chimeric.
    assert chimera_fractions(1_000, 2, 700, 250, 0.8, 0.2) == (0.7, 0.25)
    # Second reference below its minimum, or a single reference: not.
    assert chimera_fractions(1_000, 2, 700, 150, 0.8, 0.2) is None
    assert chimera_fractions(1_000, 1, 500, 0, 0.8, 0.2) is None
    assert chimera_fractions(0, 2, 0, 0, 0.8, 0.2) is None


# --- Scaffolds ---------------------------------------------------------------

def _gapped_inputs(tmp_path, t2t=True):
    left, right = _seq(900, 1), _seq(990, 2)
    seqs = {
        "a": left + "N" * 10 + right,
        "f": _seq(300, 3),
        "g": _seq(200, 4),
    }
    qfa = tmp_path / "q.fa"
    qfa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    rfa = tmp_path / "r.fa"
    rfa.write_text(f">{REF}\n" + "A" * 3000 + "\n")
    tel = dict(is_chromosome_scale=True, has_5p_telomere=True, has_3p_telomere=True) if t2t else {}
    clfs = [
        _clf("a", "chr1A", "chrom_assigned", n_gaps=1, **tel),
        _clf("f", "chr1A_f1", "chrom_fragment"),
        _clf("g", "chr1A_f2", "chrom_fragment"),
    ]
    best_ref = {c.original_name: REF for c in clfs}
    rows = [
        _row("a", "+", 0, 900, 1_000_000, 1_900_000),
        _row("a", "+", 910, 1_900, 2_300_000, 3_000_000),
        _row("f", "-", 0, 300, 1_950_000, 2_250_000),
        _row("g", "+", 0, 200, 4_000_000, 4_200_000),
    ]
    kwargs = dict(
        query_fasta=qfa, classifications=clfs,
        qr_ref_ranges={("a", REF): (1_000_000, 3_000_000), ("f", REF): (1_950_000, 2_250_000),
                       ("g", REF): (4_000_000, 4_200_000)},
        ref_fasta=rfa, ref_lengths={REF: 5_000_000},
        best_ref=best_ref, contig_refs={k: {v} for k, v in best_ref.items()},
        qr_best_chain_ident={(k, v): 0.9 for k, v in best_ref.items()},
        gap_size=5, query_gaps={"a": [(900, 910)]}, macro_block_rows=rows,
    )
    return seqs, clfs, kwargs


def test_t2t_scaffold_with_fragment_in_its_gap(tmp_path, monkeypatch):
    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs, clfs, kwargs = _gapped_inputs(tmp_path)
    filled = (seqs["a"][:900] + "N" * 5 + reverse_complement(seqs["f"]) + "N" * 5
              + seqs["a"][910:])

    scaffolded, agp, _ = scaffold_chromosomes(
        contig_orientations={"a": False, "f": False, "g": False},
        work_dir=tmp_path / "s1", **kwargs,
    )
    assert scaffolded == {"chr1A": filled, "chr1A_f2": seqs["g"]}
    rows = [line.split("\t") for line in agp if line.startswith("chr1A\t")]
    assert [(r[4], r[5], r[6], r[7], r[8]) if r[4] == "W" else (r[4], r[5]) for r in rows] == [
        ("W", "a", "1", "900", "+"), ("N", "5"), ("W", "f", "1", "300", "-"),
        ("N", "5"), ("W", "a", "911", "1900", "+"),
    ]
    assert rows[-1][2] == str(len(filled))
    # The fragment that fits no gap is written on its own, with the reason.
    assert {c.original_name: c.scaffold_exclusion for c in clfs} == {
        "a": None, "f": None, "g": "no_gap_bracket",
    }

    # A reversed anchor reverses the whole composite.
    for c in clfs:
        c.scaffold_exclusion = None
    scaffolded, agp, _ = scaffold_chromosomes(
        contig_orientations={"a": True, "f": False, "g": False},
        work_dir=tmp_path / "s2", **kwargs,
    )
    assert scaffolded["chr1A"] == reverse_complement(filled)
    assert orientations_from_agp(agp) == {"a": True, "f": False, "g": False}


def test_t2t_scaffold_keeps_chimeric_fragment_out_of_its_gap(tmp_path, monkeypatch):
    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs, clfs, kwargs = _gapped_inputs(tmp_path)
    scaffolded, agp, _ = scaffold_chromosomes(
        contig_orientations={"a": False, "f": False, "g": False},
        work_dir=tmp_path / "s", chimeric={"f"}, **kwargs,
    )
    # The gap stays open; the fragment is written on its own.
    assert scaffolded == {"chr1A": seqs["a"], "chr1A_f1": seqs["f"], "chr1A_f2": seqs["g"]}
    assert {c.original_name: c.scaffold_exclusion for c in clfs} == {
        "a": None, "f": "chimeric_fragment", "g": "no_gap_bracket",
    }


def test_gapless_t2t_is_unchanged(tmp_path, monkeypatch):
    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs, clfs, kwargs = _gapped_inputs(tmp_path)
    kwargs["query_gaps"] = {}
    scaffolded, _, _ = scaffold_chromosomes(
        contig_orientations={"a": False, "f": False, "g": False},
        work_dir=tmp_path / "s", **kwargs,
    )
    assert scaffolded == {"chr1A": seqs["a"], "chr1A_f1": seqs["f"], "chr1A_f2": seqs["g"]}
    assert all(c.scaffold_exclusion is None for c in clfs)


def test_builtin_scaffold_expands_filled_anchor(tmp_path, monkeypatch):
    """Without telomeres, the filled anchor is joined to the group's other
    fragments by reference position."""
    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs, clfs, kwargs = _gapped_inputs(tmp_path, t2t=False)
    scaffolded, agp, _ = scaffold_chromosomes(
        contig_orientations={"a": False, "f": False, "g": False},
        work_dir=tmp_path / "s", **kwargs,
    )
    filled = (seqs["a"][:900] + "N" * 5 + reverse_complement(seqs["f"]) + "N" * 5
              + seqs["a"][910:])
    assert scaffolded == {"chr1A": filled + "N" * 5 + seqs["g"]}
    assert [line.split("\t")[5] for line in agp if line.split("\t")[4] == "W"] == ["a", "f", "a", "g"]
    assert all(c.scaffold_exclusion is None for c in clfs)


def test_ragtag_agp_expands_reversed_composite(tmp_path):
    seqs = {"a": _seq(100, 5), "f": _seq(30, 6), "g": _seq(40, 7)}
    composites = {"a": [
        ("W", "a", 1, 50, "+"), ("N", "", 0, 0, ""), ("W", "f", 1, 30, "-"),
        ("N", "", 0, 0, ""), ("W", "a", 61, 100, "+"),
    ]}
    agp = tmp_path / "ragtag.scaffold.agp"
    agp.write_text(
        "chr1A_RagTag\t1\t40\t1\tW\tg\t1\t40\t+\n"
        "chr1A_RagTag\t41\t140\t2\tU\t100\tscaffold\tyes\talign_genus\n"
        "chr1A_RagTag\t141\t230\t3\tW\ta\t1\t90\t-\n"
    )
    seq, lines, unplaced = _scaffold_from_ragtag_agp(
        agp, seqs, "chr1A", 5, "chr1A", composites=composites,
    )
    assert unplaced == []
    assert seq == (seqs["g"] + "N" * 5 + reverse_complement(seqs["a"][60:]) + "N" * 5
                   + seqs["f"] + "N" * 5 + reverse_complement(seqs["a"][:50]))
    assert [(f[5], f[6], f[7], f[8]) for f in (l.split("\t") for l in lines) if f[4] == "W"] == [
        ("g", "1", "40", "+"), ("a", "61", "100", "-"), ("f", "1", "30", "+"), ("a", "1", "50", "-"),
    ]
