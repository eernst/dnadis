#!/usr/bin/env python3
"""Reference-guided scaffolding for dnadis.

Produces chromosome-scale pseudomolecules by:
1. Grouping contigs by assigned reference chromosome
2. Placing fragments into scaffold gaps (runs of >=10 N) of the group's
   chromosome contigs where the reference positions of the gap flanks
   bracket them
3. Running RagTag scaffold per chromosome (if available)
4. Falling back to built-in reference-position-based ordering otherwise
5. Producing AGP 2.0 output describing scaffold structure
"""
from __future__ import annotations

import shutil
import statistics
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from dnadis.models import ContigClassification
from dnadis.utils.io_utils import have_exe, merge_intervals
from dnadis.utils.logging_config import get_logger
from dnadis.utils.sequence_utils import (
    read_fasta_sequences,
    reverse_complement,
    write_fasta,
    write_filtered_fasta,
)

logger = get_logger("scaffolding")


# ---------------------------------------------------------------------------
# AGP output
# ---------------------------------------------------------------------------

def orientations_from_agp(agp_lines: List[str]) -> Dict[str, bool]:
    """Extract per-contig orientation from AGP component lines.

    Returns:
        Dict mapping component_id -> is_reversed (True if AGP orientation is '-').
    """
    result: Dict[str, bool] = {}
    for line in agp_lines:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[4] != "W":
            continue
        component_id = fields[5]
        orientation = fields[8]
        result[component_id] = (orientation == "-")
    return result


def write_agp(agp_lines: List[str], output_path: Path) -> None:
    """Write AGP 2.0 file.

    Args:
        agp_lines: List of tab-delimited AGP lines (without trailing newline).
        output_path: Output file path.
    """
    with output_path.open("w") as fh:
        fh.write("##agp-version\t2.0\n")
        for line in agp_lines:
            fh.write(line + "\n")


def _agp_component_line(
    scaffold_name: str,
    scaffold_start: int,
    scaffold_end: int,
    part_number: int,
    component_id: str,
    component_start: int,
    component_end: int,
    orientation: str,
) -> str:
    """Create an AGP W (component) line."""
    return "\t".join(map(str, [
        scaffold_name,
        scaffold_start,
        scaffold_end,
        part_number,
        "W",
        component_id,
        component_start,
        component_end,
        orientation,
    ]))


def _agp_gap_line(
    scaffold_name: str,
    scaffold_start: int,
    scaffold_end: int,
    part_number: int,
    gap_length: int,
) -> str:
    """Create an AGP N (gap) line."""
    return "\t".join(map(str, [
        scaffold_name,
        scaffold_start,
        scaffold_end,
        part_number,
        "N",
        gap_length,
        "scaffold",
        "yes",
        "align_genus",
    ]))


# A scaffold is laid out as a list of parts: ("W", component, beg, end,
# orientation) with 1-based inclusive component coordinates, or a gap,
# ("N", "", 0, 0, "").  A contig with fragments placed into its input gaps
# is a "composite": its parts in its own forward orientation, which replace
# its full-length W part wherever the contig is laid out.
_Part = Tuple[str, str, int, int, str]
_GAP_PART: _Part = ("N", "", 0, 0, "")
_FLIP = {"+": "-", "-": "+"}


def _flip_parts(parts: List[_Part]) -> List[_Part]:
    """Parts of the reverse complement of a layout."""
    return [
        (k, c, b, e, _FLIP[o] if k == "W" else o) for (k, c, b, e, o) in reversed(parts)
    ]


def _expand_composites(
    parts: List[_Part], composites: Dict[str, List[_Part]],
) -> List[_Part]:
    """Replace each W part of a composite contig with the composite's parts."""
    if not composites:
        return parts
    out: List[_Part] = []
    for part in parts:
        sub = composites.get(part[1]) if part[0] == "W" else None
        if sub is None:
            out.append(part)
        else:
            out.extend(sub if part[4] == "+" else _flip_parts(sub))
    return out


def _render_parts(
    scaffold_name: str,
    parts: List[_Part],
    query_seqs: Dict[str, str],
    gap_size: int,
) -> Tuple[str, List[str]]:
    """Scaffold sequence and AGP lines for a layout.

    Components missing from ``query_seqs`` are skipped; gaps are written as
    ``gap_size`` N, never at either end and never two in a row.
    """
    pieces: List[str] = []
    rows: List[Tuple[str, int, int, str, int, int, str]] = []
    for kind, comp, beg, end, orient in parts:
        if kind == "N":
            if gap_size > 0 and rows and rows[-1][0] == "W":
                rows.append(("N", 0, 0, "", 0, 0, ""))
                pieces.append("N" * gap_size)
            continue
        seq = query_seqs.get(comp, "")
        if not seq or beg > end:
            continue
        sub = seq[beg - 1:end]
        if orient == "-":
            sub = reverse_complement(sub)
        rows.append(("W", 0, 0, comp, beg, end, orient))
        pieces.append(sub)
    if rows and rows[-1][0] == "N":
        rows.pop()
        pieces.pop()

    agp_lines: List[str] = []
    pos = 1
    for part_num, (row, piece) in enumerate(zip(rows, pieces), start=1):
        end_pos = pos + len(piece) - 1
        if row[0] == "N":
            agp_lines.append(_agp_gap_line(scaffold_name, pos, end_pos, part_num, gap_size))
        else:
            _, _, _, comp, beg, end, orient = row
            agp_lines.append(_agp_component_line(
                scaffold_name, pos, end_pos, part_num, comp, beg, end, orient,
            ))
        pos = end_pos + 1
    return "".join(pieces), agp_lines


# ---------------------------------------------------------------------------
# Placement of fragments into input scaffold gaps
# ---------------------------------------------------------------------------

# Reference slack at each side of a gap's bracket.  Flank and fragment
# positions come from macro blocks, whose ends are coarse when the reference
# is divergent.
BRACKET_SLACK = 50_000

# Macro blocks shorter than this are not used as gap flanks: in divergent,
# repeat-rich stretches, short blocks often align far from their neighbours
# on the reference.
MIN_FLANK_BP = 20_000

# A fragment with at least this much of its sequence in macro blocks on other
# reference chromosomes (and in none on the group's) is not placed into a
# gap, since placing it would insert that sequence into the chromosome.  It
# equals BRACKET_SLACK: off-target sequence shorter than that is within the
# reference imprecision placement accepts at each side of a gap.
OFF_TARGET_MAX_BP = BRACKET_SLACK


@dataclass(frozen=True)
class _RefBlock:
    """A macro block of one contig on one reference chromosome."""

    qstart: int
    qend: int
    strand: str
    ref_start: int
    ref_end: int
    bp: int


def _index_ref_blocks(macro_block_rows) -> Dict[Tuple[str, str], List[_RefBlock]]:
    """(contig original name, ref_id) -> macro blocks with reference coordinates.

    Rows without reference coordinates (older 17-field tuples) are skipped.
    """
    out: Dict[Tuple[str, str], List[_RefBlock]] = defaultdict(list)
    for row in macro_block_rows or ():
        if len(row) < 19:
            continue
        ref_start, ref_end = row[17], row[18]
        if ref_start is None or ref_end is None:
            continue
        out[(row[0], row[2])].append(_RefBlock(
            qstart=int(row[7]), qend=int(row[8]), strand=row[5],
            ref_start=int(ref_start), ref_end=int(ref_end), bp=int(row[10]),
        ))
    return dict(out)


@dataclass(frozen=True)
class _GapBracket:
    """Reference positions at the two sides of an input gap in an anchor.

    ``left_ref``/``right_ref`` are the reference positions at the gap-side
    ends of the nearest flank blocks before and after the gap (in anchor
    coordinates); ``left_strand``/``right_strand`` are those blocks' strands.
    """

    anchor: str
    gap: Tuple[int, int]
    left_ref: int
    right_ref: int
    left_strand: str
    right_strand: str

    @property
    def lo(self) -> int:
        return min(self.left_ref, self.right_ref)

    @property
    def hi(self) -> int:
        return max(self.left_ref, self.right_ref)

    def fit(
        self, fs: int, fe: int, strand: str,
    ) -> Optional[Tuple[str, Tuple[int, int]]]:
        """Whether a fragment with reference footprint [fs, fe) and dominant
        strand ``strand`` fits this gap: (orientation in the anchor, sort key
        for fragments sharing the gap), or None.

        With flanks on one strand, the footprint must lie between the flanks'
        reference positions (``BRACKET_SLACK`` at each side).  With flanks on
        opposite strands (a rearrangement next to the gap), the fragment must
        continue the reference from exactly one flank: lie on the side of its
        reference position that the flank runs toward, starting within
        ``BRACKET_SLACK`` of it.  It is oriented with that flank and laid out
        next to it.
        """
        if self.left_strand == self.right_strand:
            if not (self.lo - BRACKET_SLACK <= fs and fe <= self.hi + BRACKET_SLACK):
                return None
            orient = "+" if strand == self.left_strand else "-"
            mid = fs + fe
            return orient, (0, mid if self.right_ref >= self.left_ref else -mid)

        def continues(ref: int, toward_higher: bool) -> Optional[int]:
            dist = fs - ref if toward_higher else ref - fe
            return dist if -BRACKET_SLACK <= dist <= BRACKET_SLACK else None

        # A + block runs toward higher reference positions; the sequence after
        # a left flank continues its direction, the sequence before a right
        # flank precedes it.
        d_left = continues(self.left_ref, self.left_strand == "+")
        d_right = continues(self.right_ref, self.right_strand == "-")
        if (d_left is None) == (d_right is None):
            return None
        if d_left is not None:
            return ("+" if strand == self.left_strand else "-"), (0, d_left)
        return ("+" if strand == self.right_strand else "-"), (1, -d_right)


def _gap_bracket(
    anchor: str, gap: Tuple[int, int], blocks: List[_RefBlock],
) -> Optional[_GapBracket]:
    """The bracket of one gap, or None when the gap leaves no reference
    sequence open: a flank block aligns across it, it has no flank block on
    one side, or flanks on the same strand run backwards on the reference
    by more than ``BRACKET_SLACK`` (the sequence on both sides already
    covers what lies between).  Blocks
    shorter than ``MIN_FLANK_BP`` are ignored.
    """
    gs, ge = gap
    left: Optional[_RefBlock] = None
    right: Optional[_RefBlock] = None
    for b in blocks:
        if b.bp < MIN_FLANK_BP:
            continue
        if b.qstart < gs and b.qend > ge:
            return None
        if b.qstart < gs:
            if left is None or b.qend > left.qend:
                left = b
        elif right is None or b.qstart < right.qstart:
            right = b
    if left is None or right is None:
        return None
    left_ref = left.ref_end if left.strand == "+" else left.ref_start
    right_ref = right.ref_start if right.strand == "+" else right.ref_end
    if left.strand == right.strand:
        step = right_ref - left_ref if left.strand == "+" else left_ref - right_ref
        if step < -BRACKET_SLACK:
            return None
    return _GapBracket(anchor, gap, left_ref, right_ref, left.strand, right.strand)


def _dominant_strand(blocks: List[_RefBlock]) -> str:
    plus = sum(b.bp for b in blocks if b.strand == "+")
    minus = sum(b.bp for b in blocks if b.strand == "-")
    return "+" if plus >= minus else "-"


def _off_target_bp(blocks_by_ref: Dict[str, List[_RefBlock]], ref_id: str) -> int:
    """Contig bp in macro blocks on references other than ``ref_id`` and in
    none on ``ref_id`` (block spans on the contig), from one contig's blocks
    by reference."""
    own = [(b.qstart, b.qend) for b in blocks_by_ref.get(ref_id, ())]
    other = [
        (b.qstart, b.qend)
        for r, blocks in blocks_by_ref.items() if r != ref_id
        for b in blocks
    ]
    if not other:
        return 0
    return merge_intervals(own + other)[1] - merge_intervals(own)[1]


def _place_fragments_in_gaps(
    anchors: List[str],
    fragments: List[str],
    ref_id: str,
    query_gaps: Dict[str, List[Tuple[int, int]]],
    ref_blocks: Dict[Tuple[str, str], List[_RefBlock]],
    query_lengths: Dict[str, int],
    chimeric: Optional[Set[str]] = None,
) -> Tuple[Dict[str, List[_Part]], Dict[str, str], Dict[str, str]]:
    """Place fragments into the input scaffold gaps of anchor contigs.

    Fragments in ``chimeric`` (the contig_summary ``chimeric`` flag), and
    fragments with at least ``OFF_TARGET_MAX_BP`` in macro blocks on other
    references only (``_off_target_bp``), are not placed: their sequence
    from other chromosomes would be inserted with them.

    A fragment's reference footprint on ``ref_id`` runs from its first to its
    last macro block, and its strand is the dominant one.  It is placed into
    a gap when it fits that gap (``_GapBracket.fit``) and no other gap of the
    group's anchors.  Fragments sharing a gap are laid out along the
    reference from the left flank to the right (with opposite-strand flanks,
    those continuing the left flank first); if two of them overlap on the
    reference by at least half the shorter footprint, none of that gap's
    fragments are placed.  Each filled gap becomes fragment(s) between
    gaps of the scaffolder's gap size.

    Args:
        anchors: chrom_assigned contigs of the group (original names).
        fragments: Other group members that may be placed.
        ref_id: The group's reference chromosome.
        query_gaps: Contig -> 0-based half-open N runs.
        ref_blocks: From ``_index_ref_blocks``.
        query_lengths: Contig -> length.

    Returns:
        (composites, placed, unplaced) where ``composites`` maps an anchor to
        its parts with fragments in its gaps, ``placed`` maps each placed
        fragment to a log description, and ``unplaced`` maps each fragment
        that was not placed to ``chimeric_fragment``, ``gap_ambiguous`` (it
        fits several gaps, or overlaps another fragment fitting the same gap)
        or ``no_gap_bracket``.
    """
    brackets: List[_GapBracket] = []
    for anchor in anchors:
        blocks = ref_blocks.get((anchor, ref_id), [])
        for gap in query_gaps.get(anchor, ()):
            br = _gap_bracket(anchor, gap, blocks)
            if br is not None:
                brackets.append(br)

    footprint: Dict[str, Tuple[int, int]] = {}
    unplaced: Dict[str, str] = {}
    by_gap: Dict[Tuple[str, Tuple[int, int]], List[Tuple[Tuple[int, int], str, str]]] = (
        defaultdict(list)
    )
    bracket_of: Dict[Tuple[str, Tuple[int, int]], _GapBracket] = {}
    chimeric = chimeric or set()
    frag_set = set(fragments)
    frag_blocks: Dict[str, Dict[str, List[_RefBlock]]] = defaultdict(dict)
    for (contig, ref), blocks in ref_blocks.items():
        if contig in frag_set:
            frag_blocks[contig][ref] = blocks
    for frag in fragments:
        off_target = _off_target_bp(frag_blocks.get(frag, {}), ref_id)
        if frag in chimeric or off_target >= OFF_TARGET_MAX_BP:
            unplaced[frag] = "chimeric_fragment"
            logger.info(
                f"Not placing {frag} into a gap: "
                + ("flagged chimeric" if frag in chimeric else
                   f"{off_target:,} bp aligned only to other references")
            )
            continue
        blocks = ref_blocks.get((frag, ref_id), [])
        hits = []
        if blocks:
            fs = min(b.ref_start for b in blocks)
            fe = max(b.ref_end for b in blocks)
            footprint[frag] = (fs, fe)
            strand = _dominant_strand(blocks)
            for br in brackets:
                fit = br.fit(fs, fe, strand)
                if fit is not None:
                    hits.append((br, fit))
        if len(hits) != 1:
            unplaced[frag] = "gap_ambiguous" if hits else "no_gap_bracket"
            continue
        br, (orient, key) = hits[0]
        by_gap[(br.anchor, br.gap)].append((key, frag, orient))
        bracket_of[(br.anchor, br.gap)] = br

    def _overlap(a: Tuple[int, int], b: Tuple[int, int]) -> float:
        ov = min(a[1], b[1]) - max(a[0], b[0])
        shorter = min(a[1] - a[0], b[1] - b[0])
        return ov / shorter if shorter > 0 else 0.0

    fills: Dict[str, List[Tuple[Tuple[int, int], List[Tuple[str, str]]]]] = defaultdict(list)
    placed: Dict[str, str] = {}
    for (anchor, gap), entries in sorted(by_gap.items()):
        frags = [f for _, f, _ in entries]
        if any(
            _overlap(footprint[a], footprint[b]) >= 0.5
            for i, a in enumerate(frags) for b in frags[i + 1:]
        ):
            for f in frags:
                unplaced[f] = "gap_ambiguous"
            continue
        br = bracket_of[(anchor, gap)]
        oriented = [(f, orient) for _, f, orient in sorted(entries)]
        for f, orient in oriented:
            fs, fe = footprint[f]
            placed[f] = (
                f"gap {anchor}:{gap[0] + 1}-{gap[1]} (flanks at reference "
                f"{br.left_ref + 1:,} {br.left_strand} / {br.right_ref + 1:,} "
                f"{br.right_strand}; fragment {fs + 1:,}-{fe:,}, placed {orient})"
            )
        fills[anchor].append((gap, oriented))

    composites: Dict[str, List[_Part]] = {}
    for anchor, gap_fills in fills.items():
        parts: List[_Part] = []
        prev = 0
        for (gs, ge), oriented in sorted(gap_fills):
            parts.append(("W", anchor, prev + 1, gs, "+"))
            for f, orient in oriented:
                parts.append(_GAP_PART)
                parts.append(("W", f, 1, query_lengths[f], orient))
            parts.append(_GAP_PART)
            prev = ge
        parts.append(("W", anchor, prev + 1, query_lengths[anchor], "+"))
        composites[anchor] = parts
    return composites, placed, unplaced


# ---------------------------------------------------------------------------
# Contig grouping
# ---------------------------------------------------------------------------

def _group_contigs_by_ref_chrom(
    classifications: List[ContigClassification],
    best_ref: Dict[str, str],
    contig_refs: Dict[str, Set[str]],
) -> Dict[str, List[str]]:
    """Group contigs by target reference chromosome for scaffolding.

    Chromosome-assigned contigs go to their assigned ref.
    Candidate pool (chrom_fragment only) goes to their assigned_ref_id —
    these are unique chromosomal arms too short to be full chromosomes.
    chrom_debris is excluded: it is redundant with already-placed sequence,
    so scaffolding it would duplicate sequence in the pseudomolecule.

    Returns:
        Dict mapping ref_id -> list of contig original names.
    """
    groups: Dict[str, List[str]] = defaultdict(list)

    # First pass: assigned chromosomes
    assigned_contigs: Set[str] = set()
    for clf in classifications:
        if clf.classification == "chrom_assigned" and clf.assigned_ref_id:
            groups[clf.assigned_ref_id].append(clf.original_name)
            assigned_contigs.add(clf.original_name)

    # Second pass: only chrom_fragment contigs are scaffolding candidates
    for clf in classifications:
        if clf.original_name in assigned_contigs:
            continue
        if clf.classification != "chrom_fragment":
            continue

        ref_id = clf.assigned_ref_id or best_ref.get(clf.original_name, "")
        if not ref_id:
            continue

        # Only include if this ref is already a scaffolding target
        if ref_id in groups:
            groups[ref_id].append(clf.original_name)

    return dict(groups)


def _compute_group_mean_idents(
    groups: Dict[Tuple[str, int], List[str]],
    clf_lookup: Dict[str, ContigClassification],
) -> Dict[Tuple[str, int], float]:
    """Compute mean seq_identity_vs_ref for chrom_assigned contigs in each group.

    Args:
        groups: Dict mapping (ref_id, grp) -> list of contig names.
        clf_lookup: Dict mapping contig name -> ContigClassification.

    Returns:
        Dict mapping (ref_id, grp) -> mean identity.
    """
    result: Dict[Tuple[str, int], float] = {}
    for key, contigs in groups.items():
        idents = []
        for c in contigs:
            clf = clf_lookup.get(c)
            if clf and clf.classification == "chrom_assigned" and clf.seq_identity_vs_ref is not None:
                idents.append(clf.seq_identity_vs_ref)
        if idents:
            result[key] = statistics.mean(idents)
    return result


def _own_assembly_copies(
    classifications: List[ContigClassification],
    best_ref: Dict[str, str],
) -> Dict[str, str]:
    """chrom_fragment contigs whose sequence is already in a chromosome-length
    contig of the same assembly, as original name -> evidence for the log.

    The reference-footprint check that demotes fragments to chrom_debris
    misses such copies when the reference is divergent, and RagTag then
    places them at a chromosome end with near-zero location confidence.
    Near-identical copies count against any chromosome contig.  Half-depth
    (haplotig) copies count only against a contig assigned to the fragment's
    own reference chromosome: elsewhere, half-depth single-copy loci that
    align (gametologs, paralogs of a set assembled apart from a collapsed one)
    are not haplotypes of each other.  A fragment that overlaps a chromosome
    end and extends past it is excluded too: RagTag joins contigs end to end
    without merging an overlap, so the overlapping copy would be duplicated.
    """
    from dnadis.detection.redundancy import DUPLICATE_MIN_IDENTITY, is_own_assembly_copy

    by_name = {c.new_name or c.original_name: c for c in classifications}
    out: Dict[str, str] = {}
    for clf in classifications:
        if clf.classification != "chrom_fragment" or not is_own_assembly_copy(
            clf.redundancy_cov, clf.redundancy_identity,
            clf.depth_ratio, clf.target_depth_ratio,
        ):
            continue
        target = by_name.get(clf.redundancy_target or "")
        ref_id = clf.assigned_ref_id or best_ref.get(clf.original_name, "")
        if (
            clf.redundancy_identity < DUPLICATE_MIN_IDENTITY
            and not (target and ref_id and target.assigned_ref_id == ref_id)
        ):
            continue
        out[clf.original_name] = (
            f"{clf.redundancy_target}:{clf.redundancy_region} "
            f"(coverage={clf.redundancy_cov:.2f}, "
            f"identity={clf.redundancy_identity:.4f}, "
            f"depth ratio={_fmt_ratio(clf.depth_ratio)}, "
            f"target depth ratio={_fmt_ratio(clf.target_depth_ratio)})"
        )
    return out


def _fmt_ratio(value: Optional[float]) -> str:
    return "NA" if value is None else f"{value:.2f}"


def _group_contigs_by_haplotype(
    classifications: List[ContigClassification],
    best_ref: Dict[str, str],
    contig_refs: Dict[str, Set[str]],
    qr_best_chain_ident: Dict[Tuple[str, str], float],
    include_rdna: bool = False,
) -> Dict[Tuple[str, int], List[str]]:
    """Group contigs by (ref_id, query_subgenome_grp) for haplotype-aware scaffolding.

    Pass 1: Chromosome-assigned contigs grouped by (assigned_ref_id, grp).
    Pass 2: chrom_fragment contigs join the group of their query subgenome
    set, or the nearest haplotype group by identity when they have none
    (unique arms too short to be full chromosomes; chrom_debris and extra
    copies of a region every set already spans are excluded as redundant). With
    ``include_rdna``, rDNA contigs that carry an assigned_ref_id (they passed
    the synteny gate before being moved to rDNA) are candidates too; their
    placement rests on rDNA similarity shared between NORs.

    Args:
        classifications: List of ContigClassification objects.
        best_ref: Dict mapping contig -> best ref_id.
        contig_refs: Dict mapping contig -> set of ref_ids with evidence.
        qr_best_chain_ident: Dict of (contig, ref_id) -> best chain identity.
        include_rdna: Also place rDNA contigs that have an assigned_ref_id.

    Returns:
        Dict mapping (ref_id, grp) -> list of contig original names.
    """
    groups: Dict[Tuple[str, int], List[str]] = defaultdict(list)
    clf_lookup = {clf.original_name: clf for clf in classifications}

    # Pass 1: chromosome-assigned contigs.  An extra copy of a region already
    # spanned by every set (a chrom_assigned contig shorter than the
    # chromosome can be one) would duplicate sequence in the pseudomolecule.
    assigned_contigs: Set[str] = set()
    for clf in classifications:
        if clf.subgenome_basis == "extra_copy":
            continue
        if clf.classification == "chrom_assigned" and clf.assigned_ref_id:
            grp = clf.query_subgenome_grp or 1
            groups[(clf.assigned_ref_id, grp)].append(clf.original_name)
            assigned_contigs.add(clf.original_name)

    # Build index: ref_id -> set of groups present
    ref_groups: Dict[str, Set[int]] = defaultdict(set)
    for (ref_id, grp) in groups:
        ref_groups[ref_id].add(grp)

    # Compute mean identity per group (for debris assignment)
    group_mean_idents = _compute_group_mean_idents(groups, clf_lookup)

    # Pass 2: chrom_fragment candidates (plus synteny-placed rDNA if requested)
    for clf in classifications:
        if clf.original_name in assigned_contigs:
            continue
        is_rdna_candidate = (
            include_rdna and clf.classification == "rDNA" and bool(clf.assigned_ref_id)
        )
        if clf.classification != "chrom_fragment" and not is_rdna_candidate:
            continue
        if clf.subgenome_basis == "extra_copy":
            continue

        ref_id = clf.assigned_ref_id or best_ref.get(clf.original_name, "")
        if not ref_id:
            continue

        available_grps = ref_groups.get(ref_id)
        if not available_grps:
            # No chrom_assigned contig for this reference — the chromosome survives
            # only as fragments.  Seed a scaffold group from the fragment's own
            # inferred subgenome group so reference-guided scaffolding can still
            # order these fragments along the reference (fragments of the same
            # ref and group accumulate into one pseudomolecule).
            target_grp = clf.query_subgenome_grp or 1
        elif len(available_grps) == 1:
            # Single haplotype group for this ref - assign directly
            target_grp = next(iter(available_grps))
        elif clf.query_subgenome and clf.query_subgenome_grp in available_grps:
            # Follow the set chosen by subgenome inference, so the scaffold
            # matches the contig name.
            target_grp = clf.query_subgenome_grp
        else:
            # Multiple groups - assign by closest identity
            frag_ident = qr_best_chain_ident.get((clf.original_name, ref_id), 0.0)
            target_grp = min(
                available_grps,
                key=lambda g: abs(group_mean_idents.get((ref_id, g), 0.0) - frag_ident),
            )

        groups[(ref_id, target_grp)].append(clf.original_name)

    return dict(groups)


# ---------------------------------------------------------------------------
# Built-in fallback scaffolder
# ---------------------------------------------------------------------------

def _builtin_scaffold(
    contig_names: List[str],
    sequences: Dict[str, str],
    orientations: Dict[str, bool],
    ref_ranges: Dict[Tuple[str, str], Tuple[int, int]],
    ref_id: str,
    scaffold_name: str,
    gap_size: int,
    composites: Optional[Dict[str, List[_Part]]] = None,
) -> Tuple[str, List[str]]:
    """Built-in scaffolder: sort contigs by reference position midpoint.

    Args:
        contig_names: List of contig names to scaffold.
        sequences: Dict mapping contig name -> DNA sequence.
        orientations: Dict mapping contig name -> True if reversed.
        ref_ranges: Dict mapping (contig, ref_id) -> (ref_min, ref_max).
        ref_id: Reference chromosome ID for range lookups.
        scaffold_name: Scaffold name for AGP output.
        gap_size: Number of Ns between contigs.
        composites: Contigs with fragments placed into their input gaps
            (see ``_place_fragments_in_gaps``), laid out with those fragments.

    Returns:
        Tuple of (scaffolded_sequence, agp_lines).
    """
    # Compute midpoint on reference for sorting
    def ref_midpoint(contig: str) -> float:
        rng = ref_ranges.get((contig, ref_id))
        if rng:
            return (rng[0] + rng[1]) / 2.0
        return float("inf")

    parts: List[_Part] = []
    for contig in sorted(contig_names, key=ref_midpoint):
        seq = sequences.get(contig, "")
        if not seq:
            continue
        if parts:
            parts.append(_GAP_PART)
        orient = "-" if orientations.get(contig, False) else "+"
        parts.append(("W", contig, 1, len(seq), orient))
    return _render_parts(
        scaffold_name, _expand_composites(parts, composites or {}), sequences, gap_size,
    )


# ---------------------------------------------------------------------------
# RagTag integration
# ---------------------------------------------------------------------------

def _run_ragtag_scaffold(
    ref_chr_fasta: Path,
    contigs_fasta: Path,
    work_dir: Path,
    threads: int,
) -> Optional[Path]:
    """Run ragtag.py scaffold for a single reference chromosome.

    Returns:
        Path to ragtag.scaffold.agp, or None on failure.
    """
    # Clean any stale output from previous runs.  RagTag caches intermediate
    # alignments (*.paf) and skips re-alignment when they exist, which causes
    # wrong results when the input contig set has changed.
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "ragtag.py", "scaffold",
        str(ref_chr_fasta),
        str(contigs_fasta),
        "-o", str(work_dir),
        "-t", str(threads),
        "-u",  # Include unplaced contigs
    ]
    logger.info(f"Running RagTag: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            logger.warning(f"RagTag failed (rc={result.returncode}): {result.stderr[:500]}")
            return None
    except subprocess.TimeoutExpired:
        logger.warning("RagTag timed out after 3600s")
        return None
    except FileNotFoundError:
        logger.warning("ragtag.py not found in PATH")
        return None

    # Only the AGP is needed downstream: dnadis assembles the scaffold
    # sequence itself from `query_seqs` using the AGP coordinates.  This
    # avoids ragtag.scaffold.fasta, whose writer (ragtag_agp2fa.py) crashes
    # on pysam >= 0.23 because it calls .decode() on a str.
    agp_path = work_dir / "ragtag.scaffold.agp"
    if agp_path.exists():
        return agp_path
    return None


def _scaffold_from_ragtag_agp(
    agp_path: Path,
    query_seqs: Dict[str, str],
    scaffold_name: str,
    gap_size: int,
    ref_seq_key: str,
    composites: Optional[Dict[str, List[_Part]]] = None,
) -> Tuple[Optional[str], List[str], List[str]]:
    """Build scaffold sequence + dnadis-style AGP lines from RagTag's AGP.

    Bypasses RagTag's own FASTA writer (``ragtag_agp2fa.py``), which is
    broken on pysam >= 0.23 (it calls ``.decode('utf-8')`` on a ``str``).
    Reads only the rows of the placed scaffold, ``<ref_seq_key>_RagTag``,
    and reconstructs the sequence from the in-memory contig map using the
    AGP component coordinates and orientations.  Gap rows in the RagTag
    AGP are normalized to the caller's ``gap_size`` (RagTag's default 100
    is preserved if the caller passes 100).

    Args:
        agp_path: Path to ``ragtag.scaffold.agp``.
        query_seqs: Read-only dict of all query sequences.
        scaffold_name: Output scaffold name (replaces RagTag's
            ``<refchr>_RagTag`` name).
        gap_size: Number of Ns between adjacent components.
        ref_seq_key: Name of the reference sequence RagTag scaffolded
            against.  With ``-u``, RagTag writes each contig it cannot place
            as its own scaffold, also suffixed ``_RagTag``, and input contig
            names may already end in ``_RagTag``, so only the exact placed
            scaffold name identifies placed rows.
        composites: Contigs with fragments placed into their input gaps.
            RagTag scaffolded their filled sequences; each of their rows is
            expanded into the composite's parts.

    Returns:
        Tuple of (scaffold_seq, agp_lines, unplaced).
        - ``scaffold_seq`` is None if no placed components were resolved.
        - ``agp_lines`` are dnadis-style AGP rows for the placed scaffold.
        - ``unplaced`` lists component names that RagTag emitted as their
          own scaffolds (i.e. could not place under ``-u``).
    """
    placed_name = f"{ref_seq_key}_RagTag"
    placed: List[_Part] = []
    unplaced: List[str] = []
    with agp_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            scaf, kind = fields[0], fields[4]
            if scaf != placed_name:
                if kind == "W":
                    unplaced.append(fields[5])
                continue
            if kind == "W":
                placed.append(("W", fields[5], int(fields[6]), int(fields[7]), fields[8]))
            elif kind in {"N", "U"}:
                placed.append(_GAP_PART)

    if not any(r[0] == "W" for r in placed):
        return None, [], unplaced

    seq, agp_lines = _render_parts(
        scaffold_name, _expand_composites(placed, composites or {}), query_seqs, gap_size,
    )
    if not seq:
        return None, [], unplaced
    return seq, agp_lines, unplaced


def _parse_ragtag_confidence(
    confidence_path: Path,
) -> Dict[str, Tuple[float, float, float]]:
    """Parse RagTag confidence.txt -> {contig: (grouping, location, orientation)}.

    Args:
        confidence_path: Path to ragtag.scaffold.confidence.txt.

    Returns:
        Dict mapping contig name -> (grouping_conf, location_conf, orientation_conf).
    """
    result: Dict[str, Tuple[float, float, float]] = {}
    if not confidence_path.exists():
        return result
    with confidence_path.open() as fh:
        header = fh.readline()  # skip header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            contig = fields[0]
            try:
                grouping = float(fields[1])
                location = float(fields[2])
                orientation = float(fields[3])
            except (ValueError, IndexError):
                continue
            result[contig] = (grouping, location, orientation)
    return result


def _emit_trivial_scaffold(
    scaffold_name: str,
    contig: str,
    query_seqs: Dict[str, str],
    contig_orientations: Dict[str, bool],
    all_scaffolded: Dict[str, str],
    all_agp: List[str],
    log_label: str = "trivial scaffold",
    composites: Optional[Dict[str, List[_Part]]] = None,
    gap_size: int = 100,
) -> None:
    """Emit a one-contig scaffold (with any fragments placed into its gaps):
    orient, write to dict, append AGP lines."""
    seq = query_seqs.get(contig, "")
    if not seq:
        return
    orient = "-" if contig_orientations.get(contig, False) else "+"
    parts = _expand_composites([("W", contig, 1, len(seq), orient)], composites or {})
    seq, agp_lines = _render_parts(scaffold_name, parts, query_seqs, gap_size)
    all_scaffolded[scaffold_name] = seq
    all_agp.extend(agp_lines)
    n_filled = sum(1 for p in parts if p[0] == "W" and p[1] != contig)
    extra = f" with {n_filled} fragment(s) in its gaps" if n_filled else ""
    logger.info(f"{scaffold_name}: {log_label} {contig}{extra}")


# ---------------------------------------------------------------------------
# Per-group scaffolding helper (for parallel execution)
# ---------------------------------------------------------------------------

@dataclass
class _ScaffoldGroupResult:
    """Result from scaffolding a single multi-contig group."""

    scaffold_name: str
    scaffold_seq: Optional[str] = None
    agp_lines: List[str] = field(default_factory=list)
    confidences: Dict[str, Tuple[float, float, float]] = field(default_factory=dict)
    unplaced: List[str] = field(default_factory=list)


def _scaffold_group(
    scaffold_name: str,
    ref_id: str,
    contig_names: List[str],
    query_seqs: Dict[str, str],
    ref_seqs: Dict[str, str],
    contig_orientations: Dict[str, bool],
    qr_ref_ranges: Dict[Tuple[str, str], Tuple[int, int]],
    ref_norm_to_orig: Optional[Dict[str, str]],
    have_ragtag: bool,
    work_dir: Path,
    threads: int,
    gap_size: int,
    assembly_context: str = "",
    composites: Optional[Dict[str, List[_Part]]] = None,
) -> _ScaffoldGroupResult:
    """Scaffold a single multi-contig group (RagTag or built-in fallback).

    This function is safe to call from multiple threads because:
    - query_seqs and ref_seqs are read-only shared dicts
    - Each group writes to its own work_dir / scaffold_name subdirectory
    - The built-in scaffolder is pure Python with no shared mutable state

    Args:
        scaffold_name: Output scaffold name (e.g., 'chr1A').
        ref_id: Normalized reference chromosome ID.
        contig_names: Contigs to scaffold for this group.
        query_seqs: Read-only dict of all query sequences.
        ref_seqs: Read-only dict of all reference sequences.
        contig_orientations: Dict mapping contig -> True if reversed.
        qr_ref_ranges: Dict mapping (contig, ref_id) -> (ref_min, ref_max).
        ref_norm_to_orig: Optional mapping from normalized to original ref IDs.
        have_ragtag: Whether ragtag.py is available.
        work_dir: Base working directory (per-group subdir is created).
        threads: Number of threads for RagTag.
        gap_size: Number of Ns between scaffolded contigs.
        composites: Contigs with fragments placed into their input gaps;
            RagTag scaffolds their filled sequences.

    Returns:
        _ScaffoldGroupResult with scaffold sequence, AGP lines, and confidences.
    """
    # Set thread-local assembly context for log messages (ThreadPoolExecutor
    # threads don't inherit the parent's threading.local state)
    if assembly_context:
        from dnadis.utils.logging_config import set_assembly_context
        set_assembly_context(assembly_context)

    result = _ScaffoldGroupResult(scaffold_name=scaffold_name)
    n_contigs = len(contig_names)
    logger.info(f"{scaffold_name}: scaffolding {n_contigs} contigs")

    if have_ragtag and ref_id in ref_seqs:
        # Extract reference chromosome to temp FASTA
        chr_work = work_dir / scaffold_name
        chr_work.mkdir(parents=True, exist_ok=True)
        ref_chr_fasta = chr_work / "ref_chr.fa"
        # Find original ref ID in ref_seqs
        ref_seq_key = None
        for key in ref_seqs:
            orig = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
            if key == orig or key == ref_id:
                ref_seq_key = key
                break
        if ref_seq_key is None:
            # Try to find by iterating ref_norm_to_orig
            if ref_norm_to_orig:
                orig_id = ref_norm_to_orig.get(ref_id, ref_id)
                if orig_id in ref_seqs:
                    ref_seq_key = orig_id

        if ref_seq_key and ref_seq_key in ref_seqs:
            write_fasta({ref_seq_key: ref_seqs[ref_seq_key]}, ref_chr_fasta)
            # Remove stale .fai so pysam/samtools re-indexes
            ref_fai = Path(str(ref_chr_fasta) + ".fai")
            if ref_fai.exists():
                ref_fai.unlink()

            # Extract contigs to temp FASTA
            contigs_fasta = chr_work / "contigs.fa"
            contig_seqs_subset = {}
            for cn in contig_names:
                if cn in (composites or {}):
                    contig_seqs_subset[cn] = _render_parts(
                        cn, composites[cn], query_seqs, gap_size,
                    )[0]
                elif cn in query_seqs:
                    contig_seqs_subset[cn] = query_seqs[cn]
            write_fasta(contig_seqs_subset, contigs_fasta)
            # Remove stale .fai so pysam/samtools re-indexes
            contigs_fai = Path(str(contigs_fasta) + ".fai")
            if contigs_fai.exists():
                contigs_fai.unlink()

            # Run RagTag (we use only its AGP; see _scaffold_from_ragtag_agp)
            ragtag_dir = chr_work / "ragtag_out"
            agp_path = _run_ragtag_scaffold(
                ref_chr_fasta, contigs_fasta, ragtag_dir, threads,
            )

            if agp_path:
                scaffold_seq, agp_lines, unplaced = _scaffold_from_ragtag_agp(
                    agp_path, query_seqs, scaffold_name, gap_size, ref_seq_key,
                    composites=composites,
                )
                if scaffold_seq and agp_lines:
                    if unplaced:
                        logger.warning(
                            f"{scaffold_name}: RagTag could not place "
                            f"{len(unplaced)} contig(s): {', '.join(unplaced)} "
                            f"(emitted as separate scaffolds)"
                        )
                    result.unplaced = unplaced
                    result.scaffold_seq = scaffold_seq
                    result.agp_lines = agp_lines
                    conf_path = ragtag_dir / "ragtag.scaffold.confidence.txt"
                    result.confidences = _parse_ragtag_confidence(conf_path)
                    return result
                else:
                    logger.warning(f"{scaffold_name}: RagTag AGP had no placed components, falling back")
            else:
                logger.warning(f"{scaffold_name}: RagTag failed, falling back to built-in scaffolder")
        else:
            logger.warning(f"{scaffold_name}: reference sequence not found, falling back")

    # Fallback: built-in scaffolder
    scaffold_seq, agp_lines = _builtin_scaffold(
        contig_names=contig_names,
        sequences=query_seqs,
        orientations=contig_orientations,
        ref_ranges=qr_ref_ranges,
        ref_id=ref_id,
        scaffold_name=scaffold_name,
        gap_size=gap_size,
        composites=composites,
    )
    if scaffold_seq:
        result.scaffold_seq = scaffold_seq
        result.agp_lines = agp_lines
    return result


# ---------------------------------------------------------------------------
# Main scaffolding function
# ---------------------------------------------------------------------------

def scaffold_chromosomes(
    query_fasta: Path,
    classifications: List[ContigClassification],
    contig_orientations: Dict[str, bool],
    qr_ref_ranges: Dict[Tuple[str, str], Tuple[int, int]],
    ref_fasta: Path,
    ref_lengths: Dict[str, int],
    best_ref: Dict[str, str],
    contig_refs: Dict[str, Set[str]],
    work_dir: Path,
    threads: int = 1,
    gap_size: int = 100,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
    qr_best_chain_ident: Optional[Dict[Tuple[str, str], float]] = None,
    include_rdna: bool = False,
    no_join: Optional[Set[str]] = None,
    query_gaps: Optional[Dict[str, List[Tuple[int, int]]]] = None,
    macro_block_rows=None,
    chimeric: Optional[Set[str]] = None,
) -> Tuple[Dict[str, str], List[str], Dict[str, Tuple[float, float, float]]]:
    """Produce scaffolded chromosome pseudomolecules.

    Groups contigs by (reference chromosome, haplotype group).  In each
    group, fragments are first placed into input scaffold gaps of the
    group's chrom_assigned contigs where the reference brackets them
    (``_place_fragments_in_gaps``).  Then the group:
    - Emits each T2T contig as its own scaffold, with any fragments placed
      into its gaps; the group's other fragments are not joined to it
    - Runs RagTag scaffold (if installed) on the rest
    - Falls back to built-in ordering by reference position

    Chromosome contigs not placed in any scaffold are emitted on their own,
    so the scaffold set holds every chrom_assigned/chrom_fragment contig.

    Args:
        query_fasta: Path to query assembly FASTA.
        classifications: List of ContigClassification objects.
        contig_orientations: Dict mapping contig -> True if reversed.
        qr_ref_ranges: Dict mapping (contig, ref_id) -> (ref_min, ref_max).
        ref_fasta: Path to reference FASTA.
        ref_lengths: Dict mapping ref_id -> length.
        best_ref: Dict mapping contig -> best ref_id.
        contig_refs: Dict mapping contig -> set of ref_ids with evidence.
        work_dir: Working directory for intermediate files.
        threads: Number of threads for RagTag.
        gap_size: Number of Ns between scaffolded contigs.
        ref_norm_to_orig: Optional mapping from normalized to original ref IDs.
        qr_best_chain_ident: Optional dict of (contig, ref_id) -> identity for
            haplotype-aware debris assignment.
        include_rdna: Also place rDNA contigs that have an assigned_ref_id
            (``--scaffold-rdna-contigs``).
        no_join: Contigs never joined to others; each is emitted as its own
            scaffold.  Used for sister-placed contigs, whose assignment says
            which chromosome they belong to but not where on it.  Fragments
            that copy own-chromosome sequence (``_own_assembly_copies``) are
            emitted on their own too, with ``scaffold_exclusion`` set.
        query_gaps: Contig -> 0-based half-open N runs of the input
            (``scan_fasta_gaps``).  Without it, no fragment is placed into
            a gap.
        macro_block_rows: Macro block tuples from chain parsing, which give
            the reference positions of gap flanks and fragments.
        chimeric: Contigs flagged ``chimeric`` in contig_summary; they are
            not placed into gaps.

    Returns:
        Tuple of (scaffolded_sequences, agp_lines, scaffold_confidences) where:
        - scaffolded_sequences: Dict mapping scaffold name -> DNA sequence
        - agp_lines: List of AGP lines for all scaffolds
        - scaffold_confidences: Dict mapping contig original_name ->
          (grouping, location, orientation) confidence from RagTag
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Group contigs by (ref_id, haplotype group)
    groups = _group_contigs_by_haplotype(
        classifications, best_ref, contig_refs, qr_best_chain_ident or {},
        include_rdna=include_rdna,
    )
    no_join = no_join or set()
    own_copies = _own_assembly_copies(classifications, best_ref)
    for clf in classifications:
        if clf.original_name in own_copies:
            clf.scaffold_exclusion = "own_assembly_copy"
    if any(
        c.classification == "chrom_fragment" and c.redundancy_cov is not None
        for c in classifications
    ) and not any(c.depth_ratio is not None for c in classifications):
        logger.info(
            "No read depth: only near-identical copies of chromosome sequence "
            "are kept out of scaffold joins (haplotigs need --reads)"
        )
    excluded = no_join | set(own_copies)
    if excluded:
        groups = {
            key: [c for c in members if c not in excluded]
            for key, members in groups.items()
        }
    all_confidences: Dict[str, Tuple[float, float, float]] = {}

    if not groups:
        logger.warning("No contigs to scaffold")
        return {}, [], {}

    # Capture assembly context for thread-local propagation to worker threads
    from dnadis.utils.logging_config import _assembly_context
    _asm_ctx = getattr(_assembly_context, "name", "")

    # Build classification lookup
    clf_lookup = {clf.original_name: clf for clf in classifications}

    # Check RagTag availability once
    have_ragtag = have_exe("ragtag.py")
    if not have_ragtag:
        logger.warning(
            "ragtag.py not found in PATH -- using built-in scaffolder. "
            "Install RagTag for better results: conda install -c bioconda ragtag"
        )

    # Load query sequences
    logger.info("Loading query sequences for scaffolding")
    query_seqs = read_fasta_sequences(query_fasta)

    # Load reference sequences (for RagTag)
    ref_seqs: Dict[str, str] = {}
    if have_ragtag:
        ref_seqs = read_fasta_sequences(ref_fasta)

    # Fragments placed into input scaffold gaps of the group's chromosome
    # contigs leave the group; the contigs become composites.  Fragments of a
    # group with gapped contigs that were not placed keep why in ``gap_unplaced``.
    composites: Dict[str, List[_Part]] = {}
    gap_unplaced: Dict[str, str] = {}
    if query_gaps:
        ref_blocks = _index_ref_blocks(macro_block_rows)
        query_lengths = {name: len(seq) for name, seq in query_seqs.items()}
        for (ref_id, grp), members in sorted(groups.items()):
            anchors = [
                c for c in members
                if clf_lookup.get(c) and clf_lookup[c].classification == "chrom_assigned"
            ]
            gapped = [c for c in anchors if query_gaps.get(c)]
            fragments = [c for c in members if c not in anchors and c in query_lengths]
            if not gapped or not fragments:
                continue
            fills, placed, unplaced = _place_fragments_in_gaps(
                gapped, fragments, ref_id, query_gaps, ref_blocks, query_lengths,
                chimeric=chimeric,
            )
            composites.update(fills)
            gap_unplaced.update(unplaced)
            for frag, where in placed.items():
                logger.info(f"Placed {frag} into {where}")
            if placed:
                groups[(ref_id, grp)] = [c for c in members if c not in placed]
        n_placed = sum(
            1 for parts in composites.values() for p in parts
            if p[0] == "W" and p[1] not in composites
        )
        if n_placed:
            logger.done(
                f"Placed {n_placed} fragment(s) into input scaffold gaps of "
                f"{len(composites)} chromosome contig(s)"
            )

    # Process each (ref_id, grp) — trivial scaffolds first, then parallel
    all_scaffolded: Dict[str, str] = {}
    all_agp: List[str] = []

    # ---- First pass: collect trivial scaffolds (T2T, single-contig) ----
    # These are instant and don't need parallelization.
    # Also collect multi-contig groups for the parallel second pass.
    work_list: List[Tuple[str, str, List[str]]] = []  # (scaffold_name, ref_id, contig_names)

    # References whose contigs were split into query subgenome sets; their
    # scaffolds carry the set suffix (_s1, _s2, …) like the contig names.
    split_refs = {
        clf.assigned_ref_id for clf in clf_lookup.values()
        if clf.query_subgenome and clf.assigned_ref_id
    }

    for (ref_id, grp) in sorted(groups.keys()):
        contig_names = groups[(ref_id, grp)]
        if not contig_names:
            continue

        # Determine scaffold name
        base_name = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
        if ref_id in split_refs or grp > 1:
            scaffold_name = f"{base_name}_s{grp}"
        else:
            scaffold_name = base_name

        # Telomere-to-telomere handling: emit each contig with telomeres at
        # both ends as its own scaffold, with the fragments placed into its
        # input gaps above.  The group's other fragments are not joined to
        # it; they are emitted unplaced below, since RagTag would only append
        # them at the ends.  Fragments redundant with the T2T contig were
        # already demoted to chrom_debris, so the ones left here carry
        # sequence it lacks.
        chrom_assigned_in_group = [
            c for c in contig_names
            if clf_lookup.get(c) and clf_lookup[c].classification == "chrom_assigned"
        ]
        t2t_contigs = [
            c for c in chrom_assigned_in_group
            if clf_lookup[c].is_chromosome_scale
            and clf_lookup[c].has_5p_telomere
            and clf_lookup[c].has_3p_telomere
        ]

        if t2t_contigs:
            t2t_set = set(t2t_contigs)

            for t2t_contig in t2t_contigs:
                t2t_name = clf_lookup[t2t_contig].new_name or scaffold_name
                n_gaps = clf_lookup[t2t_contig].n_gaps
                _emit_trivial_scaffold(
                    t2t_name, t2t_contig, query_seqs, contig_orientations,
                    all_scaffolded, all_agp,
                    log_label=f"T2T with {n_gaps} scaffold gap(s)" if n_gaps else "T2T",
                    composites=composites, gap_size=gap_size,
                )

            # Emit non-T2T chrom_assigned contigs as their own trivial
            # scaffolds — these are separate chromosomal copies that should
            # not be dropped just because a T2T copy exists in the same group.
            non_t2t_chrom = [
                c for c in chrom_assigned_in_group if c not in t2t_set
            ]
            dropped = [c for c in contig_names if c not in chrom_assigned_in_group]
            for contig in dropped:
                if contig in gap_unplaced:
                    clf_lookup[contig].scaffold_exclusion = gap_unplaced[contig]
            if dropped:
                n_chim = sum(1 for c in dropped if gap_unplaced.get(c) == "chimeric_fragment")
                n_gap = sum(1 for c in dropped if c in gap_unplaced) - n_chim
                reason = "T2T copy exists"
                if n_gap:
                    reason += f"; {n_gap} fit no single input gap"
                if n_chim:
                    reason += f"; {n_chim} chimeric"
                logger.info(
                    f"{scaffold_name}: leaving {len(dropped)} fragment contig(s) "
                    f"unplaced ({reason})"
                )
            for contig in non_t2t_chrom:
                clf = clf_lookup.get(contig)
                ctg_name = (clf.new_name if clf and clf.new_name else scaffold_name)
                _emit_trivial_scaffold(
                    ctg_name, contig, query_seqs, contig_orientations,
                    all_scaffolded, all_agp, log_label="non-T2T chrom copy",
                    composites=composites, gap_size=gap_size,
                )

            continue  # T2T + copies emitted; skip multi-contig scaffolding

        # Single-contig group: emit trivial AGP without invoking RagTag.
        # Uses clf.new_name to preserve fragment/copy suffixes (e.g. chr14A_f1)
        # so that pairwise FASTA headers match contig_summary names.
        if len(contig_names) == 1:
            contig = contig_names[0]
            clf = clf_lookup.get(contig)
            trivial_name = (clf.new_name if clf and clf.new_name else scaffold_name)
            _emit_trivial_scaffold(
                trivial_name, contig, query_seqs, contig_orientations,
                all_scaffolded, all_agp, log_label="single contig",
                composites=composites, gap_size=gap_size,
            )
            continue

        # Multi-contig group: needs scaffolding — add to work list
        work_list.append((scaffold_name, ref_id, contig_names))

    def _collect_group_result(result: _ScaffoldGroupResult) -> None:
        if not result.scaffold_seq:
            return
        all_scaffolded[result.scaffold_name] = result.scaffold_seq
        all_agp.extend(result.agp_lines)
        all_confidences.update(result.confidences)
        # Contigs RagTag could not place keep their contig names, like
        # single-contig groups.
        for contig in result.unplaced:
            clf = clf_lookup.get(contig)
            name = clf.new_name if clf and clf.new_name else contig
            if name in all_scaffolded:
                name = contig
            _emit_trivial_scaffold(
                name, contig, query_seqs, contig_orientations,
                all_scaffolded, all_agp, log_label="unplaced by RagTag",
                composites=composites, gap_size=gap_size,
            )

    # ---- Second pass: scaffold multi-contig groups (possibly in parallel) ----
    if work_list:
        n_jobs = len(work_list)
        if n_jobs == 1:
            # Single job: run directly with all threads, no executor overhead
            job_threads = max(1, threads)
            scaffold_name, ref_id, contig_names = work_list[0]
            logger.info(
                f"Scaffolding 1 multi-contig group with {job_threads} threads"
            )
            result = _scaffold_group(
                scaffold_name=scaffold_name,
                ref_id=ref_id,
                contig_names=contig_names,
                query_seqs=query_seqs,
                ref_seqs=ref_seqs,
                contig_orientations=contig_orientations,
                qr_ref_ranges=qr_ref_ranges,
                ref_norm_to_orig=ref_norm_to_orig,
                have_ragtag=have_ragtag,
                work_dir=work_dir,
                threads=job_threads,
                gap_size=gap_size,
                assembly_context=_asm_ctx,
                composites=composites,
            )
            _collect_group_result(result)
        else:
            # Multiple jobs: run in parallel
            max_workers = min(n_jobs, max(1, threads // 4))
            job_threads = max(1, threads // max_workers)
            logger.info(
                f"Scaffolding {n_jobs} multi-contig groups in parallel "
                f"({max_workers} workers, {job_threads} threads each)"
            )

            def _submit_job(
                item: Tuple[str, str, List[str]],
            ) -> _ScaffoldGroupResult:
                s_name, r_id, c_names = item
                return _scaffold_group(
                    scaffold_name=s_name,
                    ref_id=r_id,
                    contig_names=c_names,
                    query_seqs=query_seqs,
                    ref_seqs=ref_seqs,
                    contig_orientations=contig_orientations,
                    qr_ref_ranges=qr_ref_ranges,
                    ref_norm_to_orig=ref_norm_to_orig,
                    have_ragtag=have_ragtag,
                    work_dir=work_dir,
                    threads=job_threads,
                    gap_size=gap_size,
                    assembly_context=_asm_ctx,
                    composites=composites,
                )

            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = [pool.submit(_submit_job, item) for item in work_list]
                # Collect results in submission order (= sorted key order)
                for future in futures:
                    _collect_group_result(future.result())

    # Every chromosome contig belongs in the scaffold set.  Contigs not placed
    # above (fragments beside a T2T copy that fit no single input gap, extra
    # copies, fragments that copy own-chromosome sequence, fragments without
    # a reference, no_join contigs, groups that produced no scaffold) are
    # emitted on their own under their contig names, like contigs RagTag
    # could not place.
    placed = {
        fields[5] for fields in (line.split("\t") for line in all_agp)
        if len(fields) > 5 and fields[4] == "W"
    }
    for clf in classifications:
        if clf.classification not in ("chrom_assigned", "chrom_fragment"):
            continue
        if clf.original_name in placed:
            continue
        name = clf.new_name or clf.original_name
        if name in all_scaffolded:
            name = clf.original_name
        _emit_trivial_scaffold(
            name, clf.original_name, query_seqs, contig_orientations,
            all_scaffolded, all_agp,
            log_label=(
                "sister-placed, not joined" if clf.original_name in no_join
                else f"copy of {own_copies[clf.original_name]}, not joined:"
                if clf.original_name in own_copies
                else "unplaced"
            ),
            composites=composites, gap_size=gap_size,
        )

    logger.done(f"Scaffolded {len(all_scaffolded)} chromosomes")
    return all_scaffolded, all_agp, all_confidences
