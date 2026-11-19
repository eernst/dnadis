#!/usr/bin/env python3
"""
Organelle (chrC/chrM) detection for dnadis.

Contains functions for identifying chloroplast and mitochondrial contigs
using BLAST-based alignment to reference organelle sequences.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

from dnadis.detection.blast import run_makeblastdb, run_blastn_megablast
from dnadis.models import OrganelleCalls, OrganelleHit, OrganelleUnit
from dnadis.utils.io_utils import merge_intervals, open_maybe_gzip, replace_if_changed
from dnadis.utils.logging_config import get_logger
from dnadis.utils.reference_utils import (
    bare_organelle_names_allowed,
    find_organelle_sequences,
    organelle_ref_name,
    parse_organelle_id,
    reference_species_tags,
)
from dnadis.utils.sequence_utils import (
    write_filtered_fasta,
    is_hifiasm_circular,
    read_fasta_sequences,
    write_fasta,
)

logger = get_logger("organelle")


def read_fasta_headers(fasta: Path) -> List[Tuple[str, str, int]]:
    """(id, description, length) of every sequence, streaming the file."""
    headers: List[Tuple[str, str, int]] = []
    seq_id, desc, length = None, "", 0
    with open_maybe_gzip(fasta, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    headers.append((seq_id, desc, length))
                parts = line[1:].strip().split(None, 1)
                seq_id = parts[0] if parts else ""
                desc = parts[1] if len(parts) > 1 else ""
                length = 0
            else:
                length += len(line.strip())
    if seq_id is not None:
        headers.append((seq_id, desc, length))
    return headers


def _read_selected_sequences(fasta: Path, names: Set[str]) -> Dict[str, str]:
    """Sequences of *names* from *fasta*, streaming the file."""
    seqs: Dict[str, List[str]] = {}
    current = None
    with open_maybe_gzip(fasta, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].strip().split(None, 1)[0] if line[1:].strip() else ""
                current = current if current in names else None
                if current is not None:
                    seqs[current] = []
            elif current is not None:
                seqs[current].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _write_organelle_refs(
    chosen: Dict[Tuple[str, Optional[str]], Tuple[Path, str]],
    organelle: str,
    out_path: Path,
) -> Optional[Path]:
    """Write every reference of *organelle* to *out_path*, named chrC/chrC_A...

    *chosen* maps (organelle, tag) to (source FASTA, sequence id).  The file
    is replaced only when its content changes.
    """
    keys = sorted((k for k in chosen if k[0] == organelle), key=lambda k: (k[1] or ""))
    if not keys:
        return None
    by_file: Dict[Path, Set[str]] = defaultdict(set)
    for k in keys:
        src, seq_id = chosen[k]
        by_file[src].add(seq_id)
    sequences = {src: _read_selected_sequences(src, ids) for src, ids in by_file.items()}
    out: Dict[str, str] = {}
    for k in keys:
        src, seq_id = chosen[k]
        name = organelle_ref_name(*k)
        out[name] = sequences[src][seq_id]
        logger.info(
            f"{name} reference: {seq_id} from {src.name} ({len(out[name]):,} bp)"
        )
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    write_fasta(out, tmp)
    replace_if_changed(tmp, out_path)
    return out_path


def prepare_organelle_references(
    ref_fasta: Path,
    chrC_ref_arg: Optional[str],
    chrM_ref_arg: Optional[str],
    work_dir: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Prepare organelle reference FASTAs, one per organelle.

    Every organelle sequence recognised in *ref_fasta* is used
    (``reference_utils.find_organelle_sequences``): by name, with a species
    tag matching the nuclear subgenome suffixes (chrC_A, PtP...), or by an
    NCBI-style description.  ``--chrC-ref``/``--chrM-ref`` replace the
    reference's sequences for that organelle; every sequence in such a file
    is used, tagged when its name carries a species tag.  Sequences are
    named chrC, chrC_A, ... in the written files.

    Returns (chrC_fasta_path, chrM_fasta_path); None for an organelle with
    no reference.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    ref_headers = read_fasta_headers(ref_fasta)
    ref_ids = [h[0] for h in ref_headers]
    tags = reference_species_tags(ref_ids)
    allow_bare = bare_organelle_names_allowed(ref_ids)

    chosen: Dict[Tuple[str, Optional[str]], Tuple[Path, str]] = {
        key: (ref_fasta, seq_id)
        for key, seq_id in find_organelle_sequences(ref_headers, tags, allow_bare).items()
    }
    for organelle, arg in (("chrC", chrC_ref_arg), ("chrM", chrM_ref_arg)):
        if not arg:
            continue
        path = Path(arg)
        if not path.exists():
            logger.warning(f"{organelle} reference not found: {arg}")
            continue
        chosen = {k: v for k, v in chosen.items() if k[0] != organelle}
        user: Dict[Optional[str], List[Tuple[int, str]]] = defaultdict(list)
        for seq_id, _, length in read_fasta_headers(path):
            parsed = parse_organelle_id(seq_id, tags, allow_bare)
            tag = parsed[1] if parsed and parsed[0] == organelle else None
            user[tag].append((length, seq_id))
        for tag, seqs in user.items():
            seqs.sort(key=lambda x: (-x[0], x[1]))
            chosen[(organelle, tag)] = (path, seqs[0][1])
            if len(seqs) > 1:
                logger.warning(
                    f"{len(seqs)} sequences in {path.name} match {organelle_ref_name(organelle, tag)}: "
                    f"using the longest, {seqs[0][1]}; ignoring " + ", ".join(n for _, n in seqs[1:])
                )
        logger.info(f"Using user-provided {organelle} reference: {path}")

    paths = []
    for organelle in ("chrC", "chrM"):
        p = _write_organelle_refs(chosen, organelle, work_dir / f"{organelle}_ref.fa")
        if p is None:
            logger.info(f"No {organelle} reference found")
        paths.append(p)
    return paths[0], paths[1]


@dataclass(frozen=True)
class OrganelleThresholds:
    """Criteria for organelle calls (see ``call_organelles``).

    A contig is ``complete`` when it covers ``min_ref_cov`` of the reference
    (each contig position counted once for the ``one_to_one`` organelles),
    at least ``min_query_cov`` of the contig aligns, its length is within
    the organelle's (min, max) ratio to the reference length, and, for the
    ``closed_circle`` organelles, a linear contig's ends meet on the
    reference.  When no contig
    passes, the fewest contigs (each at ``min_query_cov``) that together cover
    ``set_min_ref_cov`` of the reference are ``partial``.  Other contigs with
    at least ``debris_min_query_cov`` aligned are ``debris``.
    """
    min_ref_cov: float = 0.95
    min_query_cov: float = 0.90
    set_min_ref_cov: float = 0.90
    debris_min_query_cov: float = 0.50
    same_genome_identity: float = 0.99
    # Organelles whose oversized contigs (tandem concatemers) may be cut to
    # one unit.  Not chrM: a long mitochondrial contig is a walk through the
    # repeat-mediated recombination graph, not a tandem array.
    cut_units: Tuple[str, ...] = ("chrC",)
    # Organelles whose reference coverage counts each contig position once
    # (one_to_one_ref_intervals()).  A plastome's inverted repeats match both
    # reference copies, so a contig with one IR copy partly missing still
    # covers the whole reference by the union of HSPs.  Not chrM: assemblies
    # collapse mitochondrial repeat copies the reference carries, and only
    # their cross-matches cover the reference.
    one_to_one: Tuple[str, ...] = ("chrC",)
    # Organelles whose linear contigs are complete only when their two ends
    # meet on the circular reference, at most max_end_gap bp apart
    # (end_join()); contigs whose ends overlap are trimmed to one copy.  A
    # plastome assembled linear can stop short of closing and still cover
    # nearly all of the reference.  Not chrM (repeat-graph walks).
    closed_circle: Tuple[str, ...] = ("chrC",)
    max_end_gap: int = 1000
    len_ratio: Tuple[Tuple[str, float, float], ...] = (
        ("chrC", 0.90, 1.50),
        ("chrM", 0.50, 3.00),
    )

    def len_bounds(self, organelle: str) -> Tuple[float, float]:
        for name, lo, hi in self.len_ratio:
            if name == organelle:
                return lo, hi
        return 0.0, float("inf")


# Complete candidates whose length differs from the reference by the same
# number of these steps (fraction of the reference length) are ranked by
# their end overlap: one-copy units trimmed from contigs of one genome differ
# by a few bp in the estimated period.
_LENGTH_STEP = 0.005

# A set member must add at least this fraction of the reference.
_SET_MIN_GAIN = 0.01

# Two contigs belong to the same organelle genome at or above
# OrganelleThresholds.same_genome_identity with at least this fraction of the
# query aligned.  Copies of one plastome in an assembly are near-identical;
# plastomes of related species are typically 94-99% identical.
_SAME_GENOME_MIN_COV = 0.50

# Without read depth, the genome with the most organelle sequence is used only
# when it has at least this many times the sequence of the next one.
_GROUP_DOMINANCE = 2.0

# (organelle, complete candidates, all contigs assigned to that organelle)
#   -> (candidate -> genome group, group -> bp of assigned contigs in it)
GenomeGrouper = Callable[[str, List[str], List[str]], Tuple[Dict[str, int], Dict[int, int]]]


@dataclass
class _PairEvidence:
    query_intervals: List[Tuple[int, int]] = field(default_factory=list)
    ref_intervals: List[Tuple[int, int]] = field(default_factory=list)
    matches: float = 0.0
    aln_len: int = 0
    # Per HSP: query start, query end, identity, alignment length, reference
    # start, reference end (half-open), and whether it aligns to the minus strand
    hsps: List[Tuple[int, int, float, int, int, int, bool]] = field(default_factory=list)


def parse_organelle_blast(blast_out: Path) -> Dict[Tuple[str, str], _PairEvidence]:
    """Collect query/reference intervals and matches per (contig, organelle)."""
    pairs: Dict[Tuple[str, str], _PairEvidence] = defaultdict(_PairEvidence)
    if not (blast_out.exists() and blast_out.stat().st_size > 0):
        return pairs
    with blast_out.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                pident = float(fields[2])
                aln_length = int(fields[3])
                q1, q2, s1, s2 = (int(fields[i]) for i in (6, 7, 8, 9))
                qstart, qend = sorted((q1, q2))
                sstart, send = sorted((s1, s2))
                minus = (s1 > s2) != (q1 > q2)
            except ValueError:
                continue
            ev = pairs[(fields[0], fields[1])]
            # BLAST coordinates are 1-based inclusive; store half-open intervals
            ev.query_intervals.append((qstart - 1, qend))
            ev.ref_intervals.append((sstart - 1, send))
            ev.matches += pident * aln_length / 100.0
            ev.aln_len += aln_length
            ev.hsps.append((qstart - 1, qend, pident / 100.0, aln_length,
                            sstart - 1, send, minus))
    return pairs


def organelle_of_ref(ref_name: str) -> str:
    """chrC/chrM of an organelle reference name (chrC, chrC_A, ...)."""
    return ref_name.split("_", 1)[0]


def _covered_bp(intervals: List[Tuple[int, int]]) -> int:
    return merge_intervals(intervals)[1] if intervals else 0


def _subtract(iv: Tuple[int, int], taken: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Parts of half-open *iv* outside every interval in *taken*."""
    parts = [iv]
    for s, e in taken:
        parts = [p for a, b in parts
                 for p in ((a, min(b, s)), (max(a, e), b)) if p[0] < p[1]]
    return parts


def one_to_one_ref_intervals(
    hsps: List[Tuple[int, int, float, int, int, int, bool]],
) -> List[Tuple[int, int]]:
    """Merged reference intervals with each query and reference position used once.

    HSPs are taken longest first.  Each keeps only the part of its query
    interval not yet taken, mapped to the reference (linearly, ignoring
    indels), and of that only the reference not yet taken.  A repeat copy in
    the contig thus covers one reference copy: the collinear HSP when it is
    longer, else whichever cross-match comes first, leaving the other
    reference copy to the contig's other repeat copy.
    """
    q_taken: List[Tuple[int, int]] = []
    r_taken: List[Tuple[int, int]] = []
    for qs, qe, _, length, ss, se, minus in sorted(hsps, key=lambda h: (-h[3], -h[2], h[0], h[4])):
        if qe <= qs or se <= ss:
            continue
        scale = (se - ss) / (qe - qs)

        def to_ref(a: int, b: int) -> Tuple[int, int]:
            if minus:
                return se - round((b - qs) * scale), se - round((a - qs) * scale)
            return ss + round((a - qs) * scale), ss + round((b - qs) * scale)

        def to_query(x: int, y: int) -> Tuple[int, int]:
            if minus:
                return qs + round((se - y) / scale), qs + round((se - x) / scale)
            return qs + round((x - ss) / scale), qs + round((y - ss) / scale)

        for a, b in _subtract((qs, qe), q_taken):
            for x, y in _subtract(to_ref(a, b), r_taken):
                r_taken.append((x, y))
                q_taken.append(to_query(x, y))
    return merge_intervals(r_taken)[0]


# HSPs starting (ending) within this many bp of the contig's first (last)
# aligned base are tested as its ends: a contig end inside an inverted repeat
# aligns to both reference copies.
_END_SLACK = 1000


def end_join(
    hsps: List[Tuple[int, int, float, int, int, int, bool]],
    qlen: int,
    ref_len: int,
) -> Optional[Tuple[int, Optional[int]]]:
    """How a linear contig's two ends meet on the circular reference.

    Returns ``(gap, period)``.  ``gap`` is the reference sequence between the
    reference position at the contig's end and the one at its start (read
    in the contig's direction), less the contig sequence outside those end
    HSPs, which may fill it (down to 0); negative when the ends overlap.
    ``period`` is the contig length of one copy of the circle when the ends
    overlap: the contig position, after its start, that aligns to the
    reference position at its start (None if no HSP holds it).  Returns
    None when no HSPs at the two ends align to the same strand.
    """
    if not hsps or qlen <= 0 or ref_len <= 0:
        return None
    first = min(h[0] for h in hsps)
    last = max(h[1] for h in hsps)
    starts = [h for h in hsps if h[0] <= first + _END_SLACK]
    ends = [h for h in hsps if h[1] >= last - _END_SLACK]
    best: Optional[Tuple[int, int, tuple, tuple]] = None
    for a in starts:
        for b in ends:
            if a[6] != b[6]:
                continue
            # Plus strand: the contig continues past b's reference end to a's
            # reference start.  Minus: past b's reference start (reading
            # down) to a's reference end.  The contig sequence before a and
            # after b may fill a gap, but does not make an overlap.
            d = (a[4] - b[5]) if not a[6] else (b[4] - a[5])
            d = (d + ref_len // 2) % ref_len - ref_len // 2
            gap = d if d <= 0 else max(0, d - a[0] - (qlen - b[1]))
            if best is None or abs(gap) < abs(best[0]):
                best = (gap, d, a, b)
    if best is None:
        return None
    gap, d, a, b = best
    period = None
    if gap < 0:
        # Of the contig positions aligned to a's first reference base, the
        # one closest to where the overlap puts the start of the next copy
        expected = b[1] + d - a[0]
        qs_a, _, _, _, ss_a, se_a, minus = a
        x = se_a - 1 if minus else ss_a
        for qs, qe, _, _, ss, se, m in hsps:
            if m != minus or not ss <= x < se or qe <= qs:
                continue
            scale = (qe - qs) / (se - ss)
            p = qs + round(((se - 1 - x) if minus else (x - ss)) * scale) - qs_a
            if p > _END_SLACK and (period is None or abs(p - expected) < abs(period - expected)):
                period = p
    return gap, period


def _window_stats(
    hsps: List[Tuple[int, int, float, int, int, int, bool]],
    start: int,
    end: int,
    ref_len: int,
    one_to_one: bool,
) -> Tuple[float, float, float, float]:
    """Reference coverage, aligned fraction, plus-strand aligned fraction and
    identity of the contig window [start, end), from HSPs clipped to it."""
    q_iv, plus_iv, clipped = [], [], []
    matches = aligned = 0.0
    for qs, qe, ident, length, ss, se, minus in hsps:
        a, b = max(qs, start), min(qe, end)
        if a >= b:
            continue
        q_iv.append((a, b))
        if minus:
            r = (max(ss, se - (b - qs)), max(ss, se - (a - qs)))
        else:
            r = (min(se, ss + (a - qs)), min(se, ss + (b - qs)))
            plus_iv.append((a, b))
        clipped.append((a, b, ident, b - a, *r, minus))
        matches += ident * (b - a)
        aligned += b - a
    if one_to_one:
        ref_cov = _covered_bp(one_to_one_ref_intervals(clipped)) / ref_len
    else:
        ref_cov = _covered_bp([h[4:6] for h in clipped if h[4] < h[5]]) / ref_len
    span = end - start
    return (ref_cov, _covered_bp(q_iv) / span, _covered_bp(plus_iv) / span,
            matches / aligned if aligned else 0.0)


def call_organelles(
    pairs: Dict[Tuple[str, str], _PairEvidence],
    query_lengths: Dict[str, int],
    ref_lengths: Dict[str, int],
    thresholds: OrganelleThresholds = OrganelleThresholds(),
    circular_contigs: Optional[Set[str]] = None,
    genome_groups: Optional[GenomeGrouper] = None,
) -> OrganelleCalls:
    """Call complete, partial and debris organelle contigs from BLAST evidence.

    Each contig is assigned to the organelle it matches best (aligned fraction
    of the contig x identity): plastid-derived stretches of the mitogenome
    otherwise label chloroplast fragments as mitochondrial.

    A linear contig of a ``closed_circle`` organelle (chrC) is complete only
    when its ends meet on the reference (``end_join()``); one whose ends
    overlap is trimmed to one copy (``hit.unit``).  Complete candidates are
    ranked circular first, then by length (of the unit, if any) closest to
    the reference, then by end overlap.  A candidate that is not selected
    stays organelle debris.

    A sample can hold the organelle genome of another species (a second
    plant in the culture, contamination), and its contigs can pass the
    complete tests.  With ``genome_groups``, candidates are grouped into
    genomes by their identity to each other.  A genome with at least twice the
    assembled sequence of the next is ranked alone (the others are reported
    in ``OrganelleCalls.other_genomes``); otherwise all candidates are ranked
    and ``repick_by_read_depth()`` chooses by read depth after depth analysis.
    """
    circular = circular_contigs if circular_contigs is not None else {
        c for c, _ in pairs if is_hifiasm_circular(c)
    }

    # Best organelle reference (organelle and species) per contig.  With
    # several species' references for an organelle, the margin is the
    # contig's identity to its best one minus that to the next.
    # contig -> [(rank key, hit, evidence, covered reference intervals)]
    by_contig: Dict[str, List[tuple]] = defaultdict(list)
    for (contig, ref_name), ev in pairs.items():
        ref_len = ref_lengths.get(ref_name, 0)
        qlen = query_lengths.get(contig, 0)
        if ref_len <= 0 or qlen <= 0:
            continue
        ref_ivs = (one_to_one_ref_intervals(ev.hsps)
                   if organelle_of_ref(ref_name) in thresholds.one_to_one
                   else merge_intervals(ev.ref_intervals)[0])
        hit = OrganelleHit(
            organelle_type=organelle_of_ref(ref_name),
            coverage=_covered_bp(ev.query_intervals) / qlen,
            identity=(ev.matches / ev.aln_len) if ev.aln_len else 0.0,
            length_ratio=qlen / ref_len,
            ref_coverage=_covered_bp(ref_ivs) / ref_len,
            ref_name=ref_name,
        )
        by_contig[contig].append(((hit.coverage * hit.identity, hit.ref_coverage), hit, ev, ref_ivs))
    best: Dict[str, Tuple[str, OrganelleHit, _PairEvidence]] = {}
    covered_ref: Dict[str, List[Tuple[int, int]]] = {}
    for contig, options in by_contig.items():
        # The organelle is the one matched best (aligned fraction x identity);
        # within it, the species reference with the highest identity over
        # non-overlapping best HSPs.  Pooled identity counts inverted-repeat
        # and repeat cross-matches, and small differences in aligned fraction
        # would outweigh the 0.5-3 point identity differences between species.
        organelle = max(options, key=lambda o: (o[0], o[1].ref_name))[1].organelle_type
        same = [(_best_hsp_identity(ev)[0], hit.coverage, hit.ref_name, hit, ev, ivs)
                for _, hit, ev, ivs in options if hit.organelle_type == organelle]
        same.sort(key=lambda o: o[:3], reverse=True)
        _, _, _, hit, ev, ivs = same[0]
        if len(same) > 1:
            hit.ref_margin = same[0][0] - same[1][0]
        best[contig] = (organelle, hit, ev)
        covered_ref[contig] = ivs

    def ref_len_of(contig: str) -> int:
        return ref_lengths[best[contig][1].ref_name]

    calls = OrganelleCalls()
    for organelle in sorted({organelle_of_ref(r) for r in ref_lengths}):
        members = {c: (h, ev) for c, (o, h, ev) in best.items() if o == organelle}
        lo, hi = thresholds.len_bounds(organelle)

        # Linear contigs whose ends do not meet on the reference; contigs
        # whose ends overlap get the one-copy unit
        open_ends: Dict[str, Optional[int]] = {}
        overlap: Dict[str, int] = {}
        if organelle in thresholds.closed_circle:
            for c, (h, ev) in members.items():
                if c in circular:
                    continue
                join = end_join(ev.hsps, query_lengths[c], ref_len_of(c))
                if join is None:
                    open_ends[c] = None
                    continue
                gap, period = join
                overlap[c] = max(0, -gap)
                fits = (period is not None and period < query_lengths[c]
                        and lo <= period / ref_len_of(c) <= hi)
                # An overlap longer than the allowed gap must hold one copy:
                # sequence at an end that aligns elsewhere (a divergent or
                # chimeric piece) can mimic one
                if gap > thresholds.max_end_gap or (gap < -thresholds.max_end_gap and not fits):
                    open_ends[c] = gap
                elif fits and h.length_ratio <= hi:
                    ref_cov, query_cov, _, identity = _window_stats(
                        ev.hsps, 0, period, ref_len_of(c), organelle in thresholds.one_to_one)
                    h.unit = OrganelleUnit(0, period, False, ref_cov, query_cov, identity)

        def passes(c: str) -> bool:
            h = members[c][0]
            ref_cov, query_cov, length = _effective(h, query_lengths[c])
            return (ref_cov >= thresholds.min_ref_cov
                    and query_cov >= thresholds.min_query_cov
                    and lo <= length / ref_len_of(c) <= hi)

        candidates = [c for c in members if passes(c) and c not in open_ends]
        for c in sorted(set(open_ends) & {c for c in members if passes(c)}):
            gap = open_ends[c]
            logger.info(
                f"{organelle}: {c} is not complete: " + (
                    f"its ends align to opposite strands of {members[c][0].ref_name}"
                    if gap is None else
                    f"its ends are {gap:,} bp apart on {members[c][0].ref_name} "
                    f"(> {thresholds.max_end_gap:,})" if gap > 0 else
                    f"its ends overlap by {-gap:,} bp on {members[c][0].ref_name}, "
                    f"but not by one copy")
            )

        def length_off(c: str) -> int:
            return abs(_effective(members[c][0], query_lengths[c])[2] - ref_len_of(c))

        # Circular first, then length closest to the reference (in steps, so
        # that trimmed copies of one genome tie), then the least sequence
        # repeated at the ends
        candidates.sort(key=lambda c: (
            c not in circular,
            round(length_off(c) / (_LENGTH_STEP * ref_len_of(c))),
            overlap.get(c, 0),
            length_off(c),
            -_effective(members[c][0], query_lengths[c])[0],
            query_lengths[c],
            c,
        ))
        all_candidates = set(candidates)
        assigned = sorted(c for c, (h, _) in members.items()
                          if h.coverage >= thresholds.debris_min_query_cov
                          or c in all_candidates)

        if genome_groups is not None and len(candidates) > 1:
            group_of, group_bp = genome_groups(organelle, candidates, assigned)
            groups = {group_of[c] for c in candidates}
            if len(groups) > 1:
                calls.genome_groups[organelle] = dict(group_of)
                calls.complete_candidates[organelle] = list(candidates)
                first_rank = {g: min(i for i, c in enumerate(candidates) if group_of[c] == g)
                              for g in groups}
                by_bp = sorted(groups, key=lambda g: (-group_bp.get(g, 0), first_rank[g]))
                main, runner_up = by_bp[0], by_bp[1]
                if group_bp.get(main, 0) >= _GROUP_DOMINANCE * group_bp.get(runner_up, 0):
                    others = [c for c in candidates if group_of[c] != main]
                    candidates = [c for c in candidates if group_of[c] == main]
                    calls.other_genomes[organelle] = others
                    logger.warning(
                        f"{organelle}: candidates come from {len(groups)} different genomes "
                        f"(< {thresholds.same_genome_identity:.1%} identity to each other). Using "
                        f"the one with the most {organelle} sequence in the assembly "
                        f"({group_bp.get(main, 0):,} vs {group_bp.get(runner_up, 0):,} bp); "
                        f"kept as organelle_debris: " + ", ".join(others)
                    )
                else:
                    logger.warning(
                        f"{organelle}: candidates come from {len(groups)} different genomes "
                        f"(< {thresholds.same_genome_identity:.1%} identity to each other) with "
                        f"similar amounts of sequence ({group_bp.get(main, 0):,} vs "
                        f"{group_bp.get(runner_up, 0):,} bp); ranking all of them. With read "
                        f"depth, the genome with the most read bases is used."
                    )

        chosen: Optional[str] = candidates[0] if candidates else None
        if chosen is None and organelle in thresholds.cut_units:
            chosen = _cut_concatemer(
                members, ref_lengths, lo, hi, thresholds, circular, query_lengths, organelle,
                genome_groups, assigned)
        partial: List[str] = []
        if chosen is None:
            # Reference coordinates differ between species' references, so a
            # set is built within one reference; the one covering most wins.
            best_set: Tuple[List[str], float, str] = ([], 0.0, "")
            for ref_name in sorted({h.ref_name for h, _ in members.values()}):
                eligible = {c: covered_ref[c] for c, (h, _) in members.items()
                            if h.ref_name == ref_name and h.coverage >= thresholds.min_query_cov}
                found, cov = _greedy_organelle_set(
                    eligible, ref_lengths[ref_name], thresholds.set_min_ref_cov,
                    circular, query_lengths,
                )
                if found and (cov, -len(found)) > (best_set[1], -len(best_set[0]) if best_set[0] else 0):
                    best_set = (found, cov, ref_name)
            partial = best_set[0]
            if partial:
                calls.partial[organelle] = partial
                calls.set_ref_coverage[organelle] = best_set[1]
                calls.ref_name[organelle] = best_set[2]
        else:
            h = members[chosen][0]
            calls.set_ref_coverage[organelle] = h.unit.ref_coverage if h.unit else h.ref_coverage
            calls.ref_name[organelle] = h.ref_name

        for contig, (hit, _) in members.items():
            if contig == chosen:
                hit.status = "complete"
            elif contig in partial:
                hit.status = "partial"
            elif hit.coverage >= thresholds.debris_min_query_cov or contig in all_candidates:
                # Candidates stay in hits for the read-depth re-pick
                calls.debris.add(contig)
            else:
                continue
            calls.hits[contig] = hit

        if organelle == "chrC":
            calls.chrC_contig = chosen
        elif organelle == "chrM":
            calls.chrM_contig = chosen

        if chosen:
            h = members[chosen][0]
            margin = f", {h.ref_margin:+.4f} identity over the next species" if h.ref_margin is not None else ""
            logger.info(
                f"Selected {organelle}: {chosen} ({query_lengths[chosen]:,} bp, "
                f"{h.length_ratio:.2f}x {h.ref_name}, ref cov={h.ref_coverage:.3f}, "
                f"query cov={h.coverage:.3f}, ident={h.identity:.4f}{margin})"
            )
            if h.unit is not None and h.unit.start == 0 and h.length_ratio <= hi:
                logger.info(
                    f"{organelle}: {chosen} ends overlap by "
                    f"{query_lengths[chosen] - h.unit.end:,} bp; writing 1-{h.unit.end:,} "
                    f"as {organelle} (ref cov={h.unit.ref_coverage:.3f})"
                )
        elif partial:
            logger.info(
                f"No single {organelle} contig passes; {len(partial)} contigs together "
                f"cover {calls.set_ref_coverage[organelle]:.3f} of the reference: "
                + ", ".join(partial)
            )
        elif members:
            top = max(members, key=lambda c: members[c][0].ref_coverage)
            h = members[top][0]
            logger.info(
                f"No {organelle} call; best contig {top} covers {h.ref_coverage:.3f} "
                f"of the reference ({h.length_ratio:.2f}x reference length, "
                f"query cov={h.coverage:.3f})"
            )
    return calls


# Reference positions used to measure the repeat period of a concatemer.
_PERIOD_ANCHORS = tuple(i / 50 for i in range(1, 50, 2))
# Fraction of a unit window that must align on the contig's main strand: a
# head-to-head (palindromic) window is half minus strand.  Below 1 because an
# inverted single-copy region (plastid flip-flop isomer) aligns minus.
_UNIT_MIN_PLUS_COV = 0.80


def _effective(hit: OrganelleHit, qlen: int) -> Tuple[float, float, int]:
    """Reference coverage, aligned fraction and length of what is written as
    the organelle: the unit when one was cut or trimmed, else the contig."""
    if hit.unit is not None:
        return hit.unit.ref_coverage, hit.unit.query_coverage, hit.unit.end - hit.unit.start
    return hit.ref_coverage, hit.coverage, qlen


def _cut_concatemer(members, ref_lengths, lo, hi, thresholds, circular, query_lengths, organelle,
                    genome_groups=None, assigned=()):
    """Cut one unit from the best oversized contig, if any passes.

    Oversized contigs cover the reference and align along their length but
    exceed the max length ratio.  Ranked circular first, then shortest.
    With ``genome_groups``, a contig is cut only if its genome holds at least
    half of the sequence assigned to the organelle (another species'
    concatemer in the sample is not cut).  Returns the contig name with
    ``hit.unit`` set, or None.
    """
    assigned_bp = sum(query_lengths.get(c, 0) for c in assigned)
    oversized = [
        c for c, (h, _) in members.items()
        if h.ref_coverage >= thresholds.min_ref_cov
        and h.coverage >= thresholds.min_query_cov
        and h.length_ratio > hi
    ]
    oversized.sort(key=lambda c: (c not in circular, query_lengths[c], c))
    for contig in oversized:
        hit, ev = members[contig]
        unit = find_concatemer_unit(
            ev.hsps, query_lengths[contig], ref_lengths[hit.ref_name], lo, hi,
            thresholds.min_ref_cov, thresholds.min_query_cov,
            one_to_one=organelle in thresholds.one_to_one,
        )
        if unit is None:
            logger.info(
                f"{organelle}: no unit cut from {contig} ({hit.length_ratio:.2f}x reference): "
                f"no consistent repeat period or no window passes"
            )
            continue
        if genome_groups is not None and assigned_bp:
            _, group_bp = genome_groups(organelle, [contig], list(assigned))
            if group_bp.get(0, 0) < 0.5 * assigned_bp:
                logger.warning(
                    f"{organelle}: not cutting {contig}: its genome holds "
                    f"{group_bp.get(0, 0):,} of {assigned_bp:,} bp assigned to {organelle}"
                )
                continue
        hit.unit = unit
        logger.info(
            f"{organelle}: cut one unit from concatemer {contig} "
            f"({query_lengths[contig]:,} bp, {hit.length_ratio:.2f}x reference): "
            f"{unit.start + 1:,}-{unit.end:,}{' (reverse)' if unit.reverse else ''}, "
            f"{unit.end - unit.start:,} bp, ref cov={unit.ref_coverage:.3f}, "
            f"ident={unit.identity:.4f}"
        )
        return contig
    return None


def find_concatemer_unit(
    hsps: List[Tuple[int, int, float, int, int, int, bool]],
    qlen: int,
    ref_len: int,
    min_len_ratio: float,
    max_len_ratio: float,
    min_ref_cov: float,
    min_query_cov: float,
    one_to_one: bool = False,
) -> Optional[OrganelleUnit]:
    """One reference-length unit of a tandem concatemer, from its HSPs.

    The period is the median distance between successive occurrences of the
    same reference positions along the contig (on its main strand).  Windows
    of one period that start at the reference start are tested like a
    complete contig: reference coverage, aligned fraction, length ratio, and
    collinearity (most of the window on the main strand).  The window with
    the highest reference coverage is returned.  With *one_to_one*, the
    window's reference coverage counts each window position once
    (``one_to_one_ref_intervals()``).
    """
    if not hsps or qlen <= 0 or ref_len <= 0:
        return None
    plus_bp = sum(h[3] for h in hsps if not h[6])
    minus_bp = sum(h[3] for h in hsps if h[6])
    reverse = minus_bp > plus_bp
    # Orient so the main strand is plus
    oriented = [
        (qlen - qe, qlen - qs, ident, length, ss, se, not minus) if reverse
        else (qs, qe, ident, length, ss, se, minus)
        for qs, qe, ident, length, ss, se, minus in hsps
    ]

    def occurrences(x: int) -> List[int]:
        return sorted({qs + (x - ss) for qs, qe, _, _, ss, se, minus in oriented
                       if not minus and ss <= x < se and qs + (x - ss) < qe})

    diffs = []
    for f in _PERIOD_ANCHORS:
        pos = occurrences(int(f * ref_len))
        diffs += [b - a for a, b in zip(pos, pos[1:])
                  if min_len_ratio * ref_len <= b - a <= max_len_ratio * ref_len]
    if len(diffs) < 2:
        return None
    diffs.sort()
    period = diffs[len(diffs) // 2]
    if sum(abs(d - period) <= 0.01 * period for d in diffs) < len(diffs) / 2:
        return None

    x0 = int(_PERIOD_ANCHORS[0] * ref_len)
    best: Optional[Tuple[Tuple[float, float], OrganelleUnit]] = None
    for occ in occurrences(x0):
        start, end = occ - x0, occ - x0 + period
        if start < 0 or end > qlen:
            continue
        ref_cov, query_cov, plus_cov, identity = _window_stats(
            oriented, start, end, ref_len, one_to_one)
        if (ref_cov < min_ref_cov or query_cov < min_query_cov
                or plus_cov < _UNIT_MIN_PLUS_COV):
            continue
        if reverse:
            start, end = qlen - end, qlen - start
        unit = OrganelleUnit(start, end, reverse, ref_cov, query_cov, identity)
        if best is None or (ref_cov, identity) > best[0]:
            best = ((ref_cov, identity), unit)
    return best[1] if best else None


def _greedy_organelle_set(
    eligible: Dict[str, List[Tuple[int, int]]],
    ref_len: int,
    min_ref_cov: float,
    circular: Set[str],
    query_lengths: Dict[str, int],
) -> Tuple[List[str], float]:
    """Fewest contigs that together cover ``min_ref_cov`` of the reference.

    Greedy set cover over each contig's covered reference intervals
    (*eligible*): each step adds the contig that covers the most reference
    bp not yet covered (ties: circular, then shorter).  Returns ([], coverage)
    when the threshold is not reached.
    """
    chosen: List[str] = []
    covered: List[Tuple[int, int]] = []
    covered_bp = 0
    remaining = {c: merge_intervals(ivs)[0] for c, ivs in eligible.items()}
    while remaining and covered_bp < min_ref_cov * ref_len:
        gains = {c: _new_bp(ivs, covered) for c, ivs in remaining.items()}
        pick = min(gains, key=lambda c: (-gains[c], c not in circular, query_lengths.get(c, 0), c))
        if gains[pick] < _SET_MIN_GAIN * ref_len:
            break
        chosen.append(pick)
        covered = merge_intervals(covered + remaining.pop(pick))[0]
        covered_bp = _covered_bp(covered)
    frac = covered_bp / ref_len if ref_len else 0.0
    return (chosen, frac) if frac >= min_ref_cov else ([], frac)


def _new_bp(intervals: List[Tuple[int, int]], covered: List[Tuple[int, int]]) -> int:
    """Bp of sorted, merged *intervals* outside sorted, merged *covered*."""
    total, j = 0, 0
    for s, e in intervals:
        while j < len(covered) and covered[j][1] <= s:
            j += 1
        pos, k = s, j
        while k < len(covered) and covered[k][0] < e:
            cs, ce = covered[k]
            if cs > pos:
                total += cs - pos
            pos = max(pos, ce)
            k += 1
        if pos < e:
            total += e - pos
    return total


def detect_organelles(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    chrC_ref: Optional[Path],
    chrM_ref: Optional[Path],
    work_dir: Path,
    threads: int,
    thresholds: OrganelleThresholds = OrganelleThresholds(),
    circular_contigs: Optional[Set[str]] = None,
) -> OrganelleCalls:
    """BLAST contigs against chrC/chrM and call organelle contigs.

    See ``call_organelles`` for the criteria.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    if not chrC_ref and not chrM_ref:
        logger.info("No organelle references available, skipping organelle detection")
        return OrganelleCalls()

    # One BLAST subject per organelle reference, named chrC or chrC_<tag>
    # (prepare_organelle_references writes these names)
    combined_seqs: Dict[str, str] = {}
    for organelle, ref in (("chrC", chrC_ref), ("chrM", chrM_ref)):
        if ref:
            seqs = read_fasta_sequences(ref)
            for i, (name, seq) in enumerate(seqs.items()):
                if name != organelle and not name.startswith(organelle + "_"):
                    name = organelle if i == 0 else f"{organelle}_{i + 1}"
                combined_seqs[name] = seq
    ref_lengths = {name: len(seq) for name, seq in combined_seqs.items()}
    combined_ref = work_dir / "organelle_refs.fa"
    write_fasta(combined_seqs, combined_ref)

    organelle_db = work_dir / "organelle_refs"
    run_makeblastdb(combined_ref, organelle_db, err_path=work_dir / "makeblastdb.err")

    # Organelle genomes (chloroplast inverted repeats, mitochondrial repeats)
    # produce many HSPs per query-subject pair.
    blast_out = work_dir / "organelle_blast.txt"
    run_blastn_megablast(
        query_fasta=query_fasta,
        db_paths=[str(organelle_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=1000,
        err_path=work_dir / "blastn.err",
    )

    return call_organelles(
        parse_organelle_blast(blast_out),
        query_lengths,
        ref_lengths,
        thresholds=thresholds,
        circular_contigs=circular_contigs,
        genome_groups=blast_genome_grouper(
            query_fasta, query_lengths, work_dir, threads,
            min_identity=thresholds.same_genome_identity,
        ),
    )


def blast_genome_grouper(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    work_dir: Path,
    threads: int,
    min_identity: float = 0.99,
) -> GenomeGrouper:
    """Group complete candidates into genomes with megablast (see call_organelles).

    Candidates at >= ``min_identity`` over half their length are linked
    (single linkage).  Every contig assigned to the organelle counts toward the group
    of the candidate it matches best, if at that identity.
    """
    def grouper(organelle: str, candidates: List[str], assigned: List[str]):
        gdir = work_dir / f"genome_groups_{organelle}"
        gdir.mkdir(parents=True, exist_ok=True)
        cand_fa, query_fa = gdir / "candidates.fa", gdir / "assigned.fa"
        for path, names in ((cand_fa, set(candidates)), (query_fa, set(assigned) | set(candidates))):
            tmp = path.with_suffix(".fa.tmp")
            write_filtered_fasta(query_fasta, tmp, names)
            replace_if_changed(tmp, path)
        db = gdir / "candidates"
        run_makeblastdb(cand_fa, db, err_path=gdir / "makeblastdb.err")
        out = gdir / "assigned_vs_candidates.txt"
        run_blastn_megablast(
            query_fasta=query_fa, db_paths=[str(db)], output_path=out, threads=threads,
            max_hsps=100, max_target_seqs=len(candidates), err_path=gdir / "blastn.err",
        )
        return group_by_identity(parse_organelle_blast(out), candidates, assigned,
                                 query_lengths, min_identity)

    return grouper


def group_by_identity(
    pairs: Dict[Tuple[str, str], _PairEvidence],
    candidates: List[str],
    assigned: List[str],
    query_lengths: Dict[str, int],
    min_identity: float = 0.99,
) -> Tuple[Dict[str, int], Dict[int, int]]:
    """Genome groups and the bp of assigned contigs in each.

    Returns a group for every candidate and every assigned contig that
    matches one, and the summed contig length per group.
    """
    def same_genome(q: str, s: str) -> Optional[float]:
        ev = pairs.get((q, s))
        qlen = query_lengths.get(q, 0)
        if ev is None or not ev.hsps or qlen <= 0:
            return None
        ident, covered = _best_hsp_identity(ev)
        return ident if ident >= min_identity and covered / qlen >= _SAME_GENOME_MIN_COV else None

    parent = {c: c for c in candidates}

    def find(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for q in candidates:
        for s in candidates:
            if q != s and same_genome(q, s) is not None:
                parent[find(q)] = find(s)
    roots = sorted({find(c) for c in candidates})
    group_of = {c: roots.index(find(c)) for c in candidates}

    group_bp: Dict[int, int] = defaultdict(int)
    for q in sorted(set(assigned) | set(candidates)):
        if q not in group_of:
            hits = [(ident, s) for s in candidates
                    if (ident := same_genome(q, s)) is not None]
            if not hits:
                continue
            group_of[q] = group_of[max(hits)[1]]
        group_bp[group_of[q]] += query_lengths.get(q, 0)
    return group_of, dict(group_bp)


def _best_hsp_identity(ev: _PairEvidence) -> Tuple[float, int]:
    """Identity and aligned bp over non-overlapping HSPs, longest first.

    Pooling every HSP would count matches between different repeat copies
    (frequent in mitogenomes), which lowers the identity of two copies of
    one genome.
    """
    taken: List[Tuple[int, int]] = []
    matches = 0.0
    aligned = 0
    for qs, qe, ident, length, *_ in sorted(ev.hsps, key=lambda h: (-h[3], -h[2])):
        if any(qs < e and s < qe for s, e in taken):
            continue
        taken.append((qs, qe))
        matches += ident * length
        aligned += length
    return (matches / aligned if aligned else 0.0), aligned


def repick_by_read_depth(
    calls: OrganelleCalls,
    depth_median: Dict[str, float],
    query_lengths: Dict[str, int],
) -> Dict[str, Tuple[Optional[str], str]]:
    """Re-choose complete organelle contigs by read depth when candidates
    come from more than one genome.

    The genome whose contigs hold the most read bases (median depth x length,
    summed over every contig assigned to it) is the sample's own; another
    species' organelle genome in the sample is at lower depth even where it
    has as much or more assembled sequence.  Updates ``calls`` in place and
    returns organelle -> (previous complete contig, new one) for changes.
    """
    swaps: Dict[str, Tuple[Optional[str], str]] = {}
    for organelle, candidates in calls.complete_candidates.items():
        group_of = calls.genome_groups.get(organelle, {})
        read_bases: Dict[int, float] = defaultdict(float)
        for contig, group in group_of.items():
            depth = depth_median.get(contig)
            if depth:
                read_bases[group] += depth * query_lengths.get(contig, 0)
        if not read_bases:
            continue
        main = max(read_bases, key=lambda g: read_bases[g])
        logger.info(
            f"{organelle}: read bases per genome: "
            + ", ".join(
                f"{'/'.join(c for c in candidates if group_of[c] == g)[:60]} "
                f"{read_bases.get(g, 0) / 1e6:,.1f} Mb"
                for g in sorted({group_of[c] for c in candidates})
            )
        )
        new = next((c for c in candidates if group_of[c] == main), None)
        if new is None:
            continue
        calls.other_genomes[organelle] = [c for c in candidates if group_of[c] != main]
        current = getattr(calls, f"{organelle}_contig")
        if new == current:
            continue
        if current is not None:
            old_hit = calls.hits[current]
            old_hit.status = "debris"
            calls.debris.add(current)
        new_hit = calls.hits[new]
        new_hit.status = "complete"
        calls.debris.discard(new)
        setattr(calls, f"{organelle}_contig", new)
        calls.set_ref_coverage[organelle] = (
            new_hit.unit.ref_coverage if new_hit.unit else new_hit.ref_coverage)
        calls.ref_name[organelle] = new_hit.ref_name
        swaps[organelle] = (current, new)
        logger.warning(
            f"{organelle}: {new} replaces {current} as the complete contig; its genome "
            f"has the most read bases ({read_bases[main] / 1e6:,.1f} Mb)"
        )
    return swaps
