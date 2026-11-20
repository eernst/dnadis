"""``dnadis misjoin identify``: candidate breaks from dnadis outputs.

Per target assembly:

1. Junctions: switches between reference chromosomes along a contig, from
   ``<asm>.macro_blocks.tsv`` (``junctions.py``).
2. Sister evidence (when the results folder, or ``--sister-results``, holds
   other assemblies): the junction neighbourhoods are aligned against every
   other assembly's chromosome contigs (``sisters.py``).
3. Read evidence (optional, ``--bam``): split-read clusters and spanning
   reads at each break (``reads.py``).
4. A call per junction and one output row per break.

Calls (``call`` column) and the default ``action``:

``misjoin`` (break)
    No sister unit carries the join and at least ``min_contradict`` units
    contradict it.
``repeat_misjoin`` (break)
    As ``misjoin``, where the sequence at the junction is shared by both
    sides (their sister chains overlap) or is, or borders, a repeat with
    ``repeat_min_copies`` or more copies per sister.  Two rows, one per edge
    of that segment: a swap through a repeat needs both edges broken.
``insertion`` (review)
    Sisters run straight across the other side's (short) run: a segment
    present here and absent from the sisters at this place.
``sister_supported`` (keep)
    At least one sister unit of another genotype carries the same join (or
    holds both sides near each other on one contig): a
    structural variant shared with that sister (or a segment placed
    differently by the reference's lineage), not a misjoin.
``unresolved`` (review)
    Sisters neither carry nor contradict the join often enough.
``reference_only`` (review)
    No sister evidence: the reference alone cannot tell a misjoin from a
    structural variant, or from a segment whose location differs between
    the reference's lineage and the assembly's.

Reads override the action when they are clear: ``min_split`` split reads at
the break with fewer than two spanning reads give ``break``; ``min_spanning``
spanning reads with fewer than ``min_split`` split reads give ``keep``.
"""
from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from dnadis.misjoin import reads as rd
from dnadis.misjoin import sisters as sis
from dnadis.misjoin.junctions import (
    Junction,
    find_junctions,
    find_partners,
    is_homeolog_switch,
    read_macro_blocks,
)
from dnadis.utils.io_utils import open_maybe_gzip
from dnadis.utils.logging_config import get_logger
from dnadis.utils.sequence_utils import reverse_complement

logger = get_logger("misjoin")

COLUMNS = [
    "assembly", "contig", "contig_len", "junction_id", "break_id",
    "break_pos", "break_lo", "break_hi", "break_kind",
    "call", "action", "confidence",
    "left_ref", "left_ref_pos", "left_strand", "left_aligned_bp",
    "left_contig_start", "left_contig_end",
    "right_ref", "right_ref_pos", "right_strand", "right_aligned_bp",
    "right_contig_start", "right_contig_end",
    "ref_gap_start", "ref_gap_end",
    "core_type", "core_start", "core_end", "core_copies",
    "n_sister_units", "n_carry", "n_nearby", "n_carry_same_genotype", "n_spans", "n_contradict",
    "n_ends", "n_absent",
    "carriers", "nearby", "same_genotype_carriers", "contradictors",
    "reads_split", "reads_split_pos", "reads_spanning",
    "partner_junction", "flags", "evidence",
]


@dataclass
class IdentifyParams:
    min_block_bp: int = 10_000
    min_side_bp: int = 50_000
    min_contradict: int = 2
    shared_min_bp: int = 1_000
    repeat_min_copies: int = 3
    partner_max_ref_dist: int = 1_000_000
    max_probe_bp: int = 150_000
    partner_identity: float = 0.99
    repeat_edge_tol: int = 2_000
    sisters: sis.SisterParams = field(default_factory=sis.SisterParams)
    reads: rd.ReadParams = field(default_factory=rd.ReadParams)


@dataclass
class Break:
    pos: int
    lo: int
    hi: int
    kind: str          # single | shared_left_edge | shared_right_edge


@dataclass
class JunctionCall:
    assembly: str
    junction: Junction
    call: str
    action: str
    confidence: str
    breaks: List[Break]
    units: List[sis.UnitEvidence] = field(default_factory=list)
    core: Optional[Tuple[str, int, int]] = None      # (type, start, end)
    core_copies: Optional[int] = None
    reads: List[Optional[rd.ReadEvidence]] = field(default_factory=list)
    partner: str = ""
    flags: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Decision
# ---------------------------------------------------------------------------
def _status(units: Sequence[sis.UnitEvidence], status: str) -> List[sis.UnitEvidence]:
    return [u for u in units if u.status == status]


def locate_breaks(
    j: Junction, units: Sequence[sis.UnitEvidence], shared_min_bp: int,
) -> Tuple[List[Break], Optional[Tuple[str, int, int]]]:
    """Break position(s) from the sister chain ends of non-carrying units.

    Each non-carrying sister shows how far the left side's chromosome runs
    (its chain end) and where the right side's begins.  A sister's chain
    stops early where it is divergent, never past the join unless the
    sequence there belongs to both sides, so the furthest left chain end and
    the earliest right chain start bound the break.  Falls back to the reference junction interval without sister chains.
    Returns (breaks, core) where core is the segment between (``gap``) or
    shared by (``shared``) the two sides' chains.
    """
    informative = [u for u in units if u.status in ("contradicts", "ends")]
    xl_all = sorted(u.x_left for u in informative if u.x_left is not None)
    xr_all = sorted(u.x_right for u in informative if u.x_right is not None)
    if not xl_all and not xr_all:
        mid = (j.start + j.end) // 2
        return [Break(mid, j.start, j.end, "single")], None
    if not xr_all:
        xl = xl_all[-1]
        return [Break(xl, min(xl, j.end), max(xl, j.end), "single")], None
    if not xl_all:
        xr = xr_all[0]
        return [Break(xr, min(j.start, xr), max(j.start, xr), "single")], None
    xl, xr = xl_all[-1], xr_all[0]
    if xr < xl and xl - xr >= shared_min_bp:
        return [
            Break(xr, xr, xr, "shared_left_edge"),
            Break(xl, xl, xl, "shared_right_edge"),
        ], ("shared", xr, xl)
    lo, hi = min(xl, xr), max(xl, xr)
    core = ("gap", lo, hi) if hi - lo >= shared_min_bp else None
    return [Break((lo + hi) // 2, lo, hi, "single")], core


def _repeat_copies(
    core: Tuple[str, int, int], repeats: Sequence[Tuple[int, int, int]], min_frac: float = 0.5,
) -> Optional[int]:
    """Copies of the repeat runs covering at least *min_frac* of *core*, or None."""
    lo, hi = core[1], core[2]
    near = [r for r in repeats if r[0] < hi and r[1] > lo]
    covered = sum(min(hi, r[1]) - max(lo, r[0]) for r in near)
    if hi > lo and covered >= min_frac * (hi - lo):
        return max(r[2] for r in near)
    return None


def _edges(lo: int, hi: int) -> List[Break]:
    return [Break(lo, lo, lo, "shared_left_edge"), Break(hi, hi, hi, "shared_right_edge")]


def place_breaks(
    j: Junction,
    units: Sequence[sis.UnitEvidence],
    repeats: Sequence[Tuple[int, int, int]],
    shared: Optional[Tuple[int, int]],
    gaps: Sequence[Tuple[int, int]],
    p: IdentifyParams,
) -> Tuple[List[Break], Optional[Tuple[str, int, int]], Optional[int]]:
    """Break positions, the segment at the junction, and its copies per sister.

    In order of precedence:

    * a scaffold gap (run of N) in or within ``repeat_edge_tol`` of the break
      interval: the join is a scaffolding join; one break at the gap;
    * *shared*: the segment this contig shares with its reciprocal partner
      at the junction (a swap through a repeat copy): breaks at both edges;
    * the sister chains (``locate_breaks``); a segment between or shared by
      them that is at least half covered by repeat runs (*repeats*) is cut
      out with a break at each edge.
    """
    breaks, core = locate_breaks(j, units, p.shared_min_bp)
    lo = min(b.lo for b in breaks)
    hi = max(b.hi for b in breaks)
    tol = p.repeat_edge_tol
    near_gaps = [g for g in gaps if g[0] < hi + tol and g[1] > lo - tol]
    if near_gaps:
        mid = (lo + hi) // 2
        g = min(near_gaps, key=lambda g: (abs((g[0] + g[1]) // 2 - mid), g[0]))
        return [Break(g[0], g[0], g[1], "scaffold_gap")], ("scaffold_gap", g[0], g[1]), None
    if shared is not None:
        core = ("partner_shared", shared[0], shared[1])
        return _edges(*shared), core, _repeat_copies(core, repeats)
    if core is None:
        return breaks, None, None
    copies = _repeat_copies(core, repeats)
    if core[0] == "gap" and copies is not None:
        core = ("repeat", core[1], core[2])
        breaks = _edges(core[1], core[2])
    return breaks, core, copies


def decide(
    asm: str,
    j: Junction,
    units: List[sis.UnitEvidence],
    p: IdentifyParams,
    repeats: Sequence[Tuple[int, int, int]] = (),
    shared: Optional[Tuple[int, int]] = None,
    gaps: Sequence[Tuple[int, int]] = (),
) -> JunctionCall:
    """Call one junction from its sister-unit evidence (see ``place_breaks``).

    *repeats* are runs of sequence around the break with a median of
    ``repeat_min_copies`` or more copies per sister: (start, end, copies).
    """
    units = sorted(units, key=sis.unit_sort_key)
    carry = _status(units, "carries") + _status(units, "nearby")
    contra = _status(units, "contradicts")
    spans = _status(units, "spans")
    present = [u for u in units if u.status != "absent"]
    breaks, core, copies = place_breaks(j, units, repeats, shared, gaps, p)

    if carry:
        call, action = "sister_supported", "keep"
        confidence = "high" if len(carry) >= 2 else "medium"
    elif not present:
        call, action, confidence = "reference_only", "review", "low"
    elif spans and len(spans) >= len(contra):
        call, action, confidence = "insertion", "review", "low"
    elif len(contra) >= p.min_contradict:
        call = "repeat_misjoin" if len(breaks) == 2 else "misjoin"
        action = "break"
        both_sides = any(u.left_continues for u in contra) and any(u.right_continues for u in contra)
        confidence = "high" if len(contra) >= 3 and both_sides and not spans else "medium"
    else:
        call, action, confidence = "unresolved", "review", "low"
    jc = JunctionCall(asm, j, call, action, confidence, breaks, units, core, copies)
    if is_homeolog_switch(j):
        jc.flags.append("homeolog_switch")
    if j.right.qend <= j.left.qend or j.left.qstart >= j.right.qstart:
        jc.flags.append("nested_runs")
    if core and core[0] == "scaffold_gap":
        jc.flags.append("scaffold_gap")
    return jc


def apply_reads(jc: JunctionCall, evs: List[Optional[rd.ReadEvidence]], p: rd.ReadParams) -> None:
    """Refine positions and override the action from read evidence."""
    jc.reads = evs
    # reads end at scaffold gaps whether or not the join is right
    seen = [e for b, e in zip(jc.breaks, evs) if e is not None and b.kind != "scaffold_gap"]
    if not seen:
        return
    split = min(e.split_reads for e in seen)
    span = max(e.spanning for e in seen)
    if split >= p.min_split and span >= 2:
        # reads both leave and cross the junction (two haplotypes, or reads
        # from another copy of a repeat): the break is not confirmed
        if jc.action == "break":
            jc.action = "review"
            jc.confidence = "low"
        jc.flags.append("reads_mixed")
        return
    for b, e in zip(jc.breaks, evs):
        if b.kind != "scaffold_gap" and e and e.split_pos is not None and e.split_reads >= p.min_split:
            b.pos = e.split_pos
            b.lo, b.hi = min(b.lo, b.pos), max(b.hi, b.pos)
    if split >= p.min_split:
        jc.action = "break"
        jc.confidence = "high"
        jc.flags.append("reads_split")
    elif span >= p.min_spanning:
        jc.action = "keep"
        jc.flags.append("reads_span")


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------
def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def _evidence_text(jc: JunctionCall) -> str:
    parts = []
    if jc.units:
        n = {s: len(_status(jc.units, s)) for s in sis._UNIT_ORDER}
        parts.append(
            f"sister units: {n['carries']} carry, {n['nearby']} nearby, "
            f"{n['carries_same_genotype']} carry (same genotype), {n['spans']} span, "
            f"{n['contradicts']} contradict, {n['ends']} end, "
            f"{n['absent']} absent"
        )
    else:
        parts.append("reference only")
    if jc.core:
        kind, lo, hi = jc.core
        what = {"shared": "shared by both sides", "gap": "between the sides",
                "repeat": "repeat between the sides",
                "partner_shared": "shared with the partner junction's contig",
                "scaffold_gap": "scaffold gap"}[kind]
        cp = f", {jc.core_copies} copies per sister" if jc.core_copies is not None else ""
        parts.append(f"{hi - lo:,} bp {what}{cp}")
    seen = [e for e in jc.reads if e is not None]
    if seen:
        parts.append("reads: " + "; ".join(
            f"{e.split_reads} split, {e.spanning} spanning" for e in seen))
    return "; ".join(parts)


def call_rows(jc: JunctionCall) -> List[Dict[str, str]]:
    j = jc.junction
    units = jc.units
    names = lambda st: ",".join(u.unit for u in _status(units, st))  # noqa: E731
    rows = []
    multi = len(jc.breaks) > 1
    for i, b in enumerate(jc.breaks):
        ev = jc.reads[i] if i < len(jc.reads) else None
        row = {
            "assembly": jc.assembly,
            "contig": j.contig,
            "contig_len": j.contig_len,
            "junction_id": j.junction_id,
            "break_id": f"{j.junction_id}{'ab'[i]}" if multi else j.junction_id,
            "break_pos": b.pos,
            "break_lo": b.lo,
            "break_hi": b.hi,
            "break_kind": b.kind,
            "call": jc.call,
            "action": jc.action,
            "confidence": jc.confidence,
            "left_ref": j.left.ref_id,
            "left_ref_pos": j.left.ref_pos,
            "left_strand": j.left.strand,
            "left_aligned_bp": j.left.aligned_bp,
            "left_contig_start": j.left.qstart,
            "left_contig_end": j.left.qend,
            "right_ref": j.right.ref_id,
            "right_ref_pos": j.right.ref_pos,
            "right_strand": j.right.strand,
            "right_aligned_bp": j.right.aligned_bp,
            "right_contig_start": j.right.qstart,
            "right_contig_end": j.right.qend,
            "ref_gap_start": j.start,
            "ref_gap_end": j.end,
            "core_type": jc.core[0] if jc.core else None,
            "core_start": jc.core[1] if jc.core else None,
            "core_end": jc.core[2] if jc.core else None,
            "core_copies": jc.core_copies,
            "n_sister_units": len(units),
            "n_carry": len(_status(units, "carries")),
            "n_nearby": len(_status(units, "nearby")),
            "n_carry_same_genotype": len(_status(units, "carries_same_genotype")),
            "n_spans": len(_status(units, "spans")),
            "n_contradict": len(_status(units, "contradicts")),
            "n_ends": len(_status(units, "ends")),
            "n_absent": len(_status(units, "absent")),
            "carriers": names("carries"),
            "nearby": names("nearby"),
            "same_genotype_carriers": names("carries_same_genotype"),
            "contradictors": names("contradicts"),
            "reads_split": ev.split_reads if ev else None,
            "reads_split_pos": ev.split_pos if ev else None,
            "reads_spanning": ev.spanning if ev else None,
            "partner_junction": jc.partner,
            "flags": ",".join(jc.flags),
            "evidence": _evidence_text(jc),
        }
        rows.append({k: _fmt(row[k]) for k in COLUMNS})
    return rows


def write_tsv(rows: List[Dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# Sequences of junction contigs
# ---------------------------------------------------------------------------
def _stream_fasta(path: Path, names) -> Dict[str, str]:
    out: Dict[str, str] = {}
    name, buf = None, []
    with open_maybe_gzip(Path(path), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in names:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name in names:
                buf.append(line.strip())
    if name in names:
        out[name] = "".join(buf)
    return out


def load_contig_sequences(
    a: sis.AssemblyOutputs, contigs: set, fasta: Optional[Path],
) -> Dict[str, str]:
    """Input-orientation sequences of *contigs* (original names).

    From the input FASTA when given; otherwise from the class FASTAs, which
    hold every contig under its output name, reverse-complemented when the
    contig summary says ``reversed``.
    """
    if fasta is not None:
        seqs = _stream_fasta(fasta, contigs)
    else:
        rows = {r["original_name"]: r for r in sis.read_contig_summary(a.contig_summary)
                if r.get("original_name") in contigs}
        by_new = {r["contig"]: orig for orig, r in rows.items()}
        seqs = {}
        for fa in a.class_fastas():
            for new, s in _stream_fasta(fa, set(by_new) - {k for k, v in by_new.items() if v in seqs}).items():
                orig = by_new[new]
                seqs[orig] = reverse_complement(s) if rows[orig].get("reversed") == "yes" else s
    missing = sorted(contigs - set(seqs))
    if missing:
        raise ValueError(f"{a.name}: sequence not found for {', '.join(missing[:5])}"
                         + (" ..." if len(missing) > 5 else ""))
    return seqs


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
@dataclass
class IdentifyInputs:
    targets: Dict[str, sis.AssemblyOutputs]
    sisters: Dict[str, sis.AssemblyOutputs]
    fastas: Dict[str, Path] = field(default_factory=dict)
    bams: Dict[str, Path] = field(default_factory=dict)
    cram_reference: Optional[Path] = None


def _run_alignments(
    plan: List[Tuple[str, Path, Path]],
    sisters: Dict[str, sis.SisterAssembly],
    threads: int,
    jobs: int,
    copies: bool,
) -> Dict[str, Path]:
    """Run (sister, query FASTA, PAF) alignments, *jobs* at a time."""
    per_job = max(1, threads // max(1, jobs))
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
        futs = {name: ex.submit(sis.align_sister, query, sisters[name].chrs_fasta, paf,
                                per_job, copies)
                for name, query, paf in plan}
        return {name: futs[name].result() for name in sorted(futs)}


def _query_files(
    sisters: Dict[str, sis.SisterAssembly],
    regions: List[sis.Region],
    seqs: Dict[Tuple[str, str], str],
    work_dir: Path,
    stem: str,
) -> List[Tuple[str, Path, Path]]:
    """Per sister: the query FASTA (regions of the other assemblies) and PAF path.

    A sister gets the regions of the other assemblies whose junctions involve
    a subgenome it has; sisters with the same set share one query file.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "paf").mkdir(exist_ok=True)
    files: Dict[Tuple[str, ...], Path] = {}
    plan = []
    for name in sorted(sisters):
        sgs = sisters[name].subgenomes
        mine = [r for r in regions if r.assembly != name and sgs.intersection(r.subgenomes)]
        if not mine:
            continue
        key = tuple(r.name for r in mine)
        if key not in files:
            path = work_dir / f"{stem}.{len(files) + 1}.fa"
            sis.write_regions_fasta(mine, seqs, path)
            files[key] = path
        plan.append((name, files[key], work_dir / "paf" / f"{stem}.{name}.paf.gz"))
    return plan


def read_gaps(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """contig -> scaffold gaps (``<asm>.gaps.bed``; empty when absent)."""
    out: Dict[str, List[Tuple[int, int]]] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3 and not line.startswith(("#", "track")):
                out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return out


def _partner_segments(
    junctions: Dict[str, List[Junction]],
    partners: Dict[str, Dict[str, str]],
    probe_of: Dict[Tuple[str, str], str],
    interval_of: Dict[Tuple[str, str], Tuple[int, int]],
    probes: Dict[str, sis.Region],
    seqs: Dict[Tuple[str, str], str],
    work_dir: Path,
    params: IdentifyParams,
) -> Dict[Tuple[str, str], Tuple[int, int]]:
    """Segment each reciprocal pair of junctions shares at the break.

    The probes (break interval +- ``probe_flank``) of the two partner
    contigs are aligned to each other; the longest alignment at
    ``partner_identity`` or more, at least ``shared_min_bp`` long and within
    ``repeat_edge_tol`` of both break intervals, is the shared segment.
    """
    out: Dict[Tuple[str, str], Tuple[int, int]] = {}
    tol = params.repeat_edge_tol
    for name in sorted(partners):
        for jid, pid in sorted(partners[name].items()):
            a, b = (name, jid), (name, pid)
            if jid > pid or a not in probe_of or b not in probe_of:
                continue
            pa, pb = probes[probe_of[a]], probes[probe_of[b]]
            hits = sis.align_pair(seqs[(name, pa.contig)][pa.start:pa.end],
                                  seqs[(name, pb.contig)][pb.start:pb.end],
                                  work_dir / "partners", f"{name}.{jid}.{pid}".replace(":", "_"))
            (alo, ahi), (blo, bhi) = interval_of[a], interval_of[b]
            best = None
            for qs, qe, ts, te, ident in hits:
                qs, qe, ts, te = qs + pa.start, qe + pa.start, ts + pb.start, te + pb.start
                if ident < params.partner_identity or qe - qs < params.shared_min_bp:
                    continue
                if not (qs < ahi + tol and qe > alo - tol and ts < bhi + tol and te > blo - tol):
                    continue
                if best is None or (qe - qs, -qs) > (best[1] - best[0], -best[0]):
                    best = (qs, qe, ts, te)
            if best:
                out[a] = (best[0], best[1])
                out[b] = (best[2], best[3])
    return out


def _sister_units(
    name: str, j: Junction, sister_asms: Dict[str, sis.SisterAssembly],
    region_hits: Dict[str, List[sis.Hit]], p: sis.SisterParams,
) -> List[sis.UnitEvidence]:
    units = []
    for sname in sorted(sister_asms):
        if sname == name:
            continue
        h = region_hits.get(sname, [])
        for label, contigs in sorted(sister_asms[sname].units(j.subgenomes).items()):
            units.append(sis.evaluate_unit(j, h, label, contigs, p))
    return units


def identify(
    inputs: IdentifyInputs,
    work_dir: Path,
    params: IdentifyParams,
    threads: int = 4,
    jobs: int = 1,
    use_sisters: bool = True,
) -> Dict[str, List[JunctionCall]]:
    """Call every junction of every target assembly."""
    junctions: Dict[str, List[Junction]] = {}
    partners: Dict[str, Dict[str, str]] = {}
    for name in sorted(inputs.targets):
        a = inputs.targets[name]
        js = find_junctions(read_macro_blocks(a.macro_blocks),
                            params.min_block_bp, params.min_side_bp)
        junctions[name] = js
        partners[name] = find_partners(js, params.partner_max_ref_dist)
        logger.info(f"{name}: {len(js)} reference-chromosome junctions in "
                    f"{len({j.contig for j in js})} contigs")

    sister_asms: Dict[str, sis.SisterAssembly] = {}
    if use_sisters:
        for name in sorted(inputs.sisters):
            s = sis.load_sister(inputs.sisters[name])
            if s is not None:
                sister_asms[name] = s
    has_sisters = any(set(sister_asms) - {t} for t in junctions)

    units: Dict[Tuple[str, str], List[sis.UnitEvidence]] = {}
    repeats: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = {}
    shared: Dict[Tuple[str, str], Tuple[int, int]] = {}
    if has_sisters and any(junctions.values()):
        sp = params.sisters
        regions = sis.build_regions(junctions, sp.flank)
        seqs: Dict[Tuple[str, str], str] = {}
        for name, js in junctions.items():
            if js:
                got = load_contig_sequences(inputs.targets[name], {j.contig for j in js},
                                            inputs.fastas.get(name))
                seqs.update({(name, c): s for c, s in got.items()})
        plan = _query_files(sister_asms, regions, seqs, work_dir, "regions")
        logger.info(f"Aligning {len(regions)} junction regions to {len(plan)} sister assemblies")
        pafs = _run_alignments(plan, sister_asms, threads, jobs, copies=False)
        by_name = {r.name: r for r in regions}
        hits = {name: sis.parse_sister_paf(paf, by_name) for name, paf in pafs.items()}

        # sister units per junction; probes around the break for repeat copies
        probes: List[sis.Region] = []
        probe_of: Dict[Tuple[str, str], str] = {}
        interval_of: Dict[Tuple[str, str], Tuple[int, int]] = {}
        for name, js in junctions.items():
            for j in js:
                region = sis.region_for(regions, name, j)
                rh = {s: hits[s].get(region.name, []) for s in hits} if region else {}
                key = (name, j.junction_id)
                units[key] = _sister_units(name, j, sister_asms, rh, sp)
                if any(u.status in ("carries", "nearby") for u in units[key]):
                    continue
                breaks, core = locate_breaks(j, units[key], params.shared_min_bp)
                lo = min([b.lo for b in breaks] + ([core[1]] if core else []))
                hi = max([b.hi for b in breaks] + ([core[2]] if core else []))
                interval_of[key] = (lo, hi)
                lo, hi = max(0, lo - sp.probe_flank), min(j.contig_len, hi + sp.probe_flank)
                if hi - lo > params.max_probe_bp:
                    continue
                r = sis.Region(f"probe{sis._SEP}{name}{sis._SEP}{j.contig}{sis._SEP}{lo}",
                               name, j.contig, lo, hi, j.subgenomes)
                if r.name not in {x.name for x in probes}:
                    probes.append(r)
                probe_of[key] = r.name
        if probes:
            plan = _query_files(sister_asms, probes, seqs, work_dir, "probes")
            logger.info(f"Counting sister copies around {len(probes)} breaks")
            pafs = _run_alignments(plan, sister_asms, threads, jobs, copies=True)
            by_name = {r.name: r for r in probes}
            probe_hits = {name: sis.parse_sister_paf(paf, by_name) for name, paf in pafs.items()}
            for key, pname in probe_of.items():
                pr = by_name[pname]
                per_sister = {s: h.get(pname, []) for s, h in probe_hits.items()
                              if s != key[0] and set(pr.subgenomes) & sister_asms[s].subgenomes}
                prof = sis.copy_profile(per_sister, pr.start, pr.end, sp)
                repeats[key] = sis.repeat_runs(prof, params.repeat_min_copies,
                                               params.shared_min_bp, sp.profile_bin)
            shared.update(_partner_segments(junctions, partners, probe_of, interval_of,
                                            by_name, seqs, work_dir, params))
        del seqs

    gaps = {name: read_gaps(inputs.targets[name].path("gaps.bed")) for name in junctions}
    out: Dict[str, List[JunctionCall]] = {}
    for name, js in junctions.items():
        bam = None
        if name in inputs.bams:
            bam = rd.open_alignments(inputs.bams[name], inputs.cram_reference)
        calls = []
        for j in js:
            key = (name, j.junction_id)
            jc = decide(name, j, units.get(key, []), params, repeats.get(key, []),
                        shared.get(key), gaps[name].get(j.contig, []))
            jc.partner = partners[name].get(j.junction_id, "")
            if bam is not None:
                evs = []
                for b in jc.breaks:
                    span = (jc.core[1], jc.core[2]) if b.kind != "single" and jc.core else (b.pos, b.pos)
                    evs.append(rd.read_evidence(bam, j.contig, b.lo, b.hi, span[0], span[1],
                                                params.reads))
                apply_reads(jc, evs, params.reads)
            calls.append(jc)
        if bam is not None:
            bam.close()
        out[name] = calls
    return out
