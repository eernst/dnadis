"""Sister-assembly evidence for reference-chromosome junctions.

Each junction's neighbourhood (the junction interval plus up to ``flank`` bp
of each side's run) is cut from the assembly and aligned to every sister
assembly's chromosome contigs (``<asm>.chrs.fasta.gz``, contig-level;
minimap2 ``asm20``, CIGAR, gap-compressed identity).  A
*sister unit* is one assembly's chromosome contigs of the junction's
subgenome(s), split by ``query_subgenome`` when the assembly holds more than
one copy (e.g. two lineages of one subgenome).  The assembly the junction is
in is never its own sister.

Per unit, the sister contig with the most sequence aligned to each side's
flank is followed along the junction contig as a chain of collinear hits
(same strand, local diagonal drift within ``chain_tol``):

``carries``
    One sister contig continues from the left flank into the right flank,
    in order and orientation: the sister has the same join.
``nearby``
    Both flanks lie on one sister contig within ``near_dist`` of each other,
    but not collinear (an inversion or rearrangement near the join): the two
    chromosomes' sequences are neighbours in the sister as well.
``carries_same_genotype``
    As ``carries`` or ``nearby``, but the sister is the same genotype (gap-compressed
    identity of the chain >= ``same_genotype_identity``): an assembly built
    from the same reads or contigs inherits a misjoin, so it cannot vouch
    for the join.
``spans``
    The sister contig of one side continues collinearly across the whole run
    of the other side, which is therefore an insertion relative to the
    sister (a short run between two runs of one chromosome).
``contradicts``
    The two flanks lie on different sister contigs (or far apart on one),
    and at least one of them continues in the sister for ``min_extend`` bp
    beyond the point where it stops matching the junction contig: the
    sister holds other sequence there.
``ends``
    The sister contigs stop at the junction too (e.g. both lineages fail to
    assemble through a long repeat).
``absent``
    Neither flank aligns.

Units of any status whose chains reach ``same_genotype_identity`` are marked
``same_genotype``: independent assemblies of one inbred line (or of the same
reads) are that close to the target as well, so a same-genotype carrier may
be a structural variant of that genotype rather than an inherited misjoin.

A unit is a sister for a junction only when its chromosome contigs
(``chrom_assigned``/``chrom_fragment``) assigned to each junction
chromosome of a subgenome it holds add up to ``min_chrom_frac`` of that
chromosome's reference length: a near-empty haplotype (a few Mb of
fragments) cannot show where a chromosome continues.

The chain ends refine the break: ``x_left`` is where the left side's sister
chain stops on the junction contig, ``x_right`` where the right side's
starts.  When ``x_left <= x_right`` the break lies between them; when the
chains overlap, the overlap is sequence shared by both sides (a repeat copy
or other duplicated segment) and both of its edges are breaks.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median_low
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from dnadis.alignment.external_tools import get_minimap2_exe
from dnadis.misjoin.junctions import Junction, read_reference_lengths
from dnadis.utils.io_utils import (
    merge_intervals,
    open_maybe_gzip,
    output_is_current,
    replace_if_changed,
    run_to_gzip,
)
from dnadis.utils.logging_config import get_logger

logger = get_logger("misjoin")

CHROM_CLASSES = {"chrom_assigned", "chrom_fragment"}
_SEP = "|"


# ---------------------------------------------------------------------------
# Results-folder layout
# ---------------------------------------------------------------------------
@dataclass
class AssemblyOutputs:
    """Paths of one assembly's dnadis outputs."""

    name: str
    directory: Path

    def path(self, suffix: str) -> Path:
        return self.directory / f"{self.name}.{suffix}"

    @property
    def contig_summary(self) -> Path:
        return self.path("contig_summary.tsv")

    @property
    def macro_blocks(self) -> Path:
        return self.path("macro_blocks.tsv")

    @property
    def chrs_fasta(self) -> Optional[Path]:
        for suffix in ("chrs.fasta.gz", "chrs.fasta"):
            p = self.path(suffix)
            if p.exists():
                return p
        return None

    def class_fastas(self) -> List[Path]:
        """Contig-level class FASTAs (together they hold every input contig)."""
        out = []
        for p in sorted(self.directory.glob(f"{self.name}.*.fasta*")):
            stem = p.name[len(self.name) + 1:]
            parts = stem.split(".")
            # <class>.fasta[.gz]; skip per-subgenome splits and scaffolds
            if len(parts) in (2, 3) and parts[1] == "fasta" and parts[0] != "rdna_consensus":
                out.append(p)
        return out


def discover_assemblies(results_dir: Path) -> Dict[str, AssemblyOutputs]:
    """Assemblies under a dnadis output folder (single- or multi-assembly).

    *results_dir* may also be one assembly's own folder.
    """
    results_dir = Path(results_dir)
    found: Dict[str, AssemblyOutputs] = {}
    candidates = [results_dir] + sorted(p for p in results_dir.iterdir() if p.is_dir())
    for d in candidates:
        a = AssemblyOutputs(d.name, d)
        if a.contig_summary.exists() and a.macro_blocks.exists():
            found[d.name] = a
    return found


def read_contig_summary(path: Path) -> List[dict]:
    with open_maybe_gzip(Path(path), "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# Sister units
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SisterContig:
    name: str            # name in the sister's chrs FASTA
    subgenome: str
    group: str           # query_subgenome ('' when the assembly has one copy)
    length: int
    ref_id: str = ""     # assigned reference chromosome


@dataclass
class SisterAssembly:
    name: str
    chrs_fasta: Path
    contigs: Dict[str, SisterContig] = field(default_factory=dict)
    ref_lengths: Dict[str, int] = field(default_factory=dict)

    @property
    def subgenomes(self) -> Set[str]:
        return {c.subgenome for c in self.contigs.values()}

    def _label(self, c: SisterContig) -> str:
        return f"{self.name}_{c.group}" if c.group else self.name

    def units(self, subgenomes: Sequence[str]) -> Dict[str, Set[str]]:
        """unit label -> contig names of the given subgenomes."""
        out: Dict[str, Set[str]] = defaultdict(set)
        for c in self.contigs.values():
            if c.subgenome in subgenomes:
                out[self._label(c)].add(c.name)
        return dict(out)

    def covers(self, label: str, sides: Sequence[Tuple[str, str]], min_frac: float) -> bool:
        """Whether unit *label* holds enough of each junction chromosome.

        *sides* are (reference id, subgenome) pairs.  A side whose subgenome
        the unit does not hold, or whose reference length is unknown (an
        older run, another reference), is not tested.
        """
        mine = [c for c in self.contigs.values() if self._label(c) == label]
        held = {c.subgenome for c in mine}
        for ref, sg in sides:
            if sg not in held or ref not in self.ref_lengths:
                continue
            bp = sum(c.length for c in mine if c.ref_id == ref)
            if bp < min_frac * self.ref_lengths[ref]:
                return False
        return True


def load_sister(a: AssemblyOutputs) -> Optional[SisterAssembly]:
    fasta = a.chrs_fasta
    if fasta is None:
        return None
    s = SisterAssembly(a.name, fasta)
    for row in read_contig_summary(a.contig_summary):
        if row.get("classification") not in CHROM_CLASSES:
            continue
        s.contigs[row["contig"]] = SisterContig(
            name=row["contig"],
            subgenome=row.get("assigned_subgenome") or "",
            group=row.get("query_subgenome") or "",
            length=int(row["length"]),
            ref_id=row.get("assigned_ref_id") or "",
        )
    s.ref_lengths, _ = read_reference_lengths(a.directory.parent)
    return s if s.contigs else None


# ---------------------------------------------------------------------------
# Regions (alignment targets)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Region:
    name: str
    assembly: str
    contig: str
    start: int
    end: int
    subgenomes: Tuple[str, ...] = ()


def junction_flanks(j: Junction, flank: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Left and right flanks: up to *flank* bp of each side's run, outside the other run.

    A flank is empty (start >= end) when one run lies within the other's extent.
    """
    left = (max(j.left.qstart, j.start - flank), min(j.left.qend, j.right.qstart))
    right = (max(j.right.qstart, j.left.qend), min(j.right.qend, j.end + flank))
    return left, right


def build_regions(junctions: Dict[str, List[Junction]], flank: int) -> List[Region]:
    """Merged neighbourhoods of all junctions, per assembly and contig."""
    spans: Dict[Tuple[str, str], List[Tuple[int, int, Tuple[str, ...]]]] = defaultdict(list)
    for asm, js in junctions.items():
        for j in js:
            (ls, _), (_, re_) = junction_flanks(j, flank)
            spans[(asm, j.contig)].append(
                (max(0, ls), min(j.contig_len, max(re_, j.end)), j.subgenomes))
    regions = []
    for (asm, contig) in sorted(spans):
        merged, _ = merge_intervals([(s, e) for s, e, _ in spans[(asm, contig)]])
        for s, e in merged:
            sgs = sorted({g for a, b, gs in spans[(asm, contig)] if a < e and b > s for g in gs})
            regions.append(Region(f"{asm}{_SEP}{contig}{_SEP}{s}", asm, contig, s, e, tuple(sgs)))
    return regions


def region_for(regions: Iterable[Region], asm: str, j: Junction) -> Optional[Region]:
    for r in regions:
        if r.assembly == asm and r.contig == j.contig and r.start <= j.start and j.end <= r.end:
            return r
    return None


def write_regions_fasta(
    regions: List[Region], sequences: Dict[Tuple[str, str], str], out: Path,
) -> None:
    tmp = out.with_name(out.name + ".tmp")
    with tmp.open("w") as fh:
        for r in regions:
            seq = sequences[(r.assembly, r.contig)][r.start:r.end]
            fh.write(f">{r.name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    replace_if_changed(tmp, out)


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------
def align_sister(
    query_fasta: Path, sister_fasta: Path, paf_gz: Path, threads: int,
    copies: bool = False,
) -> Path:
    """Align junction regions (or, with *copies*, segments) to a sister.

    Query = regions cut from the junction contigs, target = the sister's
    chromosome contigs.  For copy counting, secondary alignments are kept
    (``-N 200 -p 0.1``) so every copy of a repeat inside the probe is
    reported, even where the probe's own locus scores far higher.
    Reused while newer than both inputs.
    """
    if output_is_current(paf_gz, query_fasta, sister_fasta):
        logger.debug(f"Reusing cached sister alignment: {paf_gz}")
        return paf_gz
    mapper = get_minimap2_exe()
    if not mapper:
        raise RuntimeError("Neither mm2plus nor minimap2 found in PATH")
    extra = ["-N", "200", "-p", "0.1"] if copies else ["--secondary=no"]
    cmd = [mapper, "-c", "-x", "asm20", *extra, "-t", str(threads),
           str(sister_fasta), str(query_fasta)]
    err = paf_gz.with_name(paf_gz.name.replace(".paf.gz", ".err"))
    run_to_gzip(cmd, paf_gz, err)
    return paf_gz


def align_pair(
    seq_a: str, seq_b: str, work_dir: Path, tag: str,
) -> List[Tuple[int, int, int, int, float]]:
    """Local alignments of *seq_a* (query) to *seq_b* (target), minimap2 asm5.

    Returns [(qstart, qend, tstart, tend, gap-compressed identity)].
    """
    import subprocess

    mapper = get_minimap2_exe()
    if not mapper:
        raise RuntimeError("Neither mm2plus nor minimap2 found in PATH")
    work_dir.mkdir(parents=True, exist_ok=True)
    fa, fb = work_dir / f"{tag}.a.fa", work_dir / f"{tag}.b.fa"
    fa.write_text(f">a\n{seq_a}\n")
    fb.write_text(f">b\n{seq_b}\n")
    res = subprocess.run([mapper, "-c", "-x", "asm5", "-N", "20", "-p", "0.1", str(fb), str(fa)],
                         capture_output=True, text=True, check=True)
    out = []
    for line in res.stdout.splitlines():
        t = line.split("\t")
        de = _tag(t, "de:f:")
        ident = 1.0 - float(de) if de is not None else int(t[9]) / max(1, int(t[10]))
        out.append((int(t[2]), int(t[3]), int(t[7]), int(t[8]), ident))
    return sorted(out)


@dataclass(frozen=True)
class Hit:
    """A sister alignment, in junction-contig coordinates."""

    rstart: int
    rend: int
    sister: str
    sstart: int
    send: int
    slen: int
    strand: str
    identity: float
    primary: bool

    def diag(self) -> int:
        """Sister position at contig position 0 along the alignment's line."""
        if self.strand == "+":
            return self.sstart - self.rstart
        return self.send + self.rstart

    def sister_at_right(self) -> int:
        """Sister coordinate aligned to the hit's right end on the contig."""
        return self.send if self.strand == "+" else self.sstart

    def sister_at_left(self) -> int:
        return self.sstart if self.strand == "+" else self.send


def _tag(tokens: Sequence[str], prefix: str) -> Optional[str]:
    for t in tokens[12:]:
        if t.startswith(prefix):
            return t[len(prefix):]
    return None


def parse_sister_paf(paf_gz: Path, regions: Dict[str, Region]) -> Dict[str, List[Hit]]:
    """region name -> hits (query = region, target = sister contig)."""
    out: Dict[str, List[Hit]] = defaultdict(list)
    with open_maybe_gzip(Path(paf_gz), "rt") as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            region = regions.get(t[0])
            if region is None:
                continue
            de = _tag(t, "de:f:")
            ident = 1.0 - float(de) if de is not None else int(t[9]) / max(1, int(t[10]))
            out[region.name].append(Hit(
                rstart=region.start + int(t[2]),
                rend=region.start + int(t[3]),
                sister=t[5],
                sstart=int(t[7]),
                send=int(t[8]),
                slen=int(t[6]),
                strand=t[4],
                identity=ident,
                primary=_tag(t, "tp:A:") in (None, "P"),
            ))
    for hits in out.values():
        hits.sort(key=lambda h: (h.rstart, h.rend, h.sister, h.sstart))
    return out


# ---------------------------------------------------------------------------
# Per-unit evidence
# ---------------------------------------------------------------------------
@dataclass
class SisterParams:
    flank: int = 300_000
    min_identity: float = 0.90
    min_aln: int = 5_000
    min_flank_bp: int = 10_000
    chain_tol: int = 50_000
    carry_tol: int = 200_000
    near_dist: int = 1_000_000
    span_slack: int = 5_000
    probe_flank: int = 20_000
    profile_bin: int = 100
    min_extend: int = 50_000
    same_genotype_identity: float = 0.9995
    min_chrom_frac: float = 0.5
    repeat_identity: float = 0.95
    repeat_cov: float = 0.80


@dataclass
class UnitEvidence:
    unit: str
    status: str
    left_contig: str = ""
    right_contig: str = ""
    x_left: Optional[int] = None      # left chain end on the junction contig
    x_right: Optional[int] = None     # right chain start on the junction contig
    left_continues: bool = False
    right_continues: bool = False
    identity: Optional[float] = None
    same_genotype: bool = False


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def _best_contig(hits: List[Hit], lo: int, hi: int, min_bp: int) -> Optional[str]:
    if hi <= lo:
        return None
    per: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for h in hits:
        if _overlap(h.rstart, h.rend, lo, hi):
            per[h.sister].append((max(h.rstart, lo), min(h.rend, hi)))
    best: Optional[Tuple[int, str]] = None
    for name, ivs in per.items():
        _, bp = merge_intervals(ivs)
        if bp >= min_bp and (best is None or (-bp, name) < best):
            best = (-bp, name)
    return best[1] if best else None


def _anchor(hits: List[Hit], contig: str, lo: int, hi: int) -> Hit:
    cand = [h for h in hits if h.sister == contig]
    return max(cand, key=lambda h: (_overlap(h.rstart, h.rend, lo, hi), -h.rstart, h.sstart))


def _chain(hits: List[Hit], anchor: Hit, rightward: bool, tol: int) -> List[Hit]:
    """Hits on the anchor's sister contig collinear with it, walking one way."""
    same = [h for h in hits if h.sister == anchor.sister and h.strand == anchor.strand]
    if rightward:
        seq = sorted((h for h in same if h.rstart >= anchor.rstart and h is not anchor),
                     key=lambda h: (h.rstart, h.rend, h.sstart))
    else:
        seq = sorted((h for h in same if h.rend <= anchor.rend and h is not anchor),
                     key=lambda h: (-h.rend, -h.rstart, h.sstart))
    chain = [anchor]
    for h in seq:
        if abs(h.diag() - chain[-1].diag()) <= tol:
            chain.append(h)
    return chain


def _resumes(hits: List[Hit], chain: List[Hit], beyond: int, rightward: bool, tol: int) -> bool:
    """A hit of the chain's sister contig past *beyond*, in line with the chain end."""
    if not chain:
        return False
    end = max(chain, key=lambda h: h.rend) if rightward else min(chain, key=lambda h: h.rstart)
    for h in hits:
        if h.sister != end.sister or h.strand != end.strand:
            continue
        past = h.rstart >= beyond if rightward else h.rend <= beyond
        if past and abs(h.diag() - end.diag()) <= tol:
            return True
    return False


def _chain_identity(chain: List[Hit]) -> float:
    w = sum(h.rend - h.rstart for h in chain)
    return sum(h.identity * (h.rend - h.rstart) for h in chain) / w if w else 0.0


def evaluate_unit(
    j: Junction, hits: List[Hit], unit: str, contigs: Set[str], p: SisterParams,
) -> UnitEvidence:
    """Status of one sister unit at junction *j* (see module docstring)."""
    use = [h for h in hits if h.primary and h.sister in contigs
           and h.identity >= p.min_identity and h.rend - h.rstart >= p.min_aln]
    (l0, l1), (r0, r1) = junction_flanks(j, p.flank)
    lc = _best_contig(use, l0, l1, p.min_flank_bp)
    rc = _best_contig(use, r0, r1, p.min_flank_bp)
    ev = UnitEvidence(unit, "absent", lc or "", rc or "")
    if lc is None and rc is None:
        return ev

    lchain: List[Hit] = []
    rchain: List[Hit] = []
    if lc is not None:
        lanchor = _anchor(use, lc, l0, l1)
        lchain = _chain(use, lanchor, rightward=True, tol=p.chain_tol)
        ev.x_left = max(h.rend for h in lchain)
    if rc is not None:
        ranchor = _anchor(use, rc, r0, r1)
        rchain = _chain(use, ranchor, rightward=False, tol=p.chain_tol)
        ev.x_right = min(h.rstart for h in rchain)
    ev.identity = round(_chain_identity(lchain + rchain), 5)
    ev.same_genotype = ev.identity >= p.same_genotype_identity

    if lchain and rchain and lc == rc:
        lend = max(lchain, key=lambda h: (h.rend, h.rstart))
        rstart = min(rchain, key=lambda h: (h.rstart, h.rend))
        same_geno = ev.same_genotype
        connected = ranchor in lchain or lanchor in rchain
        collinear = lanchor.strand == ranchor.strand and (
            connected or abs(lend.diag() - rstart.diag()) <= p.carry_tol)
        near = abs(lend.sister_at_right() - rstart.sister_at_left()) <= p.near_dist
        if same_geno and (collinear or near):
            ev.status = "carries_same_genotype"
            return ev
        if collinear:
            ev.status = "carries"
            return ev
        if near:
            ev.status = "nearby"
            return ev

    # one side's sister chromosome resumes beyond the other side's whole run,
    # in line with this side: that run is an insertion relative to the sister
    if _resumes(use, lchain, j.right.qend - p.span_slack, True,
                j.right.qend - j.right.qstart + p.carry_tol) or \
            _resumes(use, rchain, j.left.qstart + p.span_slack, False,
                     j.left.qend - j.left.qstart + p.carry_tol):
        ev.status = "spans"
        return ev

    if lchain:
        end = max(lchain, key=lambda h: h.rend)
        s = end.sister_at_right()
        ev.left_continues = (end.slen - s if end.strand == "+" else s) >= p.min_extend
    if rchain:
        start = min(rchain, key=lambda h: h.rstart)
        s = start.sister_at_left()
        ev.right_continues = (s if start.strand == "+" else start.slen - s) >= p.min_extend
    ev.status = "contradicts" if (ev.left_continues or ev.right_continues) else "ends"
    return ev


def copy_profile(
    hits_by_sister: Dict[str, List[Hit]], lo: int, hi: int, p: SisterParams,
) -> List[Tuple[int, int]]:
    """Median copies per sister along contig[lo:hi], in ``profile_bin`` bins.

    *hits_by_sister* are the alignments of the probe contig[lo:hi] itself,
    secondaries included (``align_sister(copies=True)``).  Each alignment at
    ``repeat_identity`` or more is one copy of the bins it covers.  Returns
    [(bin start, median copies)].
    """
    step = p.profile_bin
    n = max(0, (hi - lo + step - 1) // step)
    if n == 0 or not hits_by_sister:
        return []
    per_sister = []
    for sister in sorted(hits_by_sister):
        depth = [0] * n
        for h in hits_by_sister[sister]:
            if h.identity < p.repeat_identity:
                continue
            a = max(0, (h.rstart - lo) // step)
            b = min(n, (h.rend - lo) // step)      # bins fully covered
            for i in range(a, b):
                depth[i] += 1
        per_sister.append(depth)
    return [(lo + i * step, median_low([d[i] for d in per_sister])) for i in range(n)]


def repeat_runs(
    profile: List[Tuple[int, int]], min_copies: int, min_len: int, step: int,
) -> List[Tuple[int, int, int]]:
    """Runs of bins at *min_copies* or more, at least *min_len* long.

    Returns [(start, end, median copies over the run)].
    """
    runs: List[Tuple[int, int, int]] = []
    cur: List[Tuple[int, int]] = []
    for pos, c in profile + [(None, 0)]:
        if pos is not None and c >= min_copies and (not cur or pos == cur[-1][0] + step):
            cur.append((pos, c))
            continue
        if cur and cur[-1][0] + step - cur[0][0] >= min_len:
            runs.append((cur[0][0], cur[-1][0] + step, median_low([x[1] for x in cur])))
        cur = [(pos, c)] if pos is not None and c >= min_copies else []
    return runs


_UNIT_ORDER = {
    "carries": 0, "nearby": 1, "carries_same_genotype": 2, "spans": 3,
    "contradicts": 4, "ends": 5, "absent": 6,
}


def unit_sort_key(ev: UnitEvidence) -> Tuple[int, str]:
    return (_UNIT_ORDER.get(ev.status, 9), ev.unit)


_NAME_RE = re.compile(re.escape(_SEP))


def split_region_name(name: str) -> Tuple[str, str, int]:
    asm, contig, start = _NAME_RE.split(name, maxsplit=2)
    return asm, contig, int(start)
