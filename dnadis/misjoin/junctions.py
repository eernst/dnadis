"""Reference-chromosome switches along a contig, from dnadis macro blocks.

A contig's macro blocks (``<asm>.macro_blocks.tsv``) are ordered along the
contig and grouped into *runs* of consecutive blocks on one reference
chromosome.  Only blocks on nuclear chromosomes count (``kind`` =
``nuclear`` in the run's ``reference.ref_lengths.tsv``; without that
column, organelle blocks are dropped): unplaced scaffolds and organelles are
not chromosomes, so blocks on them neither split a run nor make a junction.  Runs with less than ``min_side_bp`` of aligned sequence are
dropped, smallest first, and their neighbours re-merged, so short off-target
blocks (repeats, small segments placed elsewhere by a divergent reference)
do not split a run.  Each boundary between two remaining runs on different
reference chromosomes is a *junction*: a candidate break.

The junction interval runs from the end of the last block of the left run to
the start of the first block of the right run (contig coordinates, 0-based
half-open, as in PAF).  The reference cannot place the break more precisely
than this; sister assemblies (``sisters.py``) and reads (``reads.py``)
narrow it.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Collection, Dict, Iterable, List, Optional, Tuple

from dnadis.utils.io_utils import open_maybe_gzip
from dnadis.utils.reference_utils import normalize_organelle_id


@dataclass(frozen=True)
class Block:
    """One macro block: a chain of a contig on one reference chromosome."""

    contig: str
    contig_len: int
    ref_id: str
    subgenome: str
    strand: str
    qstart: int
    qend: int
    ref_start: int
    ref_end: int
    union_bp: int

    @property
    def sort_key(self) -> Tuple[int, int, str, int]:
        return (self.qstart, self.qend, self.ref_id, self.ref_start)


@dataclass
class Run:
    """Consecutive blocks of one contig on one reference chromosome."""

    ref_id: str
    blocks: List[Block] = field(default_factory=list)

    @property
    def qstart(self) -> int:
        return min(b.qstart for b in self.blocks)

    @property
    def qend(self) -> int:
        return max(b.qend for b in self.blocks)

    @property
    def aligned_bp(self) -> int:
        return sum(b.union_bp for b in self.blocks)

    def ref_span(self) -> Tuple[int, int]:
        return min(b.ref_start for b in self.blocks), max(b.ref_end for b in self.blocks)


@dataclass
class Side:
    """One side of a junction."""

    ref_id: str
    subgenome: str
    ref_pos: int          # reference position of the block end nearest the junction
    strand: str           # strand of that block
    aligned_bp: int       # aligned bp of the whole run
    qstart: int           # run extent on the contig
    qend: int
    ref_start: int        # run extent on the reference
    ref_end: int


@dataclass
class Junction:
    """A switch between reference chromosomes along a contig."""

    contig: str
    contig_len: int
    index: int            # 1-based, in contig order
    left: Side
    right: Side
    start: int            # junction interval on the contig (0-based, half-open)
    end: int

    @property
    def junction_id(self) -> str:
        return f"{self.contig}:j{self.index}"

    @property
    def subgenomes(self) -> Tuple[str, ...]:
        return tuple(sorted({self.left.subgenome, self.right.subgenome}))


def _int(row: dict, key: str) -> int:
    return int(float(row[key]))


def read_reference_lengths(results_dir: Path) -> Tuple[Dict[str, int], Optional[set]]:
    """Reference lengths and nuclear chromosome ids of a dnadis run.

    From ``<results_dir>/reference/reference.ref_lengths.tsv``; the id set
    is None when the file or its ``kind`` column is missing.
    """
    path = Path(results_dir) / "reference" / "reference.ref_lengths.tsv"
    lengths: Dict[str, int] = {}
    nuclear: Optional[set] = None
    if not path.exists():
        return lengths, nuclear
    with open_maybe_gzip(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        has_kind = "kind" in (reader.fieldnames or [])
        if has_kind:
            nuclear = set()
        for row in reader:
            lengths[row["ref_id"]] = _int(row, "ref_len")
            if has_kind and row["kind"] == "nuclear":
                nuclear.add(row["ref_id"])
    return lengths, nuclear


def read_macro_blocks(
    path: Path, chromosomes: Optional[Collection[str]] = None,
) -> Dict[str, List[Block]]:
    """contig -> its macro blocks on nuclear chromosomes, sorted along the contig.

    With *chromosomes*, only blocks on those reference sequences are kept;
    otherwise blocks on organelle sequences are dropped (organelle-derived
    insertions in nuclear sequence are expected).
    """
    out: Dict[str, List[Block]] = defaultdict(list)
    with open_maybe_gzip(Path(path), "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if chromosomes is not None:
                if row["ref_id"] not in chromosomes:
                    continue
            elif normalize_organelle_id(row["ref_id"]):
                continue
            out[row["contig"]].append(Block(
                contig=row["contig"],
                contig_len=_int(row, "contig_len"),
                ref_id=row["ref_id"],
                subgenome=row.get("subgenome") or "",
                strand=row["strand"],
                qstart=_int(row, "qstart"),
                qend=_int(row, "qend"),
                ref_start=_int(row, "ref_start"),
                ref_end=_int(row, "ref_end"),
                union_bp=_int(row, "union_bp"),
            ))
    return {c: sorted(bs, key=lambda b: b.sort_key) for c, bs in out.items()}


def _merge_adjacent(runs: List[Run]) -> List[Run]:
    merged: List[Run] = []
    for r in runs:
        if merged and merged[-1].ref_id == r.ref_id:
            merged[-1].blocks.extend(r.blocks)
        else:
            merged.append(Run(r.ref_id, list(r.blocks)))
    return merged


def build_runs(blocks: Iterable[Block], min_block_bp: int, min_side_bp: int) -> List[Run]:
    """Group blocks into runs per reference chromosome along the contig.

    Blocks under *min_block_bp* aligned bp are ignored.  Runs under
    *min_side_bp* are removed one at a time, smallest first (ties: leftmost),
    re-merging their neighbours after each removal.
    """
    kept = sorted((b for b in blocks if b.union_bp >= min_block_bp), key=lambda b: b.sort_key)
    runs = _merge_adjacent([Run(b.ref_id, [b]) for b in kept])
    while runs:
        small = [(r.aligned_bp, r.qstart, i) for i, r in enumerate(runs) if r.aligned_bp < min_side_bp]
        if not small:
            break
        _, _, i = min(small)
        runs = _merge_adjacent(runs[:i] + runs[i + 1:])
    return runs


def _side(run: Run, nearest: Block, at_end: bool) -> Side:
    if at_end:
        pos = nearest.ref_end if nearest.strand == "+" else nearest.ref_start
    else:
        pos = nearest.ref_start if nearest.strand == "+" else nearest.ref_end
    rs, re_ = run.ref_span()
    return Side(
        ref_id=run.ref_id, subgenome=nearest.subgenome, ref_pos=pos,
        strand=nearest.strand, aligned_bp=run.aligned_bp,
        qstart=run.qstart, qend=run.qend, ref_start=rs, ref_end=re_,
    )


def find_junctions(
    blocks_by_contig: Dict[str, List[Block]],
    min_block_bp: int = 10_000,
    min_side_bp: int = 50_000,
) -> List[Junction]:
    """All reference-chromosome switches, contigs in name order."""
    out: List[Junction] = []
    for contig in sorted(blocks_by_contig):
        blocks = blocks_by_contig[contig]
        runs = build_runs(blocks, min_block_bp, min_side_bp)
        for n, (left, right) in enumerate(zip(runs, runs[1:]), start=1):
            last = max(left.blocks, key=lambda b: (b.qend, b.union_bp, b.qstart))
            first = min(right.blocks, key=lambda b: (b.qstart, -b.union_bp, b.qend))
            a, b = left.qend, right.qstart
            out.append(Junction(
                contig=contig,
                contig_len=blocks[0].contig_len,
                index=n,
                left=_side(left, last, at_end=True),
                right=_side(right, first, at_end=False),
                start=min(a, b),
                end=max(a, b),
            ))
    return out


_CHROM_RE = re.compile(r"^\D*?(\d+)")


def chromosome_number(ref_id: str) -> Optional[str]:
    """Chromosome number of a reference id (``chr05A`` -> ``5``), or None."""
    m = _CHROM_RE.match(ref_id)
    if not m:
        return None
    return m.group(1).lstrip("0") or "0"


def is_homeolog_switch(j: Junction) -> bool:
    """Both sides on the same chromosome number of different subgenomes."""
    a, b = chromosome_number(j.left.ref_id), chromosome_number(j.right.ref_id)
    return a is not None and a == b


def find_partners(junctions: List[Junction], max_ref_dist: int = 1_000_000) -> Dict[str, str]:
    """Reciprocal partner of each junction among the other contigs' junctions.

    Two junctions are partners when they join the same two reference
    chromosomes, each at a reference position within *max_ref_dist* of the
    other's (the two products of one exchange, real or a repeat swap).
    Returns {junction_id: partner junction_id}; the closest partner wins.
    """
    by_pair: Dict[Tuple[str, str], List[Junction]] = defaultdict(list)
    for j in junctions:
        by_pair[tuple(sorted((j.left.ref_id, j.right.ref_id)))].append(j)
    out: Dict[str, str] = {}
    for js in by_pair.values():
        for j in js:
            pos_j = {j.left.ref_id: j.left.ref_pos, j.right.ref_id: j.right.ref_pos}
            best: Optional[Tuple[int, str]] = None
            for k in js:
                if k.contig == j.contig:
                    continue
                pos_k = {k.left.ref_id: k.left.ref_pos, k.right.ref_id: k.right.ref_pos}
                d = max(abs(pos_j[r] - pos_k[r]) for r in pos_j)
                if d <= max_ref_dist and (best is None or (d, k.junction_id) < best):
                    best = (d, k.junction_id)
            if best:
                out[j.junction_id] = best[1]
    return out
