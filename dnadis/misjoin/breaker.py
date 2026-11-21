"""``dnadis misjoin break``: split contigs at the breaks of an identify TSV.

Each contig with breaks at positions p1 < p2 < ... (0-based; a break at p
cuts between bases p-1 and p) becomes pieces ``<contig>_b1``, ``_b2``, ...
in contig order; Ns at the ends of a piece (a break at a scaffold gap) are
trimmed off.  Other contigs are written unchanged, in input order.  The
AGP (v2.1) describes every output sequence as one component of an input
contig, so it maps each piece back to its original coordinates.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

from dnadis.utils.io_utils import open_maybe_gzip, open_text_writer
from dnadis.utils.logging_config import get_logger

logger = get_logger("misjoin")

DEFAULT_ACTIONS = ("break",)


@dataclass(frozen=True)
class Piece:
    name: str
    contig: str
    start: int   # 0-based, half-open
    end: int


def read_breaks(
    tsv: Path, actions: Optional[Sequence[str]] = DEFAULT_ACTIONS,
    assembly: Optional[str] = None,
) -> Dict[str, List[int]]:
    """contig -> sorted unique break positions from an identify TSV.

    Rows whose ``action`` is not in *actions* are skipped (None: use every
    row).  With *assembly*, only that assembly's rows are used; a TSV
    holding several assemblies needs it.
    """
    out: Dict[str, Set[int]] = defaultdict(set)
    seen_asm: Set[str] = set()
    with open_maybe_gzip(Path(tsv), "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asm = row.get("assembly", "")
            if assembly is not None and asm != assembly:
                continue
            seen_asm.add(asm)
            if actions is not None and row.get("action") not in actions:
                continue
            out[row["contig"]].add(int(row["break_pos"]))
    if assembly is None and len(seen_asm) > 1:
        raise ValueError(f"{tsv} holds several assemblies ({', '.join(sorted(seen_asm))}); "
                         "choose one with --assembly")
    return {c: sorted(v) for c, v in out.items()}


def plan_pieces(
    contig: str, length: int, positions: Iterable[int], seq: Optional[str] = None,
) -> List[Piece]:
    """Pieces of one contig; positions outside (0, length) are ignored.

    With *seq*, runs of N at the ends of each piece (a break at a scaffold
    gap) are trimmed off, and pieces that are all N are dropped.
    """
    cuts = sorted({p for p in positions if 0 < p < length})
    if not cuts:
        return [Piece(contig, contig, 0, length)]
    bounds = [0] + cuts + [length]
    spans = []
    for s, e in zip(bounds, bounds[1:]):
        if seq is not None:
            while s < e and seq[s] in "Nn":
                s += 1
            while e > s and seq[e - 1] in "Nn":
                e -= 1
        if e > s:
            spans.append((s, e))
    return [Piece(f"{contig}_b{i}", contig, s, e) for i, (s, e) in enumerate(spans, start=1)]


def _fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    name, buf = None, []
    with open_maybe_gzip(Path(path), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def break_fasta(
    fasta: Path,
    breaks: Dict[str, List[int]],
    out_fasta: Path,
    out_agp: Path,
    wrap: int = 80,
    threads: int = 1,
) -> List[Piece]:
    """Write the broken assembly and its AGP; returns the broken contigs' pieces."""
    names: Set[str] = set()
    broken: List[Piece] = []
    seen: Set[str] = set()
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp_agp = out_agp.with_name(out_agp.name + ".tmp")
    with open_text_writer(out_fasta, threads=threads) as fa, tmp_agp.open("w") as agp:
        agp.write("##agp-version\t2.1\n")
        agp.write(f"# pieces of {fasta.name} split at misjoin breaks\n")
        for contig, seq in _fasta_records(fasta):
            seen.add(contig)
            pieces = plan_pieces(contig, len(seq), breaks.get(contig, []), seq)
            for pc in pieces:
                if pc.name in names:
                    raise ValueError(f"duplicate output sequence name {pc.name}")
                names.add(pc.name)
                if pc.name != contig:
                    broken.append(pc)
                sub = seq[pc.start:pc.end]
                fa.write(f">{pc.name}\n")
                for i in range(0, len(sub), wrap):
                    fa.write(sub[i:i + wrap] + "\n")
                agp.write("\t".join(map(str, [
                    pc.name, 1, pc.end - pc.start, 1, "W",
                    contig, pc.start + 1, pc.end, "+",
                ])) + "\n")
    missing = sorted(set(breaks) - seen)
    if missing:
        tmp_agp.unlink()
        raise ValueError(f"contigs in the break list are not in {fasta}: {', '.join(missing[:5])}")
    tmp_agp.replace(out_agp)
    n_contigs = len({p.contig for p in broken})
    logger.done(f"Split {n_contigs} contigs into {len(broken)} pieces: {out_fasta}, {out_agp}")
    return broken
