"""Read evidence at a candidate break: split-read clusters and spanning reads.

Reads aligned to the assembly's input contigs (BAM/CRAM with index) are
examined around each break:

* clipped alignment ends (primary or supplementary, MAPQ >= ``min_mapq``,
  clip >= ``min_clip``) inside the break interval, grouped by position
  (consecutive ends within ``bin_bp``); the largest group is the read-derived
  break position.
* spanning reads: primary alignments (MAPQ >= ``min_mapq``) that cover the
  break, or the whole shared segment of a repeat break, with ``flank`` bp of
  sequence on each side.  Uniquely placed reads through the point show that
  the sequence is contiguous in the sample.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass
class ReadParams:
    min_mapq: int = 20
    min_clip: int = 1_000
    bin_bp: int = 25
    pad: int = 5_000
    flank: int = 2_000
    min_split: int = 3
    min_spanning: int = 3


@dataclass
class ReadEvidence:
    split_reads: int = 0
    split_pos: Optional[int] = None
    spanning: int = 0


def _clip_ends(aln) -> List[Tuple[int, int]]:
    """(position, clip length) for both ends of an alignment."""
    c = aln.cigartuples or []
    if not c:
        return []
    left = c[0][1] if c[0][0] in (4, 5) else 0
    right = c[-1][1] if c[-1][0] in (4, 5) else 0
    return [(aln.reference_start, left), (aln.reference_end, right)]


def split_clusters(
    bam, contig: str, lo: int, hi: int, p: ReadParams,
) -> List[Tuple[int, List[Tuple[int, str]]]]:
    """Clusters of clipped read ends in [lo - pad, hi + pad], largest first.

    Returns [(position, [(pos, read name), ...]), ...]; position is the low
    median of the cluster.
    """
    a, b = max(0, lo - p.pad), hi + p.pad
    ends: List[Tuple[int, str]] = []
    for r in bam.fetch(contig, a, b):
        if r.is_secondary or r.is_unmapped or r.mapping_quality < p.min_mapq:
            continue
        for pos, clip in _clip_ends(r):
            if a <= pos <= b and clip >= p.min_clip:
                ends.append((pos, r.query_name))
    ends.sort()
    groups: List[List[Tuple[int, str]]] = []
    for e in ends:
        if groups and e[0] - groups[-1][-1][0] <= p.bin_bp:
            groups[-1].append(e)
        else:
            groups.append([e])
    out = []
    for g in groups:
        pos = sorted(x[0] for x in g)[(len(g) - 1) // 2]
        out.append((pos, g))
    out.sort(key=lambda x: (-len({n for _, n in x[1]}), x[0]))
    return out


def spanning_reads(bam, contig: str, lo: int, hi: int, p: ReadParams) -> int:
    """Primary reads covering [lo, hi] with *flank* bp on each side."""
    n = 0
    for r in bam.fetch(contig, max(0, lo - 1), hi + 1):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        if r.mapping_quality < p.min_mapq:
            continue
        if r.reference_start <= lo - p.flank and r.reference_end >= hi + p.flank:
            n += 1
    return n


def open_alignments(path: Path, reference: Optional[Path] = None):
    import pysam  # optional dependency, imported only when reads are given

    mode = "rc" if str(path).endswith(".cram") else "rb"
    kwargs = {"reference_filename": str(reference)} if reference and mode == "rc" else {}
    return pysam.AlignmentFile(str(path), mode, **kwargs)


def read_evidence(
    bam, contig: str, lo: int, hi: int, span_lo: int, span_hi: int, p: ReadParams,
) -> Optional[ReadEvidence]:
    """Split reads in [lo, hi]; reads spanning [span_lo, span_hi].

    Returns None when *contig* is not in the alignment file.
    """
    if contig not in bam.references:
        return None
    ev = ReadEvidence()
    clusters = split_clusters(bam, contig, lo, hi, p)
    if clusters:
        pos, group = clusters[0]
        ev.split_reads = len({n for _, n in group})
        ev.split_pos = pos
    ev.spanning = spanning_reads(bam, contig, span_lo, span_hi, p)
    return ev
