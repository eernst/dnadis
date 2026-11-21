"""Read evidence at a candidate break: split-read clusters and spanning reads.

Reads aligned to the assembly's input contigs (BAM/CRAM with index) are
examined around each break:

* clipped alignment ends (primary or supplementary, MAPQ >= ``min_mapq``,
  clip >= ``min_clip``) inside the break interval, grouped by position
  (consecutive ends within ``bin_bp``); the largest group is the read-derived
  break position.  It is evidence of a break only when its reads are at least
  ``min_split_share`` of the reads there (those ending in it plus those
  running through its position): in repeats, reads from other copies end
  where their copy differs while most reads run through.
* spanning reads: primary alignments (MAPQ >= ``min_mapq``) through the
  break with ``flank`` bp of sequence on each side.  For a single break the
  count is the fewest reads through any point of its interval, since a
  junction can lie anywhere in it; a repeat break needs reads that cover
  the whole shared segment.  Uniquely placed reads through a point show that
  the sequence is contiguous there in the sample.

Inside a tandem array (``tandem_run()``: k-mers recurring within a few kb,
from the contig sequence) reads from other arrays of the same repeat align
anywhere, are clipped where their copy differs, and get a MAPQ from small
copy differences.  There, a spanning read must cover the whole array plus
``flank``, and clipped ends inside the array (more than ``flank`` from its
edges) are not counted.
"""
from __future__ import annotations

from collections import Counter
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
    min_split_share: float = 0.5
    tandem_k: int = 21
    tandem_bin: int = 1_000
    tandem_window: int = 10_000
    tandem_min_frac: float = 0.5
    tandem_min_bp: int = 5_000
    tandem_max_bp: int = 500_000


@dataclass
class ReadEvidence:
    split_reads: int = 0
    split_pos: Optional[int] = None
    split_through: int = 0          # reads running through split_pos
    spanning: int = 0

    def split_supported(self, p: "ReadParams") -> bool:
        """Whether the clipped-end cluster is evidence of a break."""
        n = self.split_reads
        return n >= p.min_split and n >= p.min_split_share * (n + self.split_through)


def _clip_ends(aln) -> List[Tuple[int, int]]:
    """(position, clip length) for both ends of an alignment."""
    c = aln.cigartuples or []
    if not c:
        return []
    left = c[0][1] if c[0][0] in (4, 5) else 0
    right = c[-1][1] if c[-1][0] in (4, 5) else 0
    return [(aln.reference_start, left), (aln.reference_end, right)]


def _tandem_bin(seq: str, start: int, p: ReadParams) -> bool:
    """Whether most k-mers of the bin at *start* recur within ``tandem_window``."""
    k = p.tandem_k
    s, e = max(0, start - p.tandem_window), min(len(seq), start + p.tandem_bin + p.tandem_window)
    win = seq[s:e].upper()
    counts = Counter(win[i:i + k] for i in range(len(win) - k + 1))
    bs, be = start - s, min(start + p.tandem_bin, len(seq)) - s - k + 1
    kmers = [win[i:i + k] for i in range(bs, be)]
    kmers = [x for x in kmers if "N" not in x]
    if not kmers:
        return False
    return sum(counts[x] >= 2 for x in kmers) >= p.tandem_min_frac * len(kmers)


def tandem_run(seq: str, pos: int, p: ReadParams) -> Optional[Tuple[int, int]]:
    """The tandem array at *pos* (bins of ``tandem_bin`` bp), or None.

    A bin belongs to an array when at least ``tandem_min_frac`` of its
    k-mers occur again within ``tandem_window`` bp.  The run is grown from
    the bins at and next to *pos* by up to ``tandem_max_bp`` on each side and
    kept when at least ``tandem_min_bp`` long.
    """
    B = p.tandem_bin
    n_bins = (len(seq) + B - 1) // B
    here = min(pos // B, n_bins - 1)
    # the bin holding the break, else a neighbour (a break at the array edge)
    seed = next((b for b in (here, here - 1, here + 1)
                 if 0 <= b < n_bins and _tandem_bin(seq, b * B, p)), None)
    if seed is None:
        return None
    s = e = seed
    max_bins = p.tandem_max_bp // B
    while s > 0 and seed - s < max_bins and _tandem_bin(seq, (s - 1) * B, p):
        s -= 1
    while e + 1 < n_bins and e - seed < max_bins and _tandem_bin(seq, (e + 1) * B, p):
        e += 1
    start, end = s * B, min(len(seq), (e + 1) * B)
    return (start, end) if end - start >= p.tandem_min_bp else None


def split_clusters(
    bam, contig: str, lo: int, hi: int, p: ReadParams,
    repeat: Optional[Tuple[int, int]] = None,
) -> List[Tuple[int, List[Tuple[int, str]]]]:
    """Clusters of clipped read ends in [lo - pad, hi + pad], largest first.

    Ends inside *repeat* (a tandem array), more than ``flank`` from its
    edges, are left out.  Returns [(position, [(pos, read name), ...]), ...];
    position is the low median of the cluster.
    """
    a, b = max(0, lo - p.pad), hi + p.pad
    inner = (repeat[0] + p.flank, repeat[1] - p.flank) if repeat else None
    ends: List[Tuple[int, str]] = []
    for r in bam.fetch(contig, a, b):
        if r.is_secondary or r.is_unmapped or r.mapping_quality < p.min_mapq:
            continue
        for pos, clip in _clip_ends(r):
            if inner and inner[0] < pos < inner[1]:
                continue
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


def _through(bam, contig: str, lo: int, hi: int, p: ReadParams) -> List[Tuple[int, int]]:
    """For primary reads overlapping [lo, hi]: the positions each runs
    through with *flank* bp on both sides, as closed intervals."""
    out = []
    for r in bam.fetch(contig, max(0, lo - 1), hi + 1):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        if r.mapping_quality < p.min_mapq:
            continue
        a, b = r.reference_start + p.flank, r.reference_end - p.flank
        if a <= b:
            out.append((a, b))
    return out


def spanning_reads(bam, contig: str, lo: int, hi: int, p: ReadParams) -> int:
    """Primary reads covering [lo, hi] with *flank* bp on each side."""
    return sum(1 for a, b in _through(bam, contig, lo, hi, p) if a <= lo and b >= hi)


def tiled_reads(bam, contig: str, lo: int, hi: int, p: ReadParams) -> int:
    """The fewest reads running through any point of [lo, hi] (with *flank*
    bp on both sides of the point)."""
    events: Counter = Counter()
    for a, b in _through(bam, contig, lo, hi, p):
        a, b = max(a, lo), min(b, hi)
        if a <= b:
            events[a] += 1
            events[b + 1] -= 1
    keys = sorted(k for k in events if k <= hi)
    if not keys or keys[0] > lo:
        return 0
    depth, least = 0, None
    for i, k in enumerate(keys):
        depth += events[k]
        nxt = keys[i + 1] if i + 1 < len(keys) else hi + 1
        if nxt > k:                     # depth holds on [k, nxt - 1]
            least = depth if least is None else min(least, depth)
    return least or 0


def open_alignments(path: Path, reference: Optional[Path] = None):
    import pysam  # optional dependency, imported only when reads are given

    mode = "rc" if str(path).endswith(".cram") else "rb"
    kwargs = {"reference_filename": str(reference)} if reference and mode == "rc" else {}
    return pysam.AlignmentFile(str(path), mode, **kwargs)


def read_evidence(
    bam, contig: str, lo: int, hi: int, span_lo: int, span_hi: int, p: ReadParams,
    repeat: Optional[Tuple[int, int]] = None, whole: bool = False,
) -> Optional[ReadEvidence]:
    """Split reads in [lo, hi]; reads through [span_lo, span_hi].

    Spanning reads are the fewest through any point of [span_lo, span_hi],
    or with *whole* (a shared repeat segment) those covering all of it.
    With *repeat* (a tandem array at the break) they must also cover the
    whole array, and clipped ends inside it are not counted.  Returns None
    when *contig* is not in the alignment file.
    """
    if contig not in bam.references:
        return None
    ev = ReadEvidence()
    clusters = split_clusters(bam, contig, lo, hi, p, repeat)
    if clusters:
        pos, group = clusters[0]
        ev.split_reads = len({n for _, n in group})
        ev.split_pos = pos
        ev.split_through = spanning_reads(bam, contig, pos, pos, p)
    if whole:
        if repeat:
            span_lo, span_hi = min(span_lo, repeat[0]), max(span_hi, repeat[1])
        ev.spanning = spanning_reads(bam, contig, span_lo, span_hi, p)
    else:
        ev.spanning = tiled_reads(bam, contig, span_lo, span_hi, p)
        if repeat:
            ev.spanning = min(ev.spanning, spanning_reads(bam, contig, *repeat, p))
    return ev
