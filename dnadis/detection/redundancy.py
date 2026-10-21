#!/usr/bin/env python3
"""Own-assembly redundancy evidence for non-chromosome contigs.

Aligns every contig that is not chromosome-length (and not already called
organelle or rDNA) to the assembly's chromosome-length contigs with a
permissive preset (asm20), and records for each contig the target region
that best explains it: the largest collinear stretch of the contig matched
within one window of one target, with gap-compressed identity.

This separates contigs whose sequence is already in an assembled chromosome
(allelic copies, duplicated copies, assembly duplicates) from unique sequence
that only shares dispersed repeats with the chromosomes: repeat hits scatter
across a target and never fill one window.  Read depth of the contig and of
its matched region (``target_depth_ratio``) then tells those redundant cases
apart.
"""
from __future__ import annotations

import gzip
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from dnadis.alignment.external_tools import run_minimap2
from dnadis.utils.io_utils import (
    merge_intervals,
    open_maybe_gzip,
    output_is_current,
    replace_if_changed,
)
from dnadis.utils.logging_config import get_logger
from dnadis.utils.sequence_utils import write_filtered_fasta

logger = get_logger("redundancy")

# Target window for one collinear stretch: the contig's length times this
# factor plus a fixed allowance, so insertions/deletions between the copies
# fit while dispersed repeat hits on the same target do not add up.
_WINDOW_FACTOR = 1.5
_WINDOW_PAD = 20_000


@dataclass
class RedundancyHit:
    """Best-explaining region of one contig on one chromosome-length contig."""

    target: str
    coverage: float           # fraction of the contig aligned within the window
    identity: Optional[float]  # gap-compressed (1 - de), weighted by block length
    target_start: int
    target_end: int
    n_alignments: int


@dataclass
class _Aln:
    qs: int
    qe: int
    ts: int
    te: int
    blen: int
    de: Optional[float]


def _parse_paf(paf: Path) -> Dict[str, Tuple[int, Dict[str, List[_Aln]]]]:
    """query -> (query length, {target: primary alignments})."""
    out: Dict[str, Tuple[int, Dict[str, List[_Aln]]]] = {}
    with open_maybe_gzip(paf, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tags = {t[:2]: t[5:] for t in f[12:]}
            if tags.get("tp", "P") != "P":
                continue
            try:
                de = float(tags["de"]) if "de" in tags else None
                aln = _Aln(int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[10]), de)
                qlen = int(f[1])
            except ValueError:
                continue
            entry = out.setdefault(f[0], (qlen, defaultdict(list)))
            entry[1][f[5]].append(aln)
    return out


def _best_window(alns: List[_Aln], qlen: int) -> Tuple[float, List[_Aln]]:
    """Largest query coverage by alignments within one target window."""
    window = _WINDOW_FACTOR * qlen + _WINDOW_PAD
    alns = sorted(alns, key=lambda a: a.ts)
    best_cov, best_set = 0.0, []
    for i, first in enumerate(alns):
        members = [a for a in alns[i:] if a.ts < first.ts + window]
        _, covered = merge_intervals([(a.qs, a.qe) for a in members])
        cov = covered / qlen if qlen else 0.0
        if cov > best_cov:
            best_cov, best_set = cov, members
    return best_cov, best_set


def _gap_compressed_identity(alns: Iterable[_Aln]) -> Optional[float]:
    with_de = [a for a in alns if a.de is not None]
    total = sum(a.blen for a in with_de)
    if not total:
        return None
    return 1.0 - sum(a.de * a.blen for a in with_de) / total


def parse_redundancy_paf(paf: Path) -> Dict[str, Dict[str, RedundancyHit]]:
    """Per contig, the best window on each target contig it aligns to."""
    result: Dict[str, Dict[str, RedundancyHit]] = {}
    for query, (qlen, by_target) in _parse_paf(paf).items():
        per_target: Dict[str, RedundancyHit] = {}
        for target, alns in by_target.items():
            cov, members = _best_window(alns, qlen)
            if not members:
                continue
            per_target[target] = RedundancyHit(
                target=target,
                coverage=cov,
                identity=_gap_compressed_identity(members),
                target_start=min(a.ts for a in members),
                target_end=max(a.te for a in members),
                n_alignments=len(members),
            )
        if per_target:
            result[query] = per_target
    return result


def best_hit(per_target: Dict[str, RedundancyHit]) -> RedundancyHit:
    """The target window covering the most of the contig (identity breaks ties)."""
    return max(per_target.values(), key=lambda h: (h.coverage, h.identity or 0.0))


def detect_redundancy(
    query_fasta: Path,
    chromosome_contigs: Set[str],
    candidate_contigs: Set[str],
    work_dir: Path,
    threads: int = 1,
) -> Dict[str, Dict[str, RedundancyHit]]:
    """Align candidate contigs to chromosome-length contigs; return best windows.

    Input FASTAs are rewritten only when their content changes, and the PAF is
    reused only while it is newer than both, so a rerun with unchanged
    classes reuses the alignment.
    """
    if not chromosome_contigs or not candidate_contigs:
        return {}
    work_dir.mkdir(parents=True, exist_ok=True)
    anchors_fa = work_dir / "chromosome_contigs.fa"
    candidates_fa = work_dir / "candidates.fa"
    for names, final in ((chromosome_contigs, anchors_fa), (candidate_contigs, candidates_fa)):
        tmp = final.with_name(final.name + ".tmp")
        write_filtered_fasta(query_fasta, tmp, names)
        replace_if_changed(tmp, final)

    paf = work_dir / "candidates_vs_chrs.paf.gz"
    if paf.exists() and not output_is_current(paf, anchors_fa, candidates_fa):
        paf.unlink()
    ok = run_minimap2(
        ref=anchors_fa,
        qry=candidates_fa,
        paf_out=paf,
        threads=threads,
        preset="asm20",
        extra_args=["-c"],
        gzip_output=True,
        err_path=work_dir / "candidates_vs_chrs.err",
    )
    if not ok:
        logger.warning("Redundancy alignment failed; no redundancy evidence")
        return {}
    hits = parse_redundancy_paf(paf)
    logger.done(
        f"Redundancy evidence: {len(hits)}/{len(candidate_contigs)} contigs "
        f"align to {len(chromosome_contigs)} chromosome-length contigs"
    )
    return hits


def region_depth_median(
    regions_bed_gz: Path, regions: Dict[str, List[Tuple[int, int]]],
) -> Dict[Tuple[str, int, int], float]:
    """Median window depth over each (contig, start, end) region.

    Reads mosdepth ``*.regions.bed.gz`` once; windows are assigned to a
    region by their start position.
    """
    wanted: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for contig, ivs in regions.items():
        wanted[contig].extend(ivs)
    values: Dict[Tuple[str, int, int], List[float]] = defaultdict(list)
    with gzip.open(regions_bed_gz, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            ivs = wanted.get(f[0])
            if not ivs:
                continue
            ws, depth = int(f[1]), float(f[3])
            for s, e in ivs:
                if s <= ws < e:
                    values[(f[0], s, e)].append(depth)
    return {k: statistics.median(v) for k, v in values.items() if v}
