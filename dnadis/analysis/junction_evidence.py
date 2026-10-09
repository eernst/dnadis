"""Read-depth and alignment evidence at contig-internal fusion junctions.

A contig whose alignments split between two reference chromosomes is either a
real fused chromosome or an assembly misjoin.  This module locates the
junction between the two reference arms and scores three independent lines of
evidence:

- ``Aln``:   the two arms abut in the reference alignment (small unaligned gap).
- ``Depth``: long-read depth is uniform across the junction.
- ``Split``: reads span the junction rather than clipping/splitting there
  (requires the sorted BAM; opt-in).

Depth evidence needs only the persisted mosdepth ``regions.bed.gz``.
"""
from __future__ import annotations

import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from dnadis.analysis.read_depth import clipping_pileup, read_depth_windows_for_contigs
from dnadis.detection.rearrangements import unpack_macro_block_row
from dnadis.models import FusionJunction, RearrangementCall
from dnadis.utils.io_utils import merge_intervals
from dnadis.utils.logging_config import get_logger

logger = get_logger("junction_evidence")

_CONF_ORDER = ["low", "medium", "high"]

# Reads needed at the breakpoint before clipping evidence is used.
_MIN_CLIP_READS = 5
_CLIPPED_FRAC_BAD = 0.30
_CLIPPED_FRAC_OK = 0.10
_MIN_FLANK_RATIO = 0.70
_MIN_FLANK_WINDOWS = 2
_MIN_CV_WINDOWS = 3


@dataclass
class FusionThresholds:
    """Tunable thresholds for fusion-junction evidence."""
    min_fusion_frac: float = 0.30
    flank_bp: int = 5000
    pad_bp: int = 1000
    drop_ratio: float = 0.65
    spike_ratio: float = 1.5
    min_evidence_depth: float = 10.0
    max_junction_gap_bp: int = 100_000


# ---------------------------------------------------------------------------
# Junction location
# ---------------------------------------------------------------------------

def find_fusion_junctions(
    macro_block_rows: Sequence[Tuple],
    query_lengths: Dict[str, int],
    min_fusion_frac: float = 0.30,
) -> List[FusionJunction]:
    """Locate junctions between adjacent reference arms on each contig.

    A reference is an arm when its merged query span covers at least
    ``min_fusion_frac`` of the contig.  Arms are ordered by query start and
    each adjacent pair defines one junction.

    Args:
        macro_block_rows: Macro-block tuples from chain parsing.
        query_lengths: Query contig lengths (original names).
        min_fusion_frac: Minimum contig fraction for a reference to count as an arm.

    Returns:
        FusionJunction records with locations and arm spans populated.
    """
    by_contig: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: defaultdict(list))
    for row in macro_block_rows:
        b = unpack_macro_block_row(row)
        by_contig[b["contig"]][b["ref_id"]].append(b)

    junctions: List[FusionJunction] = []
    for contig, ref_blocks in by_contig.items():
        contig_len = query_lengths.get(contig) or next(
            (b["contig_len"] for bl in ref_blocks.values() for b in bl), 0)
        if contig_len <= 0:
            continue
        arms: List[Tuple[str, int, int, int]] = []  # (ref, span, q_min, q_max)
        for ref_id, blocks in ref_blocks.items():
            _, span = merge_intervals([(b["qstart"], b["qend"]) for b in blocks])
            if span / contig_len >= min_fusion_frac:
                arms.append((ref_id, span,
                             min(b["qstart"] for b in blocks),
                             max(b["qend"] for b in blocks)))
        if len(arms) < 2:
            continue
        arms.sort(key=lambda a: a[2])
        for (ref_a, span_a, _, qa_max), (ref_b, span_b, qb_min, _) in zip(arms, arms[1:]):
            if qb_min >= qa_max:
                js, je = qa_max, qb_min
                jp = (js + je) // 2
            else:  # arms overlap on the contig
                js, je = qb_min, qa_max
                jp = (js + je) // 2
            junctions.append(FusionJunction(
                contig=contig, left_ref=ref_a, right_ref=ref_b,
                junction_start=js, junction_end=je, junction_point=jp,
                contig_len=contig_len, left_span=span_a, right_span=span_b,
                aln_gap_bp=qb_min - qa_max,
            ))
    return junctions


# ---------------------------------------------------------------------------
# Depth scoring
# ---------------------------------------------------------------------------

def _cv(values: List[float]) -> Optional[float]:
    if len(values) < 2:
        return None
    mean = statistics.fmean(values)
    if mean <= 0:
        return None
    return statistics.pstdev(values) / mean


def score_depth_uniformity(
    junction: FusionJunction,
    windows_fn,
    thr: FusionThresholds,
) -> None:
    """Populate depth fields on ``junction`` from windowed depth.

    Zones are the left flank, the junction zone (the arm gap, widened to at
    least ``2 * pad_bp`` around the junction point) and the right flank.
    ``windows_fn(contig, start, end)`` returns ``(ws, we, depth)`` tuples.
    """
    jp = junction.junction_point
    j_lo = min(junction.junction_start, jp - thr.pad_bp)
    j_hi = max(junction.junction_end, jp + thr.pad_bp)
    l_lo = max(0, j_lo - thr.flank_bp)
    r_hi = min(junction.contig_len, j_hi + thr.flank_bp)

    left = [d for _, _, d in windows_fn(junction.contig, l_lo, j_lo)]
    mid = [d for _, _, d in windows_fn(junction.contig, j_lo, j_hi)]
    right = [d for _, _, d in windows_fn(junction.contig, j_hi, r_hi)]
    junction.depth_windows = len(left) + len(mid) + len(right)

    if len(left) < _MIN_FLANK_WINDOWS or len(right) < _MIN_FLANK_WINDOWS or not mid:
        junction.depth_verdict = "insufficient_depth"
        return

    l_med = statistics.median(left)
    r_med = statistics.median(right)
    m_med = statistics.median(mid)
    junction.depth_left_median = l_med
    junction.depth_right_median = r_med
    junction.depth_junction_median = m_med

    if min(l_med, r_med) < thr.min_evidence_depth:
        junction.depth_verdict = "insufficient_depth"
        return

    junction.depth_flank_ratio = min(l_med, r_med) / max(l_med, r_med)
    junction.depth_junction_ratio = m_med / statistics.fmean([l_med, r_med])
    junction.depth_flank_cv = _cv(left + right)
    junction.depth_junction_cv = _cv(mid)

    if junction.depth_junction_ratio < thr.drop_ratio:
        junction.depth_verdict = "drop"
    elif junction.depth_junction_ratio > thr.spike_ratio:
        junction.depth_verdict = "spike"
    elif (len(mid) >= _MIN_CV_WINDOWS and junction.depth_junction_cv is not None
          and junction.depth_junction_cv > 2 * (junction.depth_flank_cv or 0.0) + 0.25):
        junction.depth_verdict = "discontinuous"
    elif junction.depth_flank_ratio < _MIN_FLANK_RATIO:
        junction.depth_verdict = "asymmetric"
    else:
        junction.depth_verdict = "uniform"


# ---------------------------------------------------------------------------
# Alignment continuity (Tier A)
# ---------------------------------------------------------------------------

def assess_alignment_continuity(junction: FusionJunction, thr: FusionThresholds) -> None:
    """Classify how the two reference arms meet in the reference alignment."""
    gap = junction.aln_gap_bp
    if gap < 0:
        junction.aln_verdict = "overlapping"
    elif gap <= thr.max_junction_gap_bp:
        junction.aln_verdict = "abutting"
    else:
        junction.aln_verdict = "gapped"


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

def classify_fusion_evidence(junction: FusionJunction) -> None:
    """Set ``fusion_call``, ``fusion_confidence`` and ``evidence_codes``.

    Arms that overlap on the contig have no defined junction (the shared
    segment aligns to both references, as with repeats or organelle-derived
    insertions), so depth and read evidence there cannot discriminate fusion
    from misjoin and the call stays inconclusive.
    """
    if junction.aln_verdict == "overlapping":
        junction.fusion_call = "inconclusive"
        junction.fusion_confidence = "low"
        junction.evidence_codes = ""
        return
    depth_ok = junction.depth_verdict == "uniform"
    depth_bad = junction.depth_verdict in ("drop", "spike", "discontinuous")
    aln_ok = junction.aln_verdict == "abutting"
    aln_bad = junction.aln_verdict == "gapped"

    split_ok = split_bad = False
    if junction.clip_spanning is not None and junction.clip_clipped is not None:
        total = junction.clip_spanning + junction.clip_clipped
        if total >= _MIN_CLIP_READS:
            frac = junction.clip_clipped / total
            split_bad = frac >= _CLIPPED_FRAC_BAD
            split_ok = frac <= _CLIPPED_FRAC_OK

    asym_gap = aln_bad and junction.depth_verdict == "asymmetric"
    bad_codes = [c for c, f in (("Aln", asym_gap or (aln_bad and depth_bad)),
                                ("Depth", depth_bad), ("Split", split_bad)) if f]
    ok_codes = [c for c, f in (("Aln", aln_ok), ("Depth", depth_ok),
                               ("Split", split_ok)) if f]

    if depth_bad or split_bad or asym_gap:
        strong = sum([depth_bad, split_bad])
        conf = "high" if strong >= 2 else ("medium" if strong == 1 else "low")
        if ok_codes and conf != "low":  # conflicting evidence
            conf = _CONF_ORDER[_CONF_ORDER.index(conf) - 1]
        junction.fusion_call = "assembly_artifact"
        junction.fusion_confidence = conf
        junction.evidence_codes = ",".join(bad_codes)
    elif depth_ok or (aln_ok and split_ok):
        if depth_ok and split_ok:
            conf = "high"
        elif len(ok_codes) >= 2:
            conf = "medium"
        else:
            conf = "low"
        junction.fusion_call = "biological_fusion"
        junction.fusion_confidence = conf
        junction.evidence_codes = ",".join(ok_codes)
    else:
        junction.fusion_call = "inconclusive"
        junction.fusion_confidence = "low"
        junction.evidence_codes = ""


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def _shift_confidence(conf: str, delta: int) -> str:
    if conf not in _CONF_ORDER:
        return conf
    idx = max(0, min(len(_CONF_ORDER) - 1, _CONF_ORDER.index(conf) + delta))
    return _CONF_ORDER[idx]


def annotate_fusion_evidence(
    calls: List[RearrangementCall],
    macro_block_rows: Sequence[Tuple],
    query_lengths: Dict[str, int],
    depth_regions_path: Optional[Path],
    bam_path: Optional[Path],
    thr: FusionThresholds,
    breakpoint_tolerance: int = 50,
) -> List[FusionJunction]:
    """Score every fusion junction and write results onto the fusion calls.

    Fusion ``RearrangementCall`` objects are mutated in place: evidence fields
    are populated and ``confidence`` moves up one level for a biological
    fusion verdict and down one level for an assembly artifact.

    Args:
        calls: Rearrangement calls (non-fusion calls are ignored).
        macro_block_rows: Macro-block tuples used to locate junctions.
        query_lengths: Query contig lengths (original names).
        depth_regions_path: mosdepth regions.bed.gz, or None without ``--reads``.
        bam_path: Sorted BAM for split-read evidence, or None.
        thr: Thresholds.
        breakpoint_tolerance: Slack (bp) for the clipping pileup.

    Returns:
        Scored junctions, one per matched fusion call.
    """
    fusion_calls = [c for c in calls if c.rearrangement_type == "fusion"]
    if not fusion_calls:
        return []

    all_junctions = find_fusion_junctions(macro_block_rows, query_lengths,
                                          thr.min_fusion_frac)
    by_contig: Dict[str, List[FusionJunction]] = defaultdict(list)
    for j in all_junctions:
        by_contig[j.contig].append(j)

    # Only junctions that correspond to a reported fusion call are scored.
    matched: List[Tuple[RearrangementCall, FusionJunction]] = []
    for call in fusion_calls:
        cands = [j for j in by_contig.get(call.original_name, [])
                 if call.partner_ref_id in (j.left_ref, j.right_ref)]
        if cands:
            pick = next((j for j in cands
                         if call.assigned_ref_id in (j.left_ref, j.right_ref)),
                        cands[0])
            matched.append((call, pick))
    junctions = list({id(j): j for _, j in matched}.values())
    if not junctions:
        return []

    have_depth = depth_regions_path is not None and Path(depth_regions_path).exists()
    if have_depth:
        table = read_depth_windows_for_contigs(
            depth_regions_path, {j.contig for j in junctions})

        def windows_fn(c, s, e):
            return [w for w in table.get(c, ()) if w[1] > s and w[0] < e]
    else:
        def windows_fn(c, s, e):
            return []

    for j in junctions:
        score_depth_uniformity(j, windows_fn, thr)
        assess_alignment_continuity(j, thr)
        if bam_path is not None:
            lo, hi = ((j.junction_start, j.junction_end)
                      if j.junction_end > j.junction_start
                      else (j.junction_point, j.junction_point))
            res = clipping_pileup(bam_path, j.contig, lo, hi, breakpoint_tolerance)
            if res is not None:
                j.clip_spanning, j.clip_clipped = res
        classify_fusion_evidence(j)

    for call, pick in matched:
        call.junction_start = pick.junction_start
        call.junction_end = pick.junction_end
        call.depth_flank_ratio = pick.depth_flank_ratio
        call.depth_junction_ratio = pick.depth_junction_ratio
        call.depth_junction_cv = pick.depth_junction_cv
        call.fusion_call = pick.fusion_call
        call.fusion_confidence = pick.fusion_confidence
        call.evidence_codes = pick.evidence_codes
        if pick.fusion_call == "biological_fusion":
            call.confidence = _shift_confidence(call.confidence, +1)
        elif pick.fusion_call == "assembly_artifact":
            call.confidence = _shift_confidence(call.confidence, -1)
    return junctions
