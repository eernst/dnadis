#!/usr/bin/env python3
"""
Rearrangement detection module for final_finalizer.

Classifies structural variants (translocations, inversions, fusions, fissions)
per contig from macro_block alignment evidence.  All detection is pure Python
analysis of existing macro_block data — no external tools are required.

Detection inherits the pipeline's existing macro_block minimum span
(controlled by --min-span-bp, default 50 kb).  Any macro_block that survived
chain parsing is already credible alignment evidence.  Confidence tiers
describe trust in the biological interpretation, not the alignment quality.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.models import RearrangementCall
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("rearrangements")

# Minimum off-target span (bp) to consider as translocation evidence.
_MIN_OFF_TARGET_SPAN_BP = 50_000
# Minimum fraction of contig covered by a second ref to call a fusion.
_MIN_FUSION_FRAC = 0.30


# ---------------------------------------------------------------------------
# Helper: unpack a macro_block row tuple
# ---------------------------------------------------------------------------

def _unpack_row(row: Tuple) -> dict:
    """Unpack a macro_block_rows tuple into a labelled dict.

    Supports both 17-field (legacy) and 19-field (with ref_start/ref_end)
    tuple formats.
    """
    if len(row) >= 19:
        (
            contig, contig_len, ref_id, chrom_id, subgenome, strand, chain_id,
            qstart, qend, qspan_bp, union_bp, matches, aln_len, identity,
            score, n_segments, gene_count_chain, ref_start, ref_end,
        ) = row[:19]
    else:
        (
            contig, contig_len, ref_id, chrom_id, subgenome, strand, chain_id,
            qstart, qend, qspan_bp, union_bp, matches, aln_len, identity,
            score, n_segments, gene_count_chain,
        ) = row[:17]
        ref_start = 0
        ref_end = 0
    return {
        "contig": contig,
        "contig_len": int(contig_len),
        "ref_id": ref_id,
        "strand": strand,
        "qstart": int(qstart),
        "qend": int(qend),
        "qspan_bp": int(qspan_bp),
        "union_bp": int(union_bp),
        "ref_start": int(ref_start) if ref_start != "" else 0,
        "ref_end": int(ref_end) if ref_end != "" else 0,
        "score": float(score),
        "identity": float(identity),
    }


# ---------------------------------------------------------------------------
# Confidence scoring
# ---------------------------------------------------------------------------

def _score_confidence(span_bp: int, *, reciprocal: bool = False,
                      collinear: bool = False) -> str:
    """Assign confidence level for a rearrangement call.

    Rules (from the plan):
    - high: >1 Mb AND (reciprocal or collinear)
    - medium: 50 kb–1 Mb if reciprocal/collinear, OR >1 Mb but one-sided
    - low: everything else (near threshold or ambiguous)
    """
    if span_bp > 1_000_000:
        if reciprocal or collinear:
            return "high"
        return "medium"
    if span_bp >= _MIN_OFF_TARGET_SPAN_BP:
        if reciprocal or collinear:
            return "medium"
    return "low"


# ---------------------------------------------------------------------------
# Per-contig translocation detection
# ---------------------------------------------------------------------------

def _detect_translocations(
    contig: str,
    contig_len: int,
    assigned_ref: str,
    blocks_by_ref: Dict[str, List[dict]],
    reciprocal_partners: Set[str],
) -> List[RearrangementCall]:
    """Detect translocation events for a single contig.

    For each off-target reference with total span >= threshold:
    - Check if off-target blocks are at one end (whole-arm) vs interstitial
    - Check for reciprocal pattern
    """
    calls: List[RearrangementCall] = []

    for ref_id, blocks in blocks_by_ref.items():
        if ref_id == assigned_ref:
            continue

        total_span = sum(b["qspan_bp"] for b in blocks)
        if total_span < _MIN_OFF_TARGET_SPAN_BP:
            continue

        # Aggregate query and ref coordinates
        q_min = min(b["qstart"] for b in blocks)
        q_max = max(b["qend"] for b in blocks)
        r_min = min(b["ref_start"] for b in blocks)
        r_max = max(b["ref_end"] for b in blocks)
        strand = blocks[0]["strand"]  # majority strand

        is_reciprocal = ref_id in reciprocal_partners

        # Check if it's at one end of the contig (whole-arm)
        at_start = q_min < contig_len * 0.10
        at_end = q_max > contig_len * 0.90
        is_terminal = (at_start or at_end) and not (at_start and at_end)

        if is_reciprocal:
            rtype = "reciprocal_translocation"
            evidence = (f"Off-target blocks to {ref_id} ({total_span:,} bp) "
                        f"with reciprocal pattern from partner contig")
        elif is_terminal:
            rtype = "whole_arm_translocation"
            end_label = "5'" if at_start else "3'"
            evidence = (f"Terminal off-target blocks to {ref_id} at {end_label} end "
                        f"({total_span:,} bp)")
        else:
            rtype = "translocation"
            evidence = (f"Off-target blocks to {ref_id} ({total_span:,} bp)")

        confidence = _score_confidence(total_span, reciprocal=is_reciprocal,
                                       collinear=is_terminal)

        calls.append(RearrangementCall(
            contig=contig,
            original_name=contig,  # will be updated later if renamed
            assigned_ref_id=assigned_ref,
            rearrangement_type=rtype,
            partner_ref_id=ref_id,
            query_start=q_min,
            query_end=q_max,
            ref_start=r_min,
            ref_end=r_max,
            span_bp=total_span,
            strand=strand,
            confidence=confidence,
            evidence=evidence,
        ))

    return calls


# ---------------------------------------------------------------------------
# Inversion detection
# ---------------------------------------------------------------------------

def _detect_inversions(
    contig: str,
    assigned_ref: str,
    on_target_blocks: List[dict],
) -> List[RearrangementCall]:
    """Detect inversions as strand flips within on-target blocks.

    Groups consecutive same-strand blocks and identifies strand-flip
    boundaries as inversion breakpoints.
    """
    if len(on_target_blocks) < 2:
        return []

    # Sort by query start
    sorted_blocks = sorted(on_target_blocks, key=lambda b: b["qstart"])

    # Determine majority strand (the "expected" orientation)
    plus_bp = sum(b["qspan_bp"] for b in sorted_blocks if b["strand"] == "+")
    minus_bp = sum(b["qspan_bp"] for b in sorted_blocks if b["strand"] == "-")
    majority_strand = "+" if plus_bp >= minus_bp else "-"
    inverted_strand = "-" if majority_strand == "+" else "+"

    # Find runs of inverted-strand blocks
    calls: List[RearrangementCall] = []
    inv_run: List[dict] = []

    def _flush_inv_run() -> None:
        if not inv_run:
            return
        total_span = sum(b["qspan_bp"] for b in inv_run)
        if total_span < _MIN_OFF_TARGET_SPAN_BP:
            return
        q_min = min(b["qstart"] for b in inv_run)
        q_max = max(b["qend"] for b in inv_run)
        r_min = min(b["ref_start"] for b in inv_run)
        r_max = max(b["ref_end"] for b in inv_run)
        confidence = _score_confidence(total_span, collinear=True)
        evidence = (f"Strand flip ({inverted_strand}) within {assigned_ref} "
                    f"({total_span:,} bp, {len(inv_run)} block(s))")
        calls.append(RearrangementCall(
            contig=contig,
            original_name=contig,
            assigned_ref_id=assigned_ref,
            rearrangement_type="inversion",
            partner_ref_id=None,
            query_start=q_min,
            query_end=q_max,
            ref_start=r_min,
            ref_end=r_max,
            span_bp=total_span,
            strand=inverted_strand,
            confidence=confidence,
            evidence=evidence,
        ))

    for block in sorted_blocks:
        if block["strand"] == inverted_strand:
            inv_run.append(block)
        else:
            _flush_inv_run()
            inv_run = []
    _flush_inv_run()

    return calls


# ---------------------------------------------------------------------------
# Fusion detection
# ---------------------------------------------------------------------------

def _detect_fusions(
    contig: str,
    contig_len: int,
    assigned_ref: str,
    blocks_by_ref: Dict[str, List[dict]],
    ref_lengths: Dict[str, int],
) -> List[RearrangementCall]:
    """Detect chromosome fusions: single contig covering two refs substantially."""
    calls: List[RearrangementCall] = []

    # Collect refs with significant coverage
    significant_refs: List[Tuple[str, int, int, int]] = []  # (ref_id, span, q_min, q_max)
    for ref_id, blocks in blocks_by_ref.items():
        total_span = sum(b["qspan_bp"] for b in blocks)
        frac_of_contig = total_span / contig_len if contig_len > 0 else 0
        if frac_of_contig >= _MIN_FUSION_FRAC:
            q_min = min(b["qstart"] for b in blocks)
            q_max = max(b["qend"] for b in blocks)
            significant_refs.append((ref_id, total_span, q_min, q_max))

    if len(significant_refs) < 2:
        return calls

    # Sort by query position to check for end-to-end arrangement
    significant_refs.sort(key=lambda x: x[2])

    # Check pairs of refs for fusion pattern
    for i in range(len(significant_refs)):
        for j in range(i + 1, len(significant_refs)):
            ref_a, span_a, qa_min, qa_max = significant_refs[i]
            ref_b, span_b, qb_min, qb_max = significant_refs[j]

            # The two refs should occupy complementary regions of the contig
            # (one on the "left", one on the "right")
            total_span = span_a + span_b

            # Use the non-assigned ref as partner
            if ref_a == assigned_ref:
                partner = ref_b
            elif ref_b == assigned_ref:
                partner = ref_a
            else:
                # Neither is the assigned ref — pick the one covering more
                partner = ref_b if span_a >= span_b else ref_a

            confidence = _score_confidence(total_span, collinear=True)
            evidence = (f"Contig covers {ref_a} ({span_a:,} bp) and "
                        f"{ref_b} ({span_b:,} bp) — fusion candidate")

            calls.append(RearrangementCall(
                contig=contig,
                original_name=contig,
                assigned_ref_id=assigned_ref,
                rearrangement_type="fusion",
                partner_ref_id=partner,
                query_start=min(qa_min, qb_min),
                query_end=max(qa_max, qb_max),
                ref_start=0,
                ref_end=0,
                span_bp=total_span,
                strand="+",
                confidence=confidence,
                evidence=evidence,
            ))

    return calls


# ---------------------------------------------------------------------------
# Fission detection (cross-contig, per reference)
# ---------------------------------------------------------------------------

def _detect_fissions(
    contigs_for_ref: Dict[str, List[Tuple[str, int, int, int]]],
    ref_lengths: Dict[str, int],
    best_ref: Dict[str, str],
) -> List[RearrangementCall]:
    """Detect fissions: one reference split across multiple contigs.

    Args:
        contigs_for_ref: {ref_id: [(contig, contig_len, ref_min, ref_max), ...]}
        ref_lengths: Reference chromosome lengths.
        best_ref: Contig → assigned ref mapping.
    """
    calls: List[RearrangementCall] = []

    for ref_id, contig_ranges in contigs_for_ref.items():
        if len(contig_ranges) < 2:
            continue

        ref_len = ref_lengths.get(ref_id, 0)
        if ref_len == 0:
            continue

        # Sort by ref start
        contig_ranges.sort(key=lambda x: x[2])

        # Check that contigs cover complementary, non-overlapping ref ranges
        # and together cover a substantial fraction of the reference
        total_ref_covered = 0
        for _, _, rmin, rmax in contig_ranges:
            total_ref_covered += rmax - rmin

        ref_frac = total_ref_covered / ref_len if ref_len > 0 else 0
        if ref_frac < 0.50:
            continue

        # Each individual contig must NOT cover the full ref on its own
        max_single_frac = max(
            (rmax - rmin) / ref_len for _, _, rmin, rmax in contig_ranges
        )
        if max_single_frac > 0.90:
            # One contig already covers almost all — not a fission
            continue

        span_bp = total_ref_covered
        confidence = _score_confidence(span_bp, collinear=(ref_frac > 0.80))
        contig_names = [cr[0] for cr in contig_ranges]
        evidence = (f"{len(contig_ranges)} contigs cover {ref_id} "
                    f"({ref_frac:.0%} of ref): {', '.join(contig_names)}")

        # Emit one call per contig in the fission group
        for ctg, ctg_len, rmin, rmax in contig_ranges:
            calls.append(RearrangementCall(
                contig=ctg,
                original_name=ctg,
                assigned_ref_id=ref_id,
                rearrangement_type="fission",
                partner_ref_id=None,
                query_start=0,
                query_end=ctg_len,
                ref_start=rmin,
                ref_end=rmax,
                span_bp=rmax - rmin,
                strand="+",
                confidence=confidence,
                evidence=evidence,
            ))

    return calls


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def detect_rearrangements(
    macro_block_rows: List[Tuple],
    best_ref: Dict[str, str],
    ref_lengths: Dict[str, int],
    query_lengths: Dict[str, int],
    classifications: Optional[List] = None,
) -> List[RearrangementCall]:
    """Detect structural rearrangements from macro_block evidence.

    Analyses contigs assigned to a reference chromosome (``chrom_assigned``)
    for translocation, inversion, fusion, and fission signatures.

    Args:
        macro_block_rows: Macro-block tuples from chain parsing.
        best_ref: Mapping of original contig name → assigned reference ID.
        ref_lengths: Reference chromosome lengths (normalised IDs).
        query_lengths: Query contig lengths.
        classifications: Optional list of ContigClassification objects used to
            map renamed contig names back to originals.  When provided, the
            ``contig`` field in returned calls uses the renamed name and
            ``original_name`` is set from the classification.

    Returns:
        List of RearrangementCall instances sorted by contig then query_start.
    """
    if not macro_block_rows:
        return []

    # Build name maps from classifications if available
    orig_to_new: Dict[str, str] = {}
    if classifications:
        for clf in classifications:
            orig_to_new[clf.original_name] = clf.new_name

    # ---- Group blocks by contig and then by ref ----
    # contig → ref_id → [block_dict, ...]
    contig_ref_blocks: Dict[str, Dict[str, List[dict]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in macro_block_rows:
        b = _unpack_row(row)
        contig_ref_blocks[b["contig"]][b["ref_id"]].append(b)

    # ---- Build reciprocal partner index ----
    # For each (contig assigned to refA) that has off-target blocks to refB,
    # check whether any contig assigned to refB has off-target blocks to refA.
    # off_target_refs[contig] = set of refs with significant off-target blocks
    off_target_refs: Dict[str, Set[str]] = defaultdict(set)
    for contig, ref_blocks in contig_ref_blocks.items():
        assigned = best_ref.get(contig)
        if not assigned:
            continue
        for ref_id, blocks in ref_blocks.items():
            if ref_id == assigned:
                continue
            total_span = sum(b["qspan_bp"] for b in blocks)
            if total_span >= _MIN_OFF_TARGET_SPAN_BP:
                off_target_refs[contig].add(ref_id)

    # reciprocal_partners[contig] = set of off-target refs that show reciprocity
    reciprocal_partners: Dict[str, Set[str]] = defaultdict(set)
    for contig, off_refs in off_target_refs.items():
        assigned = best_ref.get(contig)
        if not assigned:
            continue
        for off_ref in off_refs:
            # Find contigs assigned to off_ref
            for other_contig, other_assigned in best_ref.items():
                if other_assigned != off_ref:
                    continue
                if assigned in off_target_refs.get(other_contig, set()):
                    reciprocal_partners[contig].add(off_ref)
                    break

    # ---- Per-contig detection ----
    all_calls: List[RearrangementCall] = []

    # Also collect per-ref contig ranges for fission detection
    contigs_for_ref: Dict[str, List[Tuple[str, int, int, int]]] = defaultdict(list)

    for contig, ref_blocks in contig_ref_blocks.items():
        assigned = best_ref.get(contig)
        if not assigned:
            continue

        contig_len = query_lengths.get(contig, 0)
        if contig_len == 0:
            continue

        # Collect ref coordinate range for fission detection
        if assigned in ref_blocks:
            on_blocks = ref_blocks[assigned]
            r_min = min(b["ref_start"] for b in on_blocks)
            r_max = max(b["ref_end"] for b in on_blocks)
            contigs_for_ref[assigned].append((contig, contig_len, r_min, r_max))

        # Translocation detection
        all_calls.extend(_detect_translocations(
            contig, contig_len, assigned, ref_blocks,
            reciprocal_partners.get(contig, set()),
        ))

        # Inversion detection (on-target blocks only)
        on_target = ref_blocks.get(assigned, [])
        all_calls.extend(_detect_inversions(contig, assigned, on_target))

        # Fusion detection
        all_calls.extend(_detect_fusions(
            contig, contig_len, assigned, ref_blocks, ref_lengths,
        ))

    # Fission detection (cross-contig)
    all_calls.extend(_detect_fissions(contigs_for_ref, ref_lengths, best_ref))

    # ---- Map contig names (original → renamed) ----
    if orig_to_new:
        for call in all_calls:
            renamed = orig_to_new.get(call.contig)
            if renamed:
                call.original_name = call.contig
                call.contig = renamed

    # Sort by contig, then query_start
    all_calls.sort(key=lambda c: (c.contig, c.query_start))

    if all_calls:
        type_counts: Dict[str, int] = defaultdict(int)
        for c in all_calls:
            type_counts[c.rearrangement_type] += 1
        parts = [f"{k}: {v}" for k, v in sorted(type_counts.items())]
        logger.info(f"Detected {len(all_calls)} rearrangement(s): {', '.join(parts)}")
    else:
        logger.info("No rearrangements detected")

    return all_calls
