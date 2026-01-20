#!/usr/bin/env python3
"""
Chain parsing and evidence computation for final_finalizer.

Contains functions for parsing PAF files, building chains from alignment blocks,
and computing evidence for contig-to-reference assignments.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

from final_finalizer.models import Block, Chain, ChainEvidenceResult
from final_finalizer.utils.io_utils import merge_intervals, open_maybe_gzip
from final_finalizer.utils.reference_utils import (
    normalize_ref_id,
    paf_tag_value,
    split_chrom_subgenome,
)


def _diag_for_block(b: Block) -> int:
    return (b.rs - b.qs) if b.strand == "+" else (b.rs + b.qs)


def _rpos_for_chain_compat(b: Block) -> int:
    return b.rs if b.strand == "+" else -b.re


def _finalize_chain(chain: Chain):
    merged, qbp = merge_intervals(chain.q_intervals)
    return merged, qbp, chain.matches_sum, chain.alnlen_sum


def _chain_weight(qbp: int, matches: int, alnlen: int, mode: str) -> float:
    if alnlen <= 0:
        return 0.0
    ident = matches / alnlen
    if mode == "matches":
        return float(matches)
    if mode == "qbp_ident":
        return float(qbp) * ident
    if mode == "matches_ident":
        return float(matches) * ident
    raise ValueError(f"Unknown chain score mode: {mode}")


def _tp_keep(fields: list[str], assign_tp: str) -> bool:
    if assign_tp == "ALL":
        return True
    tp = paf_tag_value(fields, "tp:A:")
    if tp is None:
        return True
    if assign_tp == "P":
        return tp == "P"
    return tp in ("P", "I")


def parse_paf_chain_evidence_and_segments(
    paf_gz_path: Path,
    contig_lengths: dict,
    assign_minlen: int,
    assign_minmapq: int,
    assign_tp: str,
    chain_q_gap: int,
    chain_r_gap: int,
    chain_diag_slop: int,
    assign_min_ident: float,
    assign_chain_topk: int,
    assign_chain_score: str,
    assign_chain_min_bp: int,
    assign_ref_score: str,
    ref_id_map: Optional[Dict[str, str]] = None,
) -> ChainEvidenceResult:
    """
    Original nucleotide-chain mode (mm2plus PAF):
      contig=query, ref_id=reference sequence (typically chrNA/chrNP/etc.).
    """
    if assign_chain_topk < 1:
        raise ValueError("--assign-chain-topk must be >= 1")

    blocks = defaultdict(list)  # (q, ref_id, strand) -> [Block,...]
    qlens_from_paf: dict[str, int] = {}

    with open_maybe_gzip(paf_gz_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            if not _tp_keep(fields, assign_tp):
                continue

            qname = fields[0]
            try:
                qlen = int(fields[1])
                qs = int(fields[2])
                qe = int(fields[3])
                strand = fields[4]
                ref_id_raw = fields[5]
                rs = int(fields[7])
                re_ = int(fields[8])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            ref_id = ref_id_map.get(ref_id_raw, normalize_ref_id(ref_id_raw)) if ref_id_map else normalize_ref_id(ref_id_raw)

            if aln_len < assign_minlen:
                continue
            if mapq < assign_minmapq:
                continue

            if qe < qs:
                qs, qe = qe, qs
            if qe <= qs:
                continue

            if re_ < rs:
                rs, re_ = re_, rs

            qlens_from_paf.setdefault(qname, qlen)
            blocks[(qname, ref_id, strand)].append(
                Block(
                    qs=qs,
                    qe=qe,
                    rs=rs,
                    re=re_,
                    matches=matches,
                    aln_len=aln_len,
                    mapq=mapq,
                    strand=strand,
                    gene_id=None,
                )
            )

    return _chains_to_evidence_and_segments(
        blocks=blocks,
        contig_lengths=contig_lengths,
        qlens_from_paf=qlens_from_paf,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
        segments_strand_from_blocks=True,
    )


def _filter_overlapping_hits_by_identity(
    blocks: dict[tuple[str, str, str], list[Block]],
) -> tuple[dict[tuple[str, str, str], list[Block]], int, int]:
    """
    Filter overlapping hits on each contig, keeping only the highest-identity hit.

    When multiple hits overlap (1+ bp) on the same query contig, only the hit with
    the highest identity (matches/aln_len) is retained. This removes "shadow" hits
    from homeologous genes that map to the same query region with lower identity.

    Args:
        blocks: Dict mapping (contig, ref_id, strand) -> list of Block objects

    Returns:
        Tuple of (filtered_blocks, total_hits_before, total_hits_removed)
    """
    # Group all blocks by contig (across all ref_ids)
    contig_hits: dict[str, list[tuple[tuple[str, str, str], Block]]] = defaultdict(list)
    for key, blk_list in blocks.items():
        contig = key[0]
        for blk in blk_list:
            contig_hits[contig].append((key, blk))

    total_before = sum(len(v) for v in contig_hits.values())
    total_removed = 0

    filtered_blocks: dict[tuple[str, str, str], list[Block]] = defaultdict(list)

    for contig, hits in contig_hits.items():
        if not hits:
            continue

        # Calculate identity for each hit and sort by (qs, -identity)
        # so higher identity comes first when starts are equal
        hits_with_ident = []
        for key, blk in hits:
            ident = blk.matches / blk.aln_len if blk.aln_len > 0 else 0.0
            hits_with_ident.append((blk.qs, -ident, key, blk, ident))

        hits_with_ident.sort(key=lambda x: (x[0], x[1]))

        # Sweep through hits, keeping track of covered intervals
        # For each hit, check if it overlaps with any kept hit that has higher identity
        kept_intervals: list[tuple[int, int, float]] = []  # (qs, qe, identity)
        kept_hits: list[tuple[tuple[str, str, str], Block]] = []

        for qs, neg_ident, key, blk, ident in hits_with_ident:
            qe = blk.qe
            # Check if this hit overlaps with any kept interval
            dominated = False
            for kept_qs, kept_qe, kept_ident in kept_intervals:
                # Check for 1+ bp overlap
                if qs < kept_qe and qe > kept_qs:
                    # Overlap exists - keep only if this has strictly higher identity
                    if ident <= kept_ident:
                        dominated = True
                        break

            if not dominated:
                # Remove any previously kept intervals that this one dominates
                new_kept_intervals = []
                new_kept_hits = []
                for i, (kept_qs, kept_qe, kept_ident) in enumerate(kept_intervals):
                    if qs < kept_qe and qe > kept_qs and ident > kept_ident:
                        # This new hit dominates the kept one - remove it
                        total_removed += 1
                    else:
                        new_kept_intervals.append((kept_qs, kept_qe, kept_ident))
                        new_kept_hits.append(kept_hits[i])

                kept_intervals = new_kept_intervals
                kept_hits = new_kept_hits

                # Add this hit
                kept_intervals.append((qs, qe, ident))
                kept_hits.append((key, blk))
            else:
                total_removed += 1

        # Add kept hits to filtered_blocks
        for key, blk in kept_hits:
            filtered_blocks[key].append(blk)

    return filtered_blocks, total_before, total_removed


def parse_miniprot_synteny_evidence_and_segments(
    miniprot_paf_gz: Path,
    contig_lengths: dict,
    tx2loc: dict,
    tx2gene: dict,
    assign_minlen: int,
    assign_minmapq: int,
    chain_q_gap: int,
    chain_r_gap: int,
    chain_diag_slop: int,
    assign_min_ident: float,
    assign_chain_topk: int,
    assign_chain_score: str,
    assign_chain_min_bp: int,
    assign_ref_score: str,
) -> ChainEvidenceResult:
    """
    Protein-anchor mode:
      - miniprot PAF has: qname=protein/transcript, tname=query contig
      - Map qname -> ref_id(ref chrom) + (start,end,strand) from GFF3 transcript features
      - Build blocks in (query-contig coordinate) vs (reference transcript genomic coordinate) space
    """
    if assign_chain_topk < 1:
        raise ValueError("--assign-chain-topk must be >= 1")

    blocks = defaultdict(list)  # (contig, ref_id, strand) -> [Block,...]
    qlens_from_paf: dict[str, int] = {}  # contig -> tlen from PAF

    with open_maybe_gzip(miniprot_paf_gz, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            prot = fields[0]
            prot = prot.replace("mRNA:", "").replace("transcript:", "")
            tname = fields[5]  # target sequence = query contig

            try:
                tlen = int(fields[6])
                ts = int(fields[7])
                te = int(fields[8])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue

            # Map protein/transcript -> reference location
            prot_keys = [prot, prot.split()[0], prot.split("|")[0]]
            loc = None
            gene_id = None
            for pk in prot_keys:
                if pk in tx2loc:
                    loc = tx2loc.get(pk)
                    gene_id = tx2gene.get(pk) or pk
                    break
            if loc is None:
                continue

            ref_id, rs0, re0, _rstrand = loc

            if te < ts:
                ts, te = te, ts
            if te <= ts:
                continue

            # In protein mode, interpret assign_minlen as minimum TARGET SPAN on the contig.
            tspan = te - ts
            if tspan < assign_minlen:
                continue
            if mapq < assign_minmapq:
                continue

            ident = (matches / aln_len) if aln_len > 0 else 0.0
            if ident < assign_min_ident:
                continue

            # Chaining in increasing "ref coordinate" space; force '+'.
            chain_strand = "+"

            qlens_from_paf.setdefault(tname, tlen)
            blocks[(tname, ref_id, chain_strand)].append(
                Block(
                    qs=ts,
                    qe=te,
                    rs=rs0,
                    re=re0,
                    matches=matches,
                    aln_len=aln_len,
                    mapq=mapq,
                    strand=chain_strand,
                    gene_id=gene_id,
                )
            )

    # Filter overlapping hits, keeping only highest-identity hit at each position
    blocks, n_before, n_removed = _filter_overlapping_hits_by_identity(blocks)
    if n_removed > 0:
        print(
            f"[info] Filtered {n_removed}/{n_before} overlapping protein hits by identity",
            file=sys.stderr,
        )

    return _chains_to_evidence_and_segments(
        blocks=blocks,
        contig_lengths=contig_lengths,
        qlens_from_paf=qlens_from_paf,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
        segments_strand_from_blocks=False,  # forced '+'
    )


def _chains_to_evidence_and_segments(
    blocks,
    contig_lengths,
    qlens_from_paf,
    chain_q_gap,
    chain_r_gap,
    chain_diag_slop,
    assign_min_ident,
    assign_chain_topk,
    assign_chain_score,
    assign_chain_min_bp,
    assign_ref_score: str,
    segments_strand_from_blocks: bool,
) -> ChainEvidenceResult:
    # Evidence across kept chains
    qr_intervals = defaultdict(list)  # (q, ref_id) -> list[(qs,qe)]
    qr_matches = defaultdict(int)
    qr_alnlen = defaultdict(int)
    qr_gene_ids = defaultdict(set)  # (q, ref_id) -> set(gene_id)

    # Assignment weights
    qr_chain_weights = defaultdict(list)
    qr_best_chain_qbp = defaultdict(int)
    qr_best_chain_ident = defaultdict(float)
    qr_best_chain_weight = defaultdict(float)
    qr_nchains_kept = defaultdict(int)
    qr_weight_all = defaultdict(float)

    chain_segments_rows = []
    macro_block_rows = []

    for (q, ref_id, strand), blks in blocks.items():
        blks.sort(key=lambda b: (b.qs, b.qe))
        chains: list[Chain] = []

        for b in blks:
            diag = _diag_for_block(b)
            rpos = _rpos_for_chain_compat(b)

            best_i = None
            best_penalty = None

            for i, ch in enumerate(chains):
                qgap = b.qs - ch.last_qe
                if qgap > chain_q_gap:
                    continue

                rgap = rpos - ch.last_r
                if rgap < -chain_r_gap or rgap > chain_r_gap:
                    continue

                if abs(diag - ch.diag0) > chain_diag_slop:
                    continue

                penalty = abs(qgap) + abs(rgap) + (abs(diag - ch.diag0) // 10)
                if best_penalty is None or penalty < best_penalty:
                    best_penalty = penalty
                    best_i = i

            if best_i is None:
                gset: set[str] = set()
                if b.gene_id:
                    gset.add(str(b.gene_id))
                chains.append(
                    Chain(
                        last_qe=b.qe,
                        last_r=rpos,
                        diag0=diag,
                        q_intervals=[(b.qs, b.qe)],
                        matches_sum=b.matches,
                        alnlen_sum=b.aln_len,
                        gene_ids=gset,
                    )
                )
            else:
                ch = chains[best_i]
                if b.gene_id:
                    ch.gene_ids.add(str(b.gene_id))
                ch.last_qe = max(ch.last_qe, b.qe)
                ch.last_r = rpos
                ch.q_intervals.append((b.qs, b.qe))
                ch.matches_sum += b.matches
                ch.alnlen_sum += b.aln_len

        qlen = int(contig_lengths.get(q, 0) or qlens_from_paf.get(q, 0) or 0)
        chrom_id, sub = split_chrom_subgenome(ref_id)

        for chain_id, ch in enumerate(chains, start=1):
            merged, qbp, msum, alnsum = _finalize_chain(ch)

            if qbp <= 0 or alnsum <= 0 or not merged:
                continue
            if qbp < assign_chain_min_bp:
                continue

            ident = msum / alnsum
            if ident < assign_min_ident:
                continue

            w = _chain_weight(qbp=qbp, matches=msum, alnlen=alnsum, mode=assign_chain_score)
            if w <= 0:
                continue

            seg_strand = strand if segments_strand_from_blocks else "+"

            # Macro block spanning the full merged chain extent on the contig
            qstart = min(s for s, _e in merged)
            qend = max(e for _s, e in merged)
            qspan = max(0, qend - qstart)
            nseg = len(merged)
            gene_count_chain = len(ch.gene_ids) if ch.gene_ids else 0

            macro_block_rows.append(
                (
                    q,
                    qlen,
                    ref_id,
                    chrom_id,
                    sub,
                    seg_strand,
                    chain_id,
                    qstart,
                    qend,
                    qspan,
                    qbp,
                    msum,
                    alnsum,
                    f"{ident:.6f}",
                    f"{w:.3f}",
                    nseg,
                    gene_count_chain,
                )
            )

            key_qr = (q, ref_id)
            if ch.gene_ids:
                qr_gene_ids[key_qr].update(ch.gene_ids)
            qr_intervals[key_qr].extend(merged)
            qr_matches[key_qr] += msum
            qr_alnlen[key_qr] += alnsum

            qr_chain_weights[key_qr].append(w)
            qr_weight_all[key_qr] += w
            qr_nchains_kept[key_qr] += 1

            if w > qr_best_chain_weight[key_qr]:
                qr_best_chain_weight[key_qr] = w
                qr_best_chain_qbp[key_qr] = qbp
                qr_best_chain_ident[key_qr] = ident

            for (s, e) in merged:
                chain_segments_rows.append((q, qlen, chrom_id, sub, seg_strand, chain_id, s, e))

    # finalize per-(q,ref_id) union bp
    qr_union_bp = defaultdict(int)
    for key, ivs in qr_intervals.items():
        _m, total = merge_intervals(ivs)
        qr_union_bp[key] = total

    contig_total = defaultdict(int)
    contig_refs = defaultdict(set)
    for (q, ref_id), ubp in qr_union_bp.items():
        contig_total[q] += ubp
        contig_refs[q].add(ref_id)

    qr_score_topk = defaultdict(float)
    for key, ws in qr_chain_weights.items():
        ws_sorted = sorted(ws, reverse=True)
        qr_score_topk[key] = float(sum(ws_sorted[:assign_chain_topk]))

    if assign_ref_score == "all":
        qr_ref_score = qr_weight_all
    elif assign_ref_score == "topk":
        qr_ref_score = qr_score_topk
    else:
        raise ValueError(f"Unknown assign_ref_score: {assign_ref_score}")

    best_ref = defaultdict(str)
    best_score = defaultdict(float)
    best_bp = defaultdict(int)

    second_ref = defaultdict(str)
    second_score = defaultdict(float)
    second_bp = defaultdict(int)

    for (q, ref_id), score in qr_ref_score.items():
        ubp = int(qr_union_bp.get((q, ref_id), 0) or 0)
        if score > best_score[q]:
            second_score[q] = best_score[q]
            second_ref[q] = best_ref[q]
            second_bp[q] = best_bp[q]

            best_score[q] = score
            best_ref[q] = ref_id
            best_bp[q] = ubp
        elif score > second_score[q]:
            second_score[q] = score
            second_ref[q] = ref_id
            second_bp[q] = ubp

    chain_summary_rows = []
    for (q, ref_id), _ubp in qr_union_bp.items():
        ubp = int(qr_union_bp.get((q, ref_id), 0) or 0)
        msum = int(qr_matches.get((q, ref_id), 0) or 0)
        alnsum = int(qr_alnlen.get((q, ref_id), 0) or 0)
        ident = (msum / alnsum) if alnsum > 0 else 0.0
        gene_count_ref = len(qr_gene_ids.get((q, ref_id), set()))

        w_all = float(qr_weight_all.get((q, ref_id), 0.0) or 0.0)
        w_topk = float(qr_score_topk.get((q, ref_id), 0.0) or 0.0)
        nch = int(qr_nchains_kept.get((q, ref_id), 0) or 0)
        bqbp = int(qr_best_chain_qbp.get((q, ref_id), 0) or 0)
        bident = float(qr_best_chain_ident.get((q, ref_id), 0.0) or 0.0)
        bw = float(qr_best_chain_weight.get((q, ref_id), 0.0) or 0.0)

        chrom_id, sub = split_chrom_subgenome(ref_id)
        chain_summary_rows.append(
            (
                q,
                ref_id,
                chrom_id,
                sub,
                ubp,
                msum,
                alnsum,
                ident,
                gene_count_ref,
                nch,
                w_topk,
                w_all,
                bw,
                bqbp,
                bident,
            )
        )

    qr_gene_count = {k: len(v) for k, v in qr_gene_ids.items()}

    return ChainEvidenceResult(
        qlens_from_paf=dict(qlens_from_paf),
        qr_union_bp=dict(qr_union_bp),
        qr_matches=dict(qr_matches),
        qr_alnlen=dict(qr_alnlen),
        qr_gene_count=qr_gene_count,
        contig_total=dict(contig_total),
        contig_refs=dict(contig_refs),
        qr_score_topk=dict(qr_score_topk),
        qr_weight_all=dict(qr_weight_all),
        qr_nchains_kept=dict(qr_nchains_kept),
        best_ref=dict(best_ref),
        best_score=dict(best_score),
        best_bp=dict(best_bp),
        second_ref=dict(second_ref),
        second_score=dict(second_score),
        second_bp=dict(second_bp),
        chain_segments_rows=chain_segments_rows,
        macro_block_rows=macro_block_rows,
        chain_summary_rows=chain_summary_rows,
    )
