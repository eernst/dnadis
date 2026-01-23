#!/usr/bin/env python3
"""
Chain parsing and evidence computation for final_finalizer.

Contains functions for parsing PAF files, building chains from alignment blocks,
and computing evidence for contig-to-reference assignments.

Performance: Uses interval trees (intervaltree package) for O(n log n) overlap
detection during block filtering, replacing the previous O(n²) pairwise comparison.
Falls back to sweep-line algorithm if intervaltree is not available.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

try:
    from intervaltree import IntervalTree
    HAVE_INTERVALTREE = True
except ImportError:
    HAVE_INTERVALTREE = False

from final_finalizer.models import Block, Chain, ChainEvidenceResult
from final_finalizer.utils.io_utils import merge_intervals, open_maybe_gzip
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.reference_utils import (
    normalize_ref_id,
    paf_tag_value,
    split_chrom_subgenome,
)

logger = get_logger("chain_parsing")


def _diag_for_block(b: Block) -> int:
    return (b.rs - b.qs) if b.strand == "+" else (b.rs + b.qs)


def _rpos_for_chain_compat(b: Block) -> int:
    return b.rs if b.strand == "+" else -b.re


def _finalize_chain(chain: Chain):
    merged, qbp = merge_intervals(chain.q_intervals)
    return merged, qbp, chain.matches_sum, chain.alnlen_sum


def _chain_weight(qbp: int, matches: int, alnlen: int, mode: str) -> float:
    """Calculate chain weight/score based on the specified scoring mode.

    All arithmetic uses float to avoid potential integer overflow on
    very large alignments.
    """
    if alnlen <= 0:
        return 0.0
    # Ensure float arithmetic throughout
    ident = float(matches) / float(alnlen)
    if mode == "matches":
        return float(matches)
    if mode == "qbp_ident":
        return float(qbp) * ident
    if mode == "matches_ident":
        return float(matches) * ident
    raise ValueError(f"Unknown chain score mode: {mode}")


def _infer_chain_strand_from_coords(blocks: list[Block]) -> str:
    """Infer the strand orientation of a chain based on query vs reference coordinates.

    For protein-anchor chains, we determine if the query contig is forward or
    reverse-complemented relative to the reference by checking if reference
    coordinates progress in the same direction as query coordinates.

    - Forward orientation ("+"):  query and reference positions correlate positively
    - Reverse orientation ("-"): query and reference positions correlate negatively

    Args:
        blocks: List of Block objects in the chain (should be sorted by query start)

    Returns:
        "+" if forward orientation, "-" if reverse orientation
    """
    if len(blocks) < 2:
        # Can't determine strand from a single block
        # Default to "+" which means no reversal
        return "+"

    # Compute correlation between query midpoints and reference midpoints
    # using a simple sign-based voting approach
    forward_votes = 0
    reverse_votes = 0

    # Sort blocks by query start
    sorted_blocks = sorted(blocks, key=lambda b: b.qs)

    for i in range(len(sorted_blocks) - 1):
        b1 = sorted_blocks[i]
        b2 = sorted_blocks[i + 1]

        # Query direction is always positive (sorted by qs)
        q_delta = (b2.qs + b2.qe) / 2 - (b1.qs + b1.qe) / 2

        # Reference direction - compute midpoints
        r1_mid = (b1.rs + b1.re) / 2
        r2_mid = (b2.rs + b2.re) / 2
        r_delta = r2_mid - r1_mid

        # If both deltas have the same sign, it's forward orientation
        # If opposite signs, it's reverse orientation
        if (q_delta > 0 and r_delta > 0) or (q_delta < 0 and r_delta < 0):
            forward_votes += 1
        elif (q_delta > 0 and r_delta < 0) or (q_delta < 0 and r_delta > 0):
            reverse_votes += 1

    return "-" if reverse_votes > forward_votes else "+"


def _infer_contig_strand_from_all_blocks(
    all_blocks: dict[tuple[str, str, str], list[Block]],
    contig: str,
    ref_id: str,
) -> str:
    """Infer strand for a contig by looking at ALL blocks mapping to a reference.

    This is used when individual chains have too few blocks to determine strand.
    By looking at the overall correlation between query and reference positions
    across all blocks for a contig-reference pair, we can determine if the
    contig is reverse-complemented.

    Args:
        all_blocks: Dict of (contig, ref_id, strand) -> list of Block
        contig: Query contig name
        ref_id: Reference chromosome ID

    Returns:
        "+" if forward orientation, "-" if reverse orientation
    """
    # Collect all blocks for this contig-ref pair (any strand key)
    blocks = []
    for key, blks in all_blocks.items():
        if key[0] == contig and key[1] == ref_id:
            blocks.extend(blks)

    if len(blocks) < 2:
        return "+"

    # Sort by query midpoint
    sorted_blocks = sorted(blocks, key=lambda b: (b.qs + b.qe) / 2)

    # Count direction votes
    forward_votes = 0
    reverse_votes = 0

    for i in range(len(sorted_blocks) - 1):
        b1 = sorted_blocks[i]
        b2 = sorted_blocks[i + 1]

        q_mid1 = (b1.qs + b1.qe) / 2
        q_mid2 = (b2.qs + b2.qe) / 2
        r_mid1 = (b1.rs + b1.re) / 2
        r_mid2 = (b2.rs + b2.re) / 2

        q_delta = q_mid2 - q_mid1
        r_delta = r_mid2 - r_mid1

        # Skip if blocks are at same position
        if abs(q_delta) < 100:  # Within 100bp, consider same position
            continue

        if (q_delta > 0 and r_delta > 0) or (q_delta < 0 and r_delta < 0):
            forward_votes += 1
        elif (q_delta > 0 and r_delta < 0) or (q_delta < 0 and r_delta > 0):
            reverse_votes += 1

    return "-" if reverse_votes > forward_votes else "+"


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
    """Parse PAF alignment file and build synteny chains for nucleotide mode.

    Nucleotide mode uses whole-genome minimap2 alignment with permissive chaining
    to create megabase-scale synteny blocks. This is ideal for chromosome-scale
    structural composition analysis and detecting features like fusions, translocations,
    or homeologous exchanges.

    Args:
        paf_gz_path: Path to PAF file (plain or gzipped)
        contig_lengths: Query contig lengths dict
        assign_minlen: Minimum alignment length (bp) to include
        assign_minmapq: Minimum mapping quality
        assign_tp: Which alignment types to keep (P/PI/ALL)
        chain_q_gap: Max query gap between blocks in chain
        chain_r_gap: Max reference gap between blocks in chain
        chain_diag_slop: Max diagonal drift allowed in chain
        assign_min_ident: Minimum identity (matches/aln_len) threshold
        assign_chain_topk: Number of top-scoring chains to consider per contig-ref pair
        assign_chain_score: Chain scoring method (matches/qbp_ident/matches_ident)
        assign_chain_min_bp: Minimum chain union-bp to keep
        assign_ref_score: Reference scoring method (topk/all)
        ref_id_map: Optional mapping for reference ID normalization

    Returns:
        ChainEvidenceResult containing chains, segments, and scoring information

    Note:
        In nucleotide mode, contig=query assembly, ref_id=reference chromosome.
        Gate filtering requires ≥1 segment (hardcoded in cli.py) because perfect
        full-length alignments produce fewer segments than fragmented protein hits.
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

    # Warn if no alignments were found
    if not blocks:
        logger.warning(
            "No valid alignments found in PAF file after filtering. "
            "Check if reference and query are compatible, or adjust filtering thresholds "
            f"(minlen={assign_minlen}, minmapq={assign_minmapq}, min_ident={assign_min_ident})."
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
) -> tuple[dict[tuple[str, str, str], list[Block]], int, int, dict[str, int]]:
    """
    Filter overlapping hits on each contig, keeping only the highest-identity hit.

    When multiple hits overlap (1+ bp) on the same query contig, only the hit with
    the highest identity (matches/aln_len) is retained. This removes "shadow" hits
    from homeologous genes that map to the same query region with lower identity.

    This is biologically important for polyploid genomes where homeologous gene
    copies may map to the same query region with different identities.

    Uses an interval tree for O(n log n) performance instead of O(n²) pairwise comparison.

    Args:
        blocks: Dict mapping (contig, ref_id, strand) -> list of Block objects

    Returns:
        Tuple of (filtered_blocks, total_hits_before, total_hits_removed, removed_per_contig)
    """
    # Group all blocks by contig (across all ref_ids)
    contig_hits: dict[str, list[tuple[tuple[str, str, str], Block]]] = defaultdict(list)
    for key, blk_list in blocks.items():
        contig = key[0]
        for blk in blk_list:
            contig_hits[contig].append((key, blk))

    total_before = sum(len(v) for v in contig_hits.values())
    total_removed = 0
    removed_per_contig: dict[str, int] = defaultdict(int)

    filtered_blocks: dict[tuple[str, str, str], list[Block]] = defaultdict(list)

    for contig, hits in contig_hits.items():
        if not hits:
            continue

        # Calculate identity for each hit and sort by identity descending (highest first)
        hits_with_ident = []
        for key, blk in hits:
            ident = blk.matches / blk.aln_len if blk.aln_len > 0 else 0.0
            hits_with_ident.append((blk.qs, blk.qe, ident, key, blk))

        # Sort by identity descending - process highest identity hits first
        hits_with_ident.sort(key=lambda x: -x[2])

        if HAVE_INTERVALTREE:
            # Use interval tree for O(n log n) performance
            tree: IntervalTree = IntervalTree()
            kept_hits: list[tuple[tuple[str, str, str], Block]] = []

            for qs, qe, ident, key, blk in hits_with_ident:
                # Ensure valid interval (qe > qs)
                if qe <= qs:
                    continue

                # Check if this hit overlaps with any existing higher-identity hit
                # Since we process in identity order, any existing interval has higher/equal identity
                overlapping = tree[qs:qe]

                if not overlapping:
                    # No overlapping higher-identity hits - keep this one
                    tree[qs:qe] = (ident, key, blk)
                    kept_hits.append((key, blk))
                else:
                    # Overlaps with higher-identity hit - filter this one out
                    total_removed += 1
                    removed_per_contig[contig] += 1
        else:
            # Fallback: O(n²) sweep line algorithm if intervaltree not available
            kept_intervals: list[tuple[int, int, float]] = []
            kept_hits = []

            for qs, qe, ident, key, blk in hits_with_ident:
                if qe <= qs:
                    continue

                # Check if this hit overlaps with any kept interval
                dominated = False
                for kept_qs, kept_qe, kept_ident in kept_intervals:
                    # Check for 1+ bp overlap
                    if qs < kept_qe and qe > kept_qs:
                        # Overlap exists - since we process in identity order,
                        # the kept interval has higher or equal identity
                        dominated = True
                        break

                if not dominated:
                    kept_intervals.append((qs, qe, ident))
                    kept_hits.append((key, blk))
                else:
                    total_removed += 1
                    removed_per_contig[contig] += 1

        # Add kept hits to filtered_blocks
        for key, blk in kept_hits:
            filtered_blocks[key].append(blk)

    return filtered_blocks, total_before, total_removed, dict(removed_per_contig)


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
    """Parse miniprot PAF file and build synteny chains for protein mode.

    Protein mode uses protein-anchored synteny based on miniprot alignments.
    Proteins are more conserved than nucleotide sequences, making this mode
    ideal for gene-level classification across distantly related species.

    Args:
        miniprot_paf_gz: Path to miniprot PAF file (plain or gzipped)
        contig_lengths: Query contig lengths dict
        tx2loc: Mapping of transcript ID -> (ref_id, start, end, strand) from GFF3
        tx2gene: Mapping of transcript ID -> gene ID from GFF3
        assign_minlen: Minimum target span (bp) on query contig
        assign_minmapq: Minimum mapping quality
        chain_q_gap: Max query gap between blocks in chain
        chain_r_gap: Max reference gap between blocks in chain
        chain_diag_slop: Max diagonal drift allowed in chain
        assign_min_ident: Minimum identity (matches/aln_len) threshold
        assign_chain_topk: Number of top-scoring chains to consider per contig-ref pair
        assign_chain_score: Chain scoring method (matches/qbp_ident/matches_ident)
        assign_chain_min_bp: Minimum chain union-bp to keep
        assign_ref_score: Reference scoring method (topk/all)

    Returns:
        ChainEvidenceResult containing chains, segments, gene counts, and scoring information

    Note:
        In protein mode, miniprot PAF has qname=protein/transcript, tname=query contig.
        We map qname to reference chromosome coordinates via GFF3 transcript features,
        then build blocks in (query contig coordinate) vs (reference genomic coordinate) space.
        Overlapping hits are filtered by identity using interval trees (O(n log n)).
        Gate filtering requires ≥5 segments, ≥3 genes (configurable via CLI) because
        protein alignments produce many fragmented hits across a chromosome.
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
    blocks, n_before, n_removed, removed_per_contig = _filter_overlapping_hits_by_identity(blocks)
    if n_removed > 0:
        logger.info(f"Filtered {n_removed}/{n_before} overlapping protein hits by identity")
        # Show contigs with most filtered hits (helps identify potential issues)
        if removed_per_contig:
            top_contigs = sorted(removed_per_contig.items(), key=lambda x: -x[1])[:3]
            if top_contigs[0][1] > 10:  # Only show if significant filtering
                top_str = ", ".join(f"{c}:{n}" for c, n in top_contigs)
                logger.info(f"  Top filtered contigs: {top_str}")

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
        segments_strand_from_blocks=False,  # Don't use raw strand from blocks
        infer_strand_from_coords=True,  # Infer strand from query vs reference coordinate correlation
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
    infer_strand_from_coords: bool = False,
) -> ChainEvidenceResult:
    # Precompute contig-level strand inference if needed
    # This uses ALL blocks for each (contig, ref_id) pair to determine orientation
    contig_ref_strand: dict[tuple[str, str], str] = {}
    if infer_strand_from_coords:
        # Get unique (contig, ref_id) pairs
        contig_ref_pairs = set()
        for (q, ref_id, _strand) in blocks.keys():
            contig_ref_pairs.add((q, ref_id))
        # Compute strand for each pair using all blocks
        for (q, ref_id) in contig_ref_pairs:
            contig_ref_strand[(q, ref_id)] = _infer_contig_strand_from_all_blocks(
                blocks, q, ref_id
            )

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
        chain_blocks: list[list[Block]] = []  # Track blocks per chain for strand inference

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
                chain_blocks.append([b])  # Start new block list for this chain
            else:
                ch = chains[best_i]
                if b.gene_id:
                    ch.gene_ids.add(str(b.gene_id))
                ch.last_qe = max(ch.last_qe, b.qe)
                ch.last_r = rpos
                ch.q_intervals.append((b.qs, b.qe))
                ch.matches_sum += b.matches
                ch.alnlen_sum += b.aln_len
                chain_blocks[best_i].append(b)  # Add block to existing chain

        qlen = int(contig_lengths.get(q, 0) or qlens_from_paf.get(q, 0) or 0)
        chrom_id, sub = split_chrom_subgenome(ref_id)

        output_chain_id = 0
        for chain_idx, ch in enumerate(chains):
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

            output_chain_id += 1

            # Determine strand for this chain
            if segments_strand_from_blocks:
                seg_strand = strand
            elif infer_strand_from_coords:
                # Use precomputed contig-level strand (based on all blocks for this contig-ref pair)
                seg_strand = contig_ref_strand.get((q, ref_id), "+")
            else:
                seg_strand = "+"

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
                    output_chain_id,
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
                chain_segments_rows.append((q, qlen, chrom_id, sub, seg_strand, output_chain_id, s, e))

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
