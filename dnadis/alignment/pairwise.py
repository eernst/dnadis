#!/usr/bin/env python3
"""Pairwise assembly-vs-assembly synteny alignment for riparian visualization.

For nucleotide synteny mode, runs minimap2 between adjacent assembly pairs
using their oriented *.chrs.fasta outputs.  For protein synteny mode, aligns
reference proteins to each pair's chrs.fasta with miniprot, then derives
asm-vs-asm synteny anchors from proteins shared between the two PAFs — the
same protein-anchored synteny philosophy used by the per-assembly phase, so
ribbons remain meaningful for assemblies too divergent for direct nucleotide
alignment.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from dnadis.alignment.chain_parsing import (
    _chains_to_evidence_and_segments,
    _filter_overlapping_hits_by_identity,
    parse_paf_chain_evidence_and_segments,
)
from dnadis.alignment.external_tools import run_minimap2_synteny, run_miniprot
from dnadis.models import Block
from dnadis.output.tsv_output import write_macro_blocks_tsv
from dnadis.utils.io_utils import file_exists_and_valid, open_maybe_gzip
from dnadis.utils.logging_config import get_logger
from dnadis.utils.sequence_utils import read_fasta_lengths

logger = get_logger("pairwise")


def compute_pairwise_synteny(
    left_fasta: Path,
    right_fasta: Path,
    left_name: str,
    right_name: str,
    outprefix: Path,
    threads: int,
    *,
    synteny_mode: str = "nucleotide",
    proteins_faa: Optional[Path] = None,
    miniprot_exe: str = "miniprot",
    miniprot_args: str = "",
    preset: str = "asm5",
    kmer: Optional[int] = None,
    window: Optional[int] = None,
    assign_minlen: int = 5000,
    assign_minmapq: int = 0,
    assign_tp: str = "PI",
    chain_q_gap: int = 500_000,
    chain_r_gap: int = 500_000,
    chain_diag_slop: int = 200_000,
    assign_min_ident: float = 0.0,
    assign_chain_topk: int = 3,
    assign_chain_score: str = "matches",
    assign_chain_min_bp: int = 50_000,
    assign_ref_score: str = "all",
) -> Optional[Path]:
    """Compute pairwise synteny between two assemblies' chrs.fasta files.

    Dispatches on ``synteny_mode``:
    - ``nucleotide``: run minimap2 directly on the two FASTAs.
    - ``protein``: align reference proteins to each FASTA with miniprot and
      derive asm-vs-asm anchors from proteins shared between the two PAFs.
      Requires ``proteins_faa``.

    All parameters are plain types (Path/str/int/float/bool) so this function
    is fully serializable and can be submitted directly to executorlib.

    Args:
        left_fasta: Path to left (target) assembly chrs.fasta.
        right_fasta: Path to right (query) assembly chrs.fasta.
        left_name: Short name for left assembly.
        right_name: Short name for right assembly.
        outprefix: Output prefix (e.g., output_dir/pairwise/asm1_vs_asm2).
        threads: Number of threads for minimap2/miniprot.
        synteny_mode: "nucleotide" or "protein".
        proteins_faa: Reference proteins FASTA (required when synteny_mode is
            "protein"; ignored otherwise).
        miniprot_exe: miniprot executable name (protein mode only).
        miniprot_args: extra args for miniprot (protein mode only).
        preset: minimap2 preset (default "asm5"; nucleotide mode only).
        kmer: minimap2 k-mer size (None = use preset default).
        window: minimap2 window size (None = use preset default).
        assign_minlen: Minimum alignment length for chain parsing.
        assign_minmapq: Minimum mapping quality.
        assign_tp: Which alignment type tags to accept ("P", "PI", or "ALL").
        chain_q_gap: Maximum query gap for chaining.
        chain_r_gap: Maximum reference gap for chaining.
        chain_diag_slop: Diagonal slop for chaining.
        assign_min_ident: Minimum alignment identity.
        assign_chain_topk: Top-K chains to consider per contig.
        assign_chain_score: Chain scoring method.
        assign_chain_min_bp: Minimum chain span in bp.
        assign_ref_score: Reference scoring method.

    Returns:
        Path to the pairwise macro_blocks TSV, or None if alignment could not
        be performed (e.g., missing chrs.fasta).
    """
    if not file_exists_and_valid(left_fasta):
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"left FASTA not found or empty: {left_fasta}"
        )
        return None
    if not file_exists_and_valid(right_fasta):
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"right FASTA not found or empty: {right_fasta}"
        )
        return None

    outprefix.parent.mkdir(parents=True, exist_ok=True)
    macro_blocks_tsv = Path(str(outprefix) + ".macro_blocks.tsv")
    if file_exists_and_valid(macro_blocks_tsv):
        # Mode-aware cache validation.  Only invalidate when the
        # intermediate file from the *other* synteny_mode is present and the
        # one from the requested mode is missing — that's positive evidence
        # the cached TSV came from the wrong mode.  If neither intermediate
        # is present (e.g. a hand-staged TSV), preserve the prior behavior
        # and reuse the cache.
        nuc_paf = Path(str(outprefix) + ".paf.gz")
        prot_left = Path(str(outprefix) + ".left.miniprot.paf.gz")
        prot_right = Path(str(outprefix) + ".right.miniprot.paf.gz")
        nuc_intermediate = file_exists_and_valid(nuc_paf)
        prot_intermediate = (
            file_exists_and_valid(prot_left)
            and file_exists_and_valid(prot_right)
        )

        if synteny_mode == "protein":
            wrong_mode_cache = nuc_intermediate and not prot_intermediate
        else:
            wrong_mode_cache = prot_intermediate and not nuc_intermediate

        if wrong_mode_cache:
            logger.info(
                f"Pairwise {left_name} vs {right_name}: cached macro_blocks "
                f"was produced by a different --synteny-mode; recomputing "
                f"for {synteny_mode} mode"
            )
            macro_blocks_tsv.unlink()
        else:
            logger.info(f"Pairwise macro_blocks exists, reusing: {macro_blocks_tsv}")
            return macro_blocks_tsv

    if synteny_mode == "protein":
        if proteins_faa is None or not file_exists_and_valid(proteins_faa):
            logger.warning(
                f"Skipping pairwise {left_name} vs {right_name} (protein mode): "
                f"proteins FASTA not found or empty: {proteins_faa}"
            )
            return None
        return _compute_pairwise_synteny_protein(
            left_fasta=left_fasta,
            right_fasta=right_fasta,
            left_name=left_name,
            right_name=right_name,
            proteins_faa=proteins_faa,
            miniprot_exe=miniprot_exe,
            miniprot_args=miniprot_args,
            outprefix=outprefix,
            threads=threads,
            assign_minlen=assign_minlen,
            assign_minmapq=assign_minmapq,
            chain_q_gap=chain_q_gap,
            chain_r_gap=chain_r_gap,
            chain_diag_slop=chain_diag_slop,
            assign_min_ident=assign_min_ident,
            assign_chain_topk=assign_chain_topk,
            assign_chain_score=assign_chain_score,
            assign_chain_min_bp=assign_chain_min_bp,
            assign_ref_score=assign_ref_score,
        )

    return _compute_pairwise_synteny_nucleotide(
        left_fasta=left_fasta,
        right_fasta=right_fasta,
        left_name=left_name,
        right_name=right_name,
        outprefix=outprefix,
        threads=threads,
        preset=preset,
        kmer=kmer,
        window=window,
        assign_minlen=assign_minlen,
        assign_minmapq=assign_minmapq,
        assign_tp=assign_tp,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
    )


def _compute_pairwise_synteny_nucleotide(
    *,
    left_fasta: Path,
    right_fasta: Path,
    left_name: str,
    right_name: str,
    outprefix: Path,
    threads: int,
    preset: str,
    kmer: Optional[int],
    window: Optional[int],
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
) -> Optional[Path]:
    """Pairwise synteny via direct minimap2 alignment between chrs.fasta files."""
    logger.info(f"Pairwise synteny (nucleotide): {left_name} vs {right_name}")

    paf_gz = Path(str(outprefix) + ".paf.gz")
    err_log = Path(str(outprefix) + ".alignment.err")
    macro_blocks_tsv = Path(str(outprefix) + ".macro_blocks.tsv")

    run_minimap2_synteny(
        ref=left_fasta,
        qry=right_fasta,
        paf_gz_out=paf_gz,
        threads=threads,
        preset=preset,
        k=kmer,
        w=window,
        err_path=err_log,
        permissive=True,
    )

    right_lengths = read_fasta_lengths(right_fasta)
    if not right_lengths:
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"no sequences in right FASTA"
        )
        return None

    # ref_id_map=None: left assembly contig names are already the renamed
    # (post-classification) names from chrs.fasta.
    ev = parse_paf_chain_evidence_and_segments(
        paf_gz_path=paf_gz,
        contig_lengths=right_lengths,
        assign_minlen=assign_minlen,
        assign_minmapq=assign_minmapq,
        assign_tp=assign_tp,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
        ref_id_map=None,
    )

    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=None)
    logger.info(
        f"Pairwise macro_blocks: {macro_blocks_tsv} "
        f"({len(ev.macro_block_rows)} blocks)"
    )
    return macro_blocks_tsv


# ---------------------------------------------------------------------------
# Protein-anchored pairwise synteny
# ---------------------------------------------------------------------------

# A single miniprot hit (protein → asm contig).
_ProteinHit = Tuple[str, int, int, int, str, int, int, int]
# (contig, contig_len, t_start, t_end, strand, matches, aln_len, mapq)


def _read_miniprot_protein_hits(
    paf_gz: Path,
) -> Dict[str, List[_ProteinHit]]:
    """Parse a miniprot PAF.gz into a dict of {protein_id: [hit, ...]}.

    Each hit captures the asm contig, its length, the protein-to-genome
    target span, the strand reported by miniprot, and the alignment quality
    metrics (matches/aln_len/mapq).
    """
    hits: Dict[str, List[_ProteinHit]] = defaultdict(list)
    with open_maybe_gzip(paf_gz, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            prot = fields[0].replace("mRNA:", "").replace("transcript:", "")
            qstrand = fields[4]
            tname = fields[5]

            try:
                tlen = int(fields[6])
                ts = int(fields[7])
                te = int(fields[8])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue

            if te < ts:
                ts, te = te, ts
            if te <= ts:
                continue

            hits[prot].append((tname, tlen, ts, te, qstrand, matches, aln_len, mapq))
    return hits


def _compute_pairwise_synteny_protein(
    *,
    left_fasta: Path,
    right_fasta: Path,
    left_name: str,
    right_name: str,
    proteins_faa: Path,
    miniprot_exe: str,
    miniprot_args: str,
    outprefix: Path,
    threads: int,
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
) -> Optional[Path]:
    """Pairwise synteny derived from shared reference protein anchors.

    Aligns ``proteins_faa`` against each pair's chrs.fasta with miniprot,
    then for every protein found in both PAFs emits a synteny "block" linking
    the two assemblies' contigs.  The resulting blocks are de-duplicated by
    identity and chained with the same logic as the per-assembly miniprot
    path, producing macro_blocks in the standard pairwise schema (left
    contigs in the ``ref_id`` column, right contigs in the ``contig``
    column).

    The ``assign_minlen`` gate is interpreted as the minimum target span (in
    bp on the contig genome) for an individual protein anchor block, so a
    smaller value than the nucleotide path's 5 kb default is appropriate
    for protein anchors (typically 100 bp – a few kb).
    """
    logger.info(f"Pairwise synteny (protein): {left_name} vs {right_name}")

    left_paf = Path(str(outprefix) + ".left.miniprot.paf.gz")
    right_paf = Path(str(outprefix) + ".right.miniprot.paf.gz")
    macro_blocks_tsv = Path(str(outprefix) + ".macro_blocks.tsv")
    left_err = Path(str(outprefix) + ".left.miniprot.err")
    right_err = Path(str(outprefix) + ".right.miniprot.err")

    run_miniprot(
        miniprot_exe, left_fasta, proteins_faa, left_paf,
        threads, miniprot_args, left_err,
    )
    run_miniprot(
        miniprot_exe, right_fasta, proteins_faa, right_paf,
        threads, miniprot_args, right_err,
    )

    left_hits = _read_miniprot_protein_hits(left_paf)
    right_hits = _read_miniprot_protein_hits(right_paf)
    shared = set(left_hits.keys()) & set(right_hits.keys())
    if not shared:
        logger.warning(
            f"Pairwise {left_name} vs {right_name} (protein): no shared "
            f"reference proteins between assemblies; no anchors emitted"
        )
        write_macro_blocks_tsv(macro_blocks_tsv, [], ref_norm_to_orig=None)
        return macro_blocks_tsv

    # Build blocks for the chain parser.  Treat the LEFT assembly as the
    # reference side (ref_id = left contig) and RIGHT as the query side
    # (q = right contig), matching the nucleotide pairwise convention.
    blocks: Dict[Tuple[str, str, str], list] = defaultdict(list)
    qlens_from_paf: Dict[str, int] = {}
    rlens_from_paf: Dict[str, int] = {}

    for prot_id in shared:
        for L in left_hits[prot_id]:
            l_contig, l_len, l_ts, l_te, l_strand, l_m, l_a, l_mapq = L
            for R in right_hits[prot_id]:
                r_contig, r_len, r_ts, r_te, r_strand, r_m, r_a, r_mapq = R

                if (r_te - r_ts) < assign_minlen:
                    continue
                if min(l_mapq, r_mapq) < assign_minmapq:
                    continue

                # Conservative pairwise quality: take the weaker of the two
                # individual hits, and combine identities multiplicatively
                # (left-to-protein × protein-to-right approximation).
                matches = min(l_m, r_m)
                aln_len = min(l_a, r_a)
                ident = (matches / aln_len) if aln_len > 0 else 0.0
                if ident < assign_min_ident:
                    continue

                # Chain-key strand kept '+' so chain logic works in
                # increasing reference-coordinate space.  The actual orient
                # is inferred from coords later (infer_strand_from_coords).
                blocks[(r_contig, l_contig, "+")].append(
                    Block(
                        qs=r_ts, qe=r_te,
                        rs=l_ts, re=l_te,
                        matches=matches,
                        aln_len=aln_len,
                        mapq=min(l_mapq, r_mapq),
                        strand="+",
                        gene_id=prot_id,
                    )
                )
                qlens_from_paf.setdefault(r_contig, r_len)
                rlens_from_paf.setdefault(l_contig, l_len)

    if not blocks:
        logger.warning(
            f"Pairwise {left_name} vs {right_name} (protein): no anchor "
            f"blocks survived gating (assign_minlen={assign_minlen}, "
            f"assign_min_ident={assign_min_ident})"
        )
        write_macro_blocks_tsv(macro_blocks_tsv, [], ref_norm_to_orig=None)
        return macro_blocks_tsv

    blocks, n_before, n_removed, _removed_per_contig = (
        _filter_overlapping_hits_by_identity(blocks)
    )
    if n_removed > 0:
        logger.info(
            f"Pairwise {left_name} vs {right_name} (protein): filtered "
            f"{n_removed}/{n_before} overlapping anchors by identity"
        )

    right_lengths = read_fasta_lengths(right_fasta)
    if not right_lengths:
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"no sequences in right FASTA"
        )
        return None

    ev = _chains_to_evidence_and_segments(
        blocks=blocks,
        contig_lengths=right_lengths,
        qlens_from_paf=qlens_from_paf,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
        segments_strand_from_blocks=False,
        infer_strand_from_coords=True,
        rlens_from_paf=rlens_from_paf,
    )

    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=None)
    logger.info(
        f"Pairwise macro_blocks: {macro_blocks_tsv} "
        f"({len(ev.macro_block_rows)} blocks from {len(shared)} shared proteins)"
    )
    return macro_blocks_tsv
