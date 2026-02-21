#!/usr/bin/env python3
"""Pairwise assembly-vs-assembly synteny alignment for riparian visualization.

Runs minimap2 between adjacent assembly pairs using their oriented *.chrs.fasta
outputs, parses alignments into chains/macro_blocks, and writes pairwise TSVs
for the comparison report synteny plot.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from final_finalizer.alignment.chain_parsing import parse_paf_chain_evidence_and_segments
from final_finalizer.alignment.external_tools import run_minimap2_synteny
from final_finalizer.output.tsv_output import write_macro_blocks_tsv
from final_finalizer.utils.io_utils import file_exists_and_valid
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.sequence_utils import read_fasta_lengths

logger = get_logger("pairwise")


def compute_pairwise_synteny(
    left_fasta: Path,
    right_fasta: Path,
    left_name: str,
    right_name: str,
    outprefix: Path,
    threads: int,
    *,
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
    """Run minimap2 + chain parsing between two assemblies' chrs.fasta files.

    All parameters are plain types (Path/str/int/float/bool) so this function
    is fully serializable and can be submitted directly to executorlib.

    Args:
        left_fasta: Path to left (target) assembly chrs.fasta.
        right_fasta: Path to right (query) assembly chrs.fasta.
        left_name: Short name for left assembly.
        right_name: Short name for right assembly.
        outprefix: Output prefix (e.g., output_dir/pairwise/asm1_vs_asm2).
        threads: Number of threads for minimap2.
        preset: minimap2 preset (default "asm5").
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

    logger.info(f"Pairwise synteny: {left_name} vs {right_name}")

    # Output paths
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    paf_gz = Path(str(outprefix) + ".paf.gz")
    err_log = Path(str(outprefix) + ".alignment.err")
    macro_blocks_tsv = Path(str(outprefix) + ".macro_blocks.tsv")

    # Check if macro_blocks already exists (reuse cached result)
    if file_exists_and_valid(macro_blocks_tsv):
        logger.info(f"Pairwise macro_blocks exists, reusing: {macro_blocks_tsv}")
        return macro_blocks_tsv

    # Run minimap2 with permissive chaining (same as ref→query nucleotide mode)
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

    # Read right assembly chrs.fasta lengths (query in minimap2 terms)
    right_lengths = read_fasta_lengths(right_fasta)

    if not right_lengths:
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"no sequences in right FASTA"
        )
        return None

    # Parse PAF into chains and macro_blocks
    # ref_id_map=None: keep left assembly contig names as-is (they're already
    # the new_name values from classification)
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

    # Write pairwise macro_blocks TSV (no ref normalization needed)
    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=None)
    logger.info(
        f"Pairwise macro_blocks: {macro_blocks_tsv} "
        f"({len(ev.macro_block_rows)} blocks)"
    )

    return macro_blocks_tsv
