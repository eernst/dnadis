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
from final_finalizer.models import AssemblyResult
from final_finalizer.output.tsv_output import write_macro_blocks_tsv
from final_finalizer.utils.io_utils import file_exists_and_valid
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.sequence_utils import read_fasta_lengths

logger = get_logger("pairwise")


def compute_pairwise_synteny(
    left_result: AssemblyResult,
    right_result: AssemblyResult,
    outprefix: Path,
    threads: int,
    args,
) -> Optional[Path]:
    """Run minimap2 + chain parsing between two assemblies' chrs.fasta files.

    Uses the oriented *.chrs.fasta outputs (sequences already reverse-complemented
    and filtered to chrom_assigned contigs) for fast, clean pairwise alignment.

    Args:
        left_result: AssemblyResult for the left (target) assembly.
        right_result: AssemblyResult for the right (query) assembly.
        outprefix: Output prefix (e.g., output_dir/pairwise/asm1_vs_asm2).
        threads: Number of threads for minimap2.
        args: Parsed CLI arguments (for chain parsing parameters).

    Returns:
        Path to the pairwise macro_blocks TSV, or None if alignment could not
        be performed (e.g., missing chrs.fasta).
    """
    left_name = left_result.assembly_name
    right_name = right_result.assembly_name

    # Derive chrs.fasta paths from assembly outprefix
    left_chrs = Path(str(left_result.outprefix) + ".chrs.fasta")
    right_chrs = Path(str(right_result.outprefix) + ".chrs.fasta")

    if not file_exists_and_valid(left_chrs):
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"left chrs.fasta not found or empty: {left_chrs}"
        )
        return None
    if not file_exists_and_valid(right_chrs):
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"right chrs.fasta not found or empty: {right_chrs}"
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
        ref=left_chrs,
        qry=right_chrs,
        paf_gz_out=paf_gz,
        threads=threads,
        preset=args.preset,
        k=args.kmer,
        w=args.window,
        err_path=err_log,
        permissive=True,
    )

    # Read right assembly chrs.fasta lengths (query in minimap2 terms)
    right_lengths = read_fasta_lengths(right_chrs)

    if not right_lengths:
        logger.warning(
            f"Skipping pairwise {left_name} vs {right_name}: "
            f"no sequences in right chrs.fasta"
        )
        return None

    # Parse PAF into chains and macro_blocks
    # ref_id_map=None: keep left assembly contig names as-is (they're already
    # the new_name values from classification)
    ev = parse_paf_chain_evidence_and_segments(
        paf_gz_path=paf_gz,
        contig_lengths=right_lengths,
        assign_minlen=args.assign_minlen,
        assign_minmapq=args.assign_minmapq,
        assign_tp=args.assign_tp,
        chain_q_gap=args.chain_q_gap,
        chain_r_gap=args.chain_r_gap,
        chain_diag_slop=args.chain_diag_slop,
        assign_min_ident=args.assign_min_ident,
        assign_chain_topk=args.assign_chain_topk,
        assign_chain_score=args.assign_chain_score,
        assign_chain_min_bp=args.assign_chain_min_bp,
        assign_ref_score=args.assign_ref_score,
        ref_id_map=None,
    )

    # Write pairwise macro_blocks TSV (no ref normalization needed)
    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=None)
    logger.info(
        f"Pairwise macro_blocks: {macro_blocks_tsv} "
        f"({len(ev.macro_block_rows)} blocks)"
    )

    return macro_blocks_tsv
