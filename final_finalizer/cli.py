#!/usr/bin/env python3
"""Command-line interface for final_finalizer."""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set, Tuple


def _validate_input_path(path_str: str) -> Path:
    """Validate that an input file path exists and is readable.

    Used as argparse type for input file arguments.
    """
    path = Path(path_str).resolve()
    if not path.exists():
        raise argparse.ArgumentTypeError(f"Input file not found: {path}")
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Input path is not a file: {path}")
    return path


def _positive_int(value: str) -> int:
    """Validate that a value is a positive integer.

    Used as argparse type for thread count and similar arguments.
    """
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: {value}")
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"Value must be positive, got {value}")
    return ivalue

from final_finalizer.utils.io_utils import file_exists_and_valid
from final_finalizer.utils.sequence_utils import (
    read_fasta_lengths,
    write_filtered_fasta,
)
from final_finalizer.utils.reference_utils import (
    compile_ref_id_patterns,
    filter_gff3_by_ref,
    get_min_nuclear_chrom_length,
    parse_gff3_transcript_coords,
    read_fasta_lengths_with_map,
    set_ref_id_patterns,
    write_ref_lengths_tsv,
)
from final_finalizer.alignment.external_tools import (
    run_gffread_extract_proteins,
    run_minimap2_synteny,
    run_miniprot,
)
from final_finalizer.alignment.stats import (
    parse_paf_primary,
    write_stats_tsv,
)
from final_finalizer.alignment.chain_parsing import (
    parse_miniprot_synteny_evidence_and_segments,
)
from final_finalizer.detection.organelle import (
    prepare_organelle_references,
    detect_organelles,
)
from final_finalizer.detection.rdna import (
    prepare_rdna_reference,
    detect_rdna_contigs,
)
from final_finalizer.detection.debris import detect_chromosome_debris
from final_finalizer.detection.contaminant import detect_contaminants
from final_finalizer.classification.classifier import (
    classify_all_contigs,
    classify_debris_and_unclassified,
    compute_mean_genes_per_Mbp,
    count_genes_per_ref_chrom,
    determine_contig_orientations,
)
from final_finalizer.output.fasta_output import write_classified_fastas
from final_finalizer.output.tsv_output import (
    build_segment_support_from_rows,
    write_chain_segments_tsv,
    write_chain_summary_tsv,
    write_contig_summary_tsv,
    write_macro_blocks_tsv,
)
from final_finalizer.output.plotting import run_plot


def main():
    p = argparse.ArgumentParser(
        description="Genome-wide subgenome distance analysis using minimap2 chains or miniprot protein-anchored synteny blocks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # =========================================================================
    # Required arguments
    # =========================================================================
    req = p.add_argument_group("Required arguments")
    req.add_argument("-r", "--ref", type=_validate_input_path, required=True, help="Reference assembly FASTA (plain or gzipped; chrNA/chrNP/etc.)")
    req.add_argument("-q", "--query", type=_validate_input_path, required=True, help="Query assembly FASTA (plain or gzipped; contigs/scaffolds)")
    req.add_argument("-o", "--outprefix", required=True, help="Output prefix (no extension)")
    req.add_argument(
        "--ref-gff3",
        type=_validate_input_path,
        required=True,
        help="Reference GFF3 with protein-coding gene annotations. "
        "Used to extract proteins (gffread) and run protein-anchor synteny blocks (miniprot).",
    )

    # =========================================================================
    # Common options
    # =========================================================================
    common = p.add_argument_group("Common options")
    common.add_argument("-t", "--threads", type=_positive_int, default=8, help="Threads for minimap2/miniprot [8]")
    common.add_argument("--plot", action="store_true", help="Generate overview plots with R/ggplot2")
    common.add_argument("--plot-html", action="store_true", help="Also generate interactive HTML plot (ggiraph)")
    common.add_argument(
        "-C", "--chr-like-minlen", type=int, default=None,
        help="Minimum contig length (bp) to be considered chromosome-like. Default: 80% of smallest nuclear ref chromosome.",
    )
    common.add_argument(
        "--add-subgenome-suffix", type=str, default=None,
        help="Add this suffix to contig names when reference has no subgenome identifiers (e.g., 'A' -> chr5A)",
    )
    common.add_argument(
        "--ref-id-pattern",
        action="append",
        default=None,
        help=(
            "Regex for reference IDs with named group 'chrom' and optional 'sg'. "
            "May be provided multiple times; tried in order. "
            "Example: --ref-id-pattern '^(?P<chrom>chr\\d+)(?P<sg>[AB])$' for chr1A/chr1B style. "
            "If unset, presets handle chrN[A-Z], chrN(At|Dt), and chrN patterns."
        ),
    )

    # =========================================================================
    # Pipeline phase toggles
    # =========================================================================
    toggles = p.add_argument_group("Pipeline phase toggles")
    toggles.add_argument(
        "--nt-synteny", action="store_true",
        help="Run nucleotide-based synteny alignments (minimap2/mm2plus) for QA/diagnostics. Optional; classification uses protein-anchored synteny only.",
    )
    toggles.add_argument("--skip-organelles", action="store_true", help="Skip organelle detection.")
    toggles.add_argument("--skip-rdna", action="store_true", help="Skip rDNA detection.")
    toggles.add_argument("--skip-contaminants", action="store_true", help="Skip contaminant detection.")

    # =========================================================================
    # External tool paths
    # =========================================================================
    tools = p.add_argument_group("External tool paths")
    tools.add_argument("--gffread", default="gffread", help="gffread executable [gffread]")
    tools.add_argument("--miniprot", default="miniprot", help="miniprot executable [miniprot]")
    tools.add_argument("--miniprot-args", default="", help='Extra args for miniprot (e.g. "-G 200k" to increase max intron size for plants with large introns) [none]')

    # =========================================================================
    # Reference inputs for classification
    # =========================================================================
    refs = p.add_argument_group("Reference inputs for classification")
    refs.add_argument(
        "--chrC-ref", type=str, default=None,
        help="Chloroplast reference FASTA. Default: extract chrC from --ref if present.",
    )
    refs.add_argument(
        "--chrM-ref", type=str, default=None,
        help="Mitochondrion reference FASTA. Default: extract chrM from --ref if present.",
    )
    refs.add_argument(
        "--rdna-ref", type=str, default=None,
        help="rDNA reference FASTA, or 'default' for bundled Arabidopsis 45S reference.",
    )
    refs.add_argument(
        "--centrifuger-idx", type=str, default=None,
        help="Path prefix to centrifuger index for contaminant screening.",
    )

    # =========================================================================
    # Thresholds: Chromosome assignment
    # =========================================================================
    chr_thresh = p.add_argument_group("Thresholds: Chromosome assignment")
    chr_thresh.add_argument(
        "--assign-min-frac", type=float, default=0.10,
        help="Min union-bp/contig_len to accept assignment [0.10]",
    )
    chr_thresh.add_argument(
        "--assign-min-ratio", type=float, default=1.25,
        help="Min best/second score ratio for unambiguous assignment [1.25]",
    )
    chr_thresh.add_argument(
        "-P", "--chimera-primary-frac", type=float, default=0.8,
        help="Primary fraction threshold for chimera detection [0.8]",
    )
    chr_thresh.add_argument(
        "-S", "--chimera-secondary-frac", type=float, default=0.2,
        help="Secondary fraction threshold for chimera detection [0.2]",
    )
    chr_thresh.add_argument(
        "--low-ref-span-threshold", type=float, default=0.75,
        help="Fraction of ref chrom that must be spanned by syntenic blocks [0.75]",
    )

    # =========================================================================
    # Thresholds: Synteny block building
    # =========================================================================
    block_thresh = p.add_argument_group("Thresholds: Synteny block building")
    block_thresh.add_argument(
        "--assign-minlen", type=int, default=10_000,
        help="Min alignment length for block building [10000]",
    )
    block_thresh.add_argument(
        "--assign-minmapq", type=int, default=0,
        help="Min MAPQ for block building [0]",
    )
    block_thresh.add_argument(
        "--assign-min-ident", type=float, default=0.0,
        help="Min identity proxy (matches/aln_len) [0.0]",
    )
    block_thresh.add_argument(
        "--assign-tp", choices=["P", "PI", "ALL"], default="PI",
        help="(nt-synteny) Which tp tags: P=primary, PI=primary+supp, ALL=all [PI]",
    )

    # =========================================================================
    # Thresholds: Chain scoring
    # =========================================================================
    chain_thresh = p.add_argument_group("Thresholds: Chain scoring")
    chain_thresh.add_argument(
        "--chain-q-gap", type=int, default=200_000,
        help="Max query gap between blocks in chain [200000]",
    )
    chain_thresh.add_argument(
        "--chain-r-gap", type=int, default=400_000,
        help="Max ref gap between blocks in chain [400000]",
    )
    chain_thresh.add_argument(
        "--chain-diag-slop", type=int, default=150_000,
        help="Max diagonal drift in chain [150000]",
    )
    chain_thresh.add_argument(
        "--assign-chain-min-bp", type=int, default=0,
        help="Min chain union-bp to keep [0]",
    )
    chain_thresh.add_argument(
        "--assign-chain-score", choices=["matches", "qbp_ident", "matches_ident"], default="matches",
        help="How to score each chain [matches]",
    )
    chain_thresh.add_argument(
        "--assign-chain-topk", type=int, default=3,
        help="Sum top-K chain scores per contig/ref [3]",
    )
    chain_thresh.add_argument(
        "--assign-ref-score", choices=["topk", "all"], default="all",
        help="Score per ref: topk or all kept chains [all]",
    )

    # =========================================================================
    # Thresholds: Protein-anchor (miniprot) gates
    # =========================================================================
    prot_thresh = p.add_argument_group("Thresholds: Protein-anchor (miniprot) gates")
    prot_thresh.add_argument(
        "--miniprot-min-genes", type=int, default=3,
        help="Min unique genes for best ref assignment [3]",
    )
    prot_thresh.add_argument(
        "--miniprot-min-segments", type=int, default=5,
        help="Min number of synteny segments [5]",
    )
    prot_thresh.add_argument(
        "--miniprot-min-span-frac", type=float, default=0.20,
        help="Min span fraction of contig [0.20]",
    )
    prot_thresh.add_argument(
        "--miniprot-min-span-bp", type=int, default=50_000,
        help="Min absolute span in bp [50000]",
    )

    # =========================================================================
    # Thresholds: Organelle detection
    # =========================================================================
    org_thresh = p.add_argument_group("Thresholds: Organelle detection")
    org_thresh.add_argument(
        "--organelle-min-cov", type=float, default=0.80,
        help="Min query coverage for organelle candidate [0.80]",
    )
    org_thresh.add_argument(
        "--chrC-len-tolerance", type=float, default=0.05,
        help="Length tolerance for chrC candidates (fraction) [0.05]",
    )
    org_thresh.add_argument(
        "--chrM-len-tolerance", type=float, default=0.20,
        help="Length tolerance for chrM candidates (fraction) [0.20]",
    )

    # =========================================================================
    # Thresholds: rDNA detection
    # =========================================================================
    rdna_thresh = p.add_argument_group("Thresholds: rDNA detection")
    rdna_thresh.add_argument(
        "--rdna-min-cov", type=float, default=0.50,
        help="Min query coverage for rDNA classification [0.50]",
    )

    # =========================================================================
    # Thresholds: Chromosome debris detection
    # =========================================================================
    chr_debris_thresh = p.add_argument_group("Thresholds: Chromosome debris detection")
    chr_debris_thresh.add_argument(
        "--chr-debris-min-cov", type=float, default=0.80,
        help="Min query coverage for detecting duplicates of assembled chromosomes [0.80]",
    )
    chr_debris_thresh.add_argument(
        "--chr-debris-min-identity", type=float, default=0.90,
        help="Min alignment identity vs assembled chromosomes [0.90]",
    )

    # =========================================================================
    # Thresholds: Contaminant detection
    # =========================================================================
    contam_thresh = p.add_argument_group("Thresholds: Contaminant detection")
    contam_thresh.add_argument(
        "--contaminant-min-score", type=int, default=150,
        help="Min centrifuger score for contaminant classification [150]",
    )

    # =========================================================================
    # Thresholds: Reference-based debris detection
    # =========================================================================
    debris_thresh = p.add_argument_group("Thresholds: Reference-based debris detection")
    debris_thresh.add_argument(
        "--debris-min-cov", type=float, default=0.50,
        help="Min alignment coverage vs reference for reference-based debris detection [0.50]",
    )
    debris_thresh.add_argument(
        "--debris-min-protein-hits", type=int, default=2,
        help="Min miniprot protein hits for debris [2]",
    )

    # =========================================================================
    # Minimap2 tuning (for --nt-synteny mode)
    # =========================================================================
    mm2_tune = p.add_argument_group("Minimap2 tuning (for --nt-synteny mode)")
    mm2_tune.add_argument("-x", "--preset", default="asm20", help="minimap2 preset [asm20]")
    mm2_tune.add_argument("-k", "--kmer", type=int, default=None, help="minimap2 k-mer size")
    mm2_tune.add_argument("-w", "--window", type=int, default=None, help="minimap2 minimizer window size")
    mm2_tune.add_argument(
        "-L", "--aln-minlen", type=int, default=10_000,
        help="Min alignment length in PRIMARY-only QA stats [10000]",
    )

    args = p.parse_args()

    if args.ref_id_pattern:
        set_ref_id_patterns(compile_ref_id_patterns(args.ref_id_pattern))

    # --- Phase 1: Reading reference and query FASTA ---
    print("[info] Phase 1: Reading reference and query FASTA", file=sys.stderr)
    ref = args.ref  # Already a Path from _validate_input_path
    qry = args.query  # Already a Path from _validate_input_path
    outprefix = Path(args.outprefix)
    ref_lengths_norm, ref_orig_to_norm, ref_norm_to_orig = read_fasta_lengths_with_map(ref)
    ref_ids_raw = set(read_fasta_lengths(ref).keys())
    qry_lengths = read_fasta_lengths(qry)
    print(f"[info] Query contigs: {len(qry_lengths)}", file=sys.stderr)

    # Outputs (shared)
    ref_lengths_tsv = Path(str(outprefix) + ".ref_lengths.tsv")
    segments_tsv = Path(str(outprefix) + ".segments.tsv")
    chain_summary_tsv = Path(str(outprefix) + ".evidence_summary.tsv")
    macro_blocks_tsv = Path(str(outprefix) + ".macro_blocks.tsv")

    # minimap2 QA outputs
    paf_gz = Path(str(outprefix) + ".all_vs_ref.paf.gz")
    err_log = Path(str(outprefix) + ".alignment.err")
    stats_tsv = Path(str(outprefix) + ".contig_ref.stats.tsv")

    # miniprot outputs
    proteins_faa = Path(str(outprefix) + ".ref_proteins.fa")
    gffread_err = Path(str(outprefix) + ".gffread.err")
    miniprot_paf_gz = Path(str(outprefix) + ".miniprot.paf.gz")
    miniprot_err = Path(str(outprefix) + ".miniprot.err")

    if not file_exists_and_valid(ref_lengths_tsv):
        print(f"[info] Writing reference lengths: {ref_lengths_tsv}", file=sys.stderr)
        write_ref_lengths_tsv(ref_lengths_tsv, ref)
    else:
        print(f"[info] Reference lengths TSV exists, reusing: {ref_lengths_tsv}", file=sys.stderr)

    # --- Phase 2: Nucleotide synteny analysis (optional QA) ---
    print("[info] Phase 2: Nucleotide synteny analysis (optional QA)", file=sys.stderr)
    if args.nt_synteny:
        run_minimap2_synteny(ref, qry, paf_gz, args.threads, args.preset, args.kmer, args.window, err_log)
        (
            aln_sum_p,
            match_sum_p,
            _alnlen_sum_p,
            _ctotal_p,
            _crefs_p,
            _cqlen_p,
            _best_aln_p,
            _best_ref_p,
            _second_aln_p,
            _second_ref_p,
        ) = parse_paf_primary(paf_gz, args.aln_minlen)
        write_stats_tsv(stats_tsv, aln_sum_p, match_sum_p)
        print(f"[done] Stats (primary QA): {stats_tsv}", file=sys.stderr)
    else:
        print("[info] --nt-synteny not set: skipping nucleotide alignments + primary-only QA stats.", file=sys.stderr)

    # --- Phase 3: Protein-anchor synteny analysis ---
    print("[info] Phase 3: Protein-anchor synteny analysis", file=sys.stderr)
    ref_gff3 = args.ref_gff3  # Already a Path from _validate_input_path
    filtered_ref_gff3 = Path(str(outprefix) + ".ref_gff3.filtered.gff3")
    ref_gff3 = filter_gff3_by_ref(
        ref_gff3,
        ref_ids_raw,
        filtered_ref_gff3,
        ref_norm_to_orig=ref_norm_to_orig,
    )
    print(f"[info] Protein-anchor mode enabled (ref GFF3: {ref_gff3})", file=sys.stderr)

    tx2loc, tx2gene = parse_gff3_transcript_coords(ref_gff3, ref_id_map=ref_orig_to_norm)
    if not (tx2loc and tx2gene):
        raise RuntimeError(f"No transcript features found in GFF3 (mRNA/transcript with ID=...). File: {ref_gff3}")

    run_gffread_extract_proteins(args.gffread, ref, ref_gff3, proteins_faa, gffread_err)
    run_miniprot(args.miniprot, qry, proteins_faa, miniprot_paf_gz, args.threads, args.miniprot_args, miniprot_err)

    ev = parse_miniprot_synteny_evidence_and_segments(
        miniprot_paf_gz=miniprot_paf_gz,
        contig_lengths=qry_lengths,
        tx2loc=tx2loc,
        tx2gene=tx2gene,
        assign_minlen=args.assign_minlen,
        assign_minmapq=args.assign_minmapq,
        chain_q_gap=args.chain_q_gap,
        chain_r_gap=args.chain_r_gap,
        chain_diag_slop=args.chain_diag_slop,
        assign_min_ident=args.assign_min_ident,
        assign_chain_topk=args.assign_chain_topk,
        assign_chain_score=args.assign_chain_score,
        assign_chain_min_bp=args.assign_chain_min_bp,
        assign_ref_score=args.assign_ref_score,
    )

    # --- Protein-anchor synteny gates: filter weak assignments ---
    seg_count, span_bp = build_segment_support_from_rows(ev.chain_segments_rows)

    def _passes_gate(
        q: str,
        ref_id: str,
        clen: int,
        min_segments: int,
        min_span_bp: int,
        min_span_frac: float,
        min_genes: int,
    ) -> tuple[bool, int, int, float, int]:
        key = (q, ref_id)
        nseg = int(seg_count.get(key, 0) or 0)
        spbp = int(span_bp.get(key, 0) or 0)
        spfrac = (spbp / clen) if clen > 0 else 0.0
        ngen = int(ev.qr_gene_count.get(key, 0) or 0)

        ok = (nseg >= min_segments) and (spbp >= min_span_bp) and (spfrac >= min_span_frac) and (ngen >= min_genes)
        return ok, nseg, spbp, spfrac, ngen

    # Gate-aware reranking over candidate refs
    # Make mutable copies for reranking
    best_ref = dict(ev.best_ref)
    best_score = dict(ev.best_score)
    best_bp = dict(ev.best_bp)
    second_ref = dict(ev.second_ref)
    second_score = dict(ev.second_score)
    second_bp = dict(ev.second_bp)

    n_demoted = 0
    n_switched = 0
    qr_ref_score = ev.qr_weight_all if args.assign_ref_score == "all" else ev.qr_score_topk

    candidates_by_contig = defaultdict(list)
    for (q, ref_id), sc in qr_ref_score.items():
        if sc > 0:
            candidates_by_contig[q].append((float(sc), ref_id))

    print(f"[info] Assignment score per ref: {args.assign_ref_score}", file=sys.stderr)

    for q in qry_lengths.keys():
        clen = int(qry_lengths.get(q, 0) or ev.qlens_from_paf.get(q, 0) or 0)
        if clen <= 0:
            continue

        cands = candidates_by_contig.get(q, [])
        if not cands:
            best_ref[q] = ""
            best_score[q] = 0.0
            best_bp[q] = 0
            second_ref[q] = ""
            second_score[q] = 0.0
            second_bp[q] = 0
            continue

        cands.sort(key=lambda x: (x[0], int(ev.qr_union_bp.get((q, x[1]), 0) or 0)), reverse=True)
        orig_best = cands[0][1]

        passing = []
        for sc, ref_id in cands:
            ok, _nseg, _spbp, _spfrac, _ngen = _passes_gate(
                q,
                ref_id,
                clen,
                args.miniprot_min_segments,
                args.miniprot_min_span_bp,
                args.miniprot_min_span_frac,
                args.miniprot_min_genes,
            )
            if ok:
                passing.append((sc, ref_id))
            if len(passing) >= 2:
                break

        if not passing:
            best_ref[q] = ""
            best_score[q] = 0.0
            best_bp[q] = 0
            second_ref[q] = ""
            second_score[q] = 0.0
            second_bp[q] = 0
            n_demoted += 1
            continue

        chosen_score, chosen_ref_id = passing[0]
        best_ref[q] = chosen_ref_id
        best_score[q] = float(chosen_score)
        best_bp[q] = int(ev.qr_union_bp.get((q, chosen_ref_id), 0) or 0)

        if len(passing) > 1:
            second_score_val, second_ref_id = passing[1]
            second_ref[q] = second_ref_id
            second_score[q] = float(second_score_val)
            second_bp[q] = int(ev.qr_union_bp.get((q, second_ref_id), 0) or 0)
        else:
            second_ref[q] = ""
            second_score[q] = 0.0
            second_bp[q] = 0

        if chosen_ref_id != orig_best:
            n_switched += 1

    print(f"[info] Protein-anchor gate-aware rerank: demoted={n_demoted} switched={n_switched}", file=sys.stderr)

    plot_suffix = "protein-anchor synteny"
    print(f"[done] Proteins:         {proteins_faa}", file=sys.stderr)
    print(f"[done] miniprot PAF.gz:  {miniprot_paf_gz}", file=sys.stderr)
    print(f"[done] miniprot err:     {miniprot_err}", file=sys.stderr)

    # Write segments + evidence summary TSVs
    write_chain_segments_tsv(segments_tsv, ev.chain_segments_rows, ref_norm_to_orig=ref_norm_to_orig)
    write_chain_summary_tsv(chain_summary_tsv, ev.chain_summary_rows, ref_norm_to_orig=ref_norm_to_orig)
    print(f"[done] Segments TSV:      {segments_tsv}", file=sys.stderr)
    print(f"[done] Evidence TSV:      {chain_summary_tsv}", file=sys.stderr)

    # Write macro blocks TSV
    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=ref_norm_to_orig)
    print(f"[done] Macro blocks TSV:  {macro_blocks_tsv}", file=sys.stderr)

    # --- Compute chr_like_minlen if not provided ---
    ref_lengths = ref_lengths_norm
    if args.chr_like_minlen is None:
        min_nuclear = get_min_nuclear_chrom_length(ref_lengths)
        chr_like_minlen = int(min_nuclear * 0.8) if min_nuclear > 0 else 1_000_000
        print(f"[info] chr_like_minlen not specified, using 80% of smallest nuclear chrom: {chr_like_minlen}", file=sys.stderr)
    else:
        chr_like_minlen = args.chr_like_minlen

    # --- Identify chromosome contigs early (before BLAST phases) ---
    # This allows creating a filtered FASTA for faster organelle/rDNA/contaminant detection
    chromosome_contigs: Set[str] = set()
    for contig, ref_id in best_ref.items():
        if ref_id and qry_lengths.get(contig, 0) >= chr_like_minlen:
            chromosome_contigs.add(contig)

    # Create work directory for classification outputs
    work_dir = outprefix.parent / f"{outprefix.name}_classification"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Create filtered FASTA with only non-chromosome contigs for BLAST-based detection
    # This significantly speeds up organelle/rDNA/contaminant BLAST searches
    non_chrom_contigs = set(qry_lengths.keys()) - chromosome_contigs
    non_chrom_fasta = work_dir / "non_chromosome_contigs.fa"
    if non_chrom_contigs:
        print(f"[info] Creating filtered FASTA for BLAST ({len(non_chrom_contigs)} non-chromosome contigs)", file=sys.stderr)
        write_filtered_fasta(qry, non_chrom_fasta, non_chrom_contigs)
        blast_query_fasta = non_chrom_fasta
        blast_query_lengths = {k: v for k, v in qry_lengths.items() if k in non_chrom_contigs}
    else:
        # Edge case: all contigs are chromosomes
        blast_query_fasta = qry
        blast_query_lengths = qry_lengths

    # --- Phase 4: Organelle detection ---
    chrC_contig: Optional[str] = None
    chrM_contig: Optional[str] = None
    organelle_debris: Set[str] = set()

    if not args.skip_organelles:
        print("[info] Phase 4: Organelle detection", file=sys.stderr)
        chrC_ref, chrM_ref = prepare_organelle_references(
            ref_fasta=ref,
            chrC_ref_arg=getattr(args, 'chrC_ref', None),
            chrM_ref_arg=getattr(args, 'chrM_ref', None),
            work_dir=work_dir / "organelles",
            ref_norm_to_orig=ref_norm_to_orig,
        )

        if chrC_ref or chrM_ref:
            chrC_contig, chrM_contig, organelle_debris = detect_organelles(
                query_fasta=blast_query_fasta,
                query_lengths=blast_query_lengths,
                chrC_ref=chrC_ref,
                chrM_ref=chrM_ref,
                work_dir=work_dir / "organelles",
                threads=args.threads,
                min_coverage=args.organelle_min_cov,
                chrC_len_tol=getattr(args, 'chrC_len_tolerance', 0.05),
                chrM_len_tol=getattr(args, 'chrM_len_tolerance', 0.20),
            )
    else:
        print("[info] Phase 4: Skipping organelle detection (--skip-organelles)", file=sys.stderr)

    # --- Phase 5: rDNA detection ---
    rdna_contigs: Set[str] = set()
    already_classified = {chrC_contig, chrM_contig} | organelle_debris
    already_classified.discard(None)

    if not args.skip_rdna:
        print("[info] Phase 5: rDNA detection", file=sys.stderr)
        script_dir = Path(__file__).resolve().parent
        rdna_ref = prepare_rdna_reference(args.rdna_ref, script_dir)

        if rdna_ref:
            rdna_contigs = detect_rdna_contigs(
                query_fasta=blast_query_fasta,
                query_lengths=blast_query_lengths,
                rdna_ref=rdna_ref,
                work_dir=work_dir / "rdna",
                threads=args.threads,
                min_coverage=args.rdna_min_cov,
                exclude_contigs=already_classified,
            )
    else:
        print("[info] Phase 5: Skipping rDNA detection (--skip-rdna)", file=sys.stderr)

    # --- Phase 6: Chromosome debris detection ---
    print("[info] Phase 6: Chromosome debris detection", file=sys.stderr)
    already_classified = already_classified | rdna_contigs

    chromosome_debris: Set[str] = set()
    if chromosome_contigs:
        chromosome_debris = detect_chromosome_debris(
            query_fasta=qry,
            query_lengths=qry_lengths,
            chromosome_contigs=chromosome_contigs,
            work_dir=work_dir / "chr_debris",
            threads=args.threads,
            min_coverage=args.chr_debris_min_cov,
            min_identity=args.chr_debris_min_identity,
            exclude_contigs=already_classified,
        )
        already_classified = already_classified | chromosome_debris

    # --- Phase 7: Contaminant detection ---
    contaminants: Dict[str, Tuple[int, str]] = {}

    if not args.skip_contaminants and args.centrifuger_idx:
        print("[info] Phase 7: Contaminant detection", file=sys.stderr)
        residual_contigs = set(qry_lengths.keys()) - already_classified - chromosome_contigs
        if residual_contigs:
            residual_fasta = work_dir / "residual_for_contaminant_screen.fa"
            write_filtered_fasta(qry, residual_fasta, residual_contigs)
            residual_lengths = {k: v for k, v in qry_lengths.items() if k in residual_contigs}

            contaminants = detect_contaminants(
                query_fasta=residual_fasta,
                query_lengths=residual_lengths,
                centrifuger_idx=args.centrifuger_idx,
                work_dir=work_dir / "contaminants",
                threads=args.threads,
                min_score=args.contaminant_min_score,
                exclude_contigs=set(),
            )
    else:
        print("[info] Phase 7: Skipping contaminant detection", file=sys.stderr)

    # --- Phase 8: Debris/unclassified classification ---
    print("[info] Phase 8: Debris/unclassified classification", file=sys.stderr)

    already_classified = already_classified | set(contaminants.keys()) | chromosome_contigs
    remaining_contigs = set(qry_lengths.keys()) - already_classified

    additional_debris, _unclassified = classify_debris_and_unclassified(
        remaining_contigs=remaining_contigs,
        query_fasta=qry,
        query_lengths=qry_lengths,
        ref_fasta=ref,
        chrs_fasta=None,
        qr_gene_count=ev.qr_gene_count,
        work_dir=work_dir / "debris",
        threads=args.threads,
        min_coverage=args.debris_min_cov,
        min_protein_hits=args.debris_min_protein_hits,
    )
    other_debris = additional_debris

    # --- Phase 9: Gene count statistics ---
    print("[info] Phase 9: Gene count statistics", file=sys.stderr)
    ref_gene_counts = count_genes_per_ref_chrom(ref_gff3, ref_id_map=ref_orig_to_norm)
    compute_mean_genes_per_Mbp(
        qr_gene_count=ev.qr_gene_count,
        query_lengths=qry_lengths,
        chromosome_contigs=chromosome_contigs,
        best_ref=best_ref,
    )

    # --- Phase 10: Orientation determination ---
    print("[info] Phase 10: Orientation determination", file=sys.stderr)
    contig_orientations = determine_contig_orientations(
        macro_block_rows=ev.macro_block_rows,
        best_ref=best_ref,
        chromosome_contigs=chromosome_contigs,
    )

    # --- Phase 11: Classify all contigs ---
    print("[info] Phase 11: Classifying all contigs", file=sys.stderr)
    classifications = classify_all_contigs(
        query_fasta=qry,
        query_lengths=qry_lengths,
        best_ref=best_ref,
        chr_like_minlen=chr_like_minlen,
        ev=ev,
        ref_gene_counts=ref_gene_counts,
        chrC_contig=chrC_contig,
        chrM_contig=chrM_contig,
        organelle_debris=organelle_debris,
        rdna_contigs=rdna_contigs,
        contaminants=contaminants,
        chromosome_debris=chromosome_debris,
        other_debris=other_debris,
        add_subgenome_suffix=args.add_subgenome_suffix,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    for clf in classifications:
        clf.reversed = contig_orientations.get(clf.original_name, False)

    # --- Phase 12: Write FASTA outputs ---
    print("[info] Phase 12: Writing FASTA outputs", file=sys.stderr)
    write_classified_fastas(
        query_fasta=qry,
        classifications=classifications,
        contig_orientations=contig_orientations,
        output_prefix=outprefix,
    )

    # --- Phase 13: Write enhanced summary TSV ---
    print("[info] Phase 13: Writing enhanced summary TSV", file=sys.stderr)
    summary_tsv = Path(str(outprefix) + ".contig_summary.tsv")

    write_contig_summary_tsv(
        output_path=summary_tsv,
        classifications=classifications,
        contig_orientations=contig_orientations,
        all_contig_lengths=qry_lengths,
        contig_total=ev.contig_total,
        contig_refs=ev.contig_refs,
        qlens_from_paf=ev.qlens_from_paf,
        qr_union_bp=ev.qr_union_bp,
        qr_matches=ev.qr_matches,
        qr_alnlen=ev.qr_alnlen,
        qr_gene_count=ev.qr_gene_count,
        best_ref=best_ref,
        best_score=best_score,
        best_bp=best_bp,
        second_ref=second_ref,
        second_score=second_score,
        second_bp=second_bp,
        chr_like_minlen=chr_like_minlen,
        chimera_primary_frac=args.chimera_primary_frac,
        chimera_secondary_frac=args.chimera_secondary_frac,
        assign_min_frac=args.assign_min_frac,
        assign_min_ratio=args.assign_min_ratio,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    print(f"[done] Summary:           {summary_tsv}", file=sys.stderr)
    print(f"[done] Ref lengths:       {ref_lengths_tsv}", file=sys.stderr)

    clf_counts: Dict[str, int] = defaultdict(int)
    for clf in classifications:
        clf_counts[clf.classification] += 1
    print("[done] Classification summary:", file=sys.stderr)
    for cat, count in sorted(clf_counts.items()):
        print(f"       {cat}: {count}", file=sys.stderr)

    if args.plot:
        run_plot(
            summary_tsv,
            ref_lengths_tsv,
            segments_tsv,
            chain_summary_tsv,
            macro_blocks_tsv,
            outprefix,
            chr_like_minlen,
            plot_suffix,
            args.plot_html,
        )


if __name__ == "__main__":
    main()
