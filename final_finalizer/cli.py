#!/usr/bin/env python3
"""Command-line interface for final_finalizer.

This module provides the main entry point for the final_finalizer tool,
which classifies genome assembly contigs into biological categories using
synteny evidence (protein-anchored or nucleotide whole-genome alignment),
organelle alignments, and taxonomic screening.

Each contig receives:
- A classification category (chrom_assigned, organelle_complete, debris, etc.)
- A confidence level (high/medium/low) based on evidence strength
- A new name reflecting its assignment (e.g., chr5A, chrC, contig_1)
- Optional read depth metrics (if --reads provided)

Classification confidence is determined by multiple lines of evidence:
- Gene proportion: fraction of reference genes aligned to the contig
- GC content deviation: how many standard deviations from expected baseline
- Alignment coverage and identity
- Protein hit counts

Read depth analysis (optional):
- Auto-detects read format (FASTQ, BAM, CRAM)
- Downsamples to target coverage before alignment (default 20X)
- Aligns with minimap2 using appropriate preset for read type
- Computes depth statistics with mosdepth (mean, median, std, breadth)

See README.md for full documentation and examples.
"""
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
    calculate_gc_content_fasta,
    calculate_gc_stats,
    check_gc_cache,
    read_fasta_lengths,
    read_gc_content_tsv,
    write_filtered_fasta,
    write_gc_content_tsv,
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
    parse_paf_chain_evidence_and_segments,
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
    compute_mean_gene_proportion,
    count_genes_per_ref_chrom,
    determine_contig_orientations,
)
from final_finalizer.output.fasta_output import write_classified_fastas
from final_finalizer.output.tsv_output import (
    build_segment_support_from_rows,
    write_chain_segments_tsv,
    write_chain_summary_tsv,
    write_contig_summary_tsv,
    write_contaminant_summary_tsv,
    write_macro_blocks_tsv,
)
from final_finalizer.output.plotting import (
    run_plot,
    run_depth_plot,
    run_contaminant_plot,
    run_contaminant_table,
)


def main():
    # Handle --dump-config early, before argument parsing
    # This allows generating config template without providing required args
    if '--dump-config' in sys.argv:
        from final_finalizer.utils.config import dump_config_template
        # Create a minimal namespace with defaults
        defaults = argparse.Namespace(
            ref=None, query=None, outprefix=None, ref_gff3=None,
            threads=8, plot=False, plot_html=False, verbose=False, quiet=False,
            log_file=None, config=None, dump_config=True, chr_like_minlen=None,
            add_subgenome_suffix=None, ref_id_pattern=None, reads=None,
            reads_type="lrhq", skip_depth=False, depth_window_size=1000,
            depth_target_coverage=0, keep_depth_bam=False, synteny_mode="protein",
            skip_organelles=False, skip_rdna=False, skip_contaminants=False,
            gffread="gffread", miniprot="miniprot", miniprot_args="",
            chrC_ref=None, chrM_ref=None, rdna_ref=None, centrifuger_idx=None,
            assign_min_frac=0.10, assign_min_ratio=1.25, chimera_primary_frac=0.8,
            chimera_secondary_frac=0.2, low_ref_span_threshold=0.75,
            assign_minlen=10000, assign_minlen_protein=150, assign_minmapq=0,
            assign_min_ident=0.0, assign_tp="PI", chain_q_gap=200000,
            chain_r_gap=400000, chain_diag_slop=150000, assign_chain_min_bp=0,
            assign_chain_score="matches", assign_chain_topk=3, assign_ref_score="all",
            miniprot_min_genes=3, miniprot_min_segments=5, min_span_frac=0.20,
            min_span_bp=50000, organelle_min_cov=0.80, chrC_len_tolerance=0.05,
            chrM_len_tolerance=0.20, rdna_min_cov=0.50, chr_debris_min_cov=0.80,
            chr_debris_min_identity=0.90, contaminant_min_score=1000, contaminant_min_coverage=0.50,
            debris_min_cov=0.50, debris_min_protein_hits=2, preset="asm20", kmer=None, window=None, aln_minlen=10000,
        )
        print(dump_config_template(defaults))
        sys.exit(0)

    p = argparse.ArgumentParser(
        description="Genome assembly finalization tool for contig classification using synteny evidence (protein-anchored or nucleotide whole-genome alignment).",
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
        required=False,  # Only required for --synteny-mode protein
        help="Reference GFF3 with protein-coding gene annotations. "
        "Required for --synteny-mode protein. "
        "Used to extract proteins (gffread) and run protein-anchor synteny blocks (miniprot).",
    )

    # =========================================================================
    # Common options
    # =========================================================================
    common = p.add_argument_group("Common options")
    common.add_argument("-t", "--threads", type=_positive_int, default=8, help="Threads for minimap2/miniprot [8]")
    common.add_argument("--plot", action="store_true", help="Generate overview plots with R/ggplot2")
    common.add_argument("--plot-html", action="store_true", help="Also generate interactive HTML plot (ggiraph)")
    common.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG level) logging")
    common.add_argument("--quiet", action="store_true", help="Suppress INFO messages (only show warnings and errors)")
    common.add_argument("--log-file", type=str, default=None, metavar="PATH", help="Also write logs to this file")
    common.add_argument("--config", type=Path, default=None, metavar="TOML", help="Load configuration from TOML file (CLI args override)")
    common.add_argument("--dump-config", action="store_true", help="Print TOML config template and exit")
    common.add_argument(
        "-C", "--chr-like-minlen", type=int, default=None,
        help="Minimum contig length (bp) to be considered chromosome-like. Default: 80%% of smallest nuclear ref chromosome.",
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
    # Read depth analysis
    # =========================================================================
    depth_grp = p.add_argument_group("Read depth analysis")
    depth_grp.add_argument(
        "--reads", type=_validate_input_path, default=None,
        help="Reads for depth analysis (FASTQ/BAM/CRAM). Auto-detects format and alignment status.",
    )
    depth_grp.add_argument(
        "--reads-type", choices=["lrhq", "r9", "sr"], default="lrhq",
        help="Read type for minimap2 alignment: lrhq (lr:hqae for HiFi/ONT Q20+), r9 (map-ont), sr (sr) [lrhq]",
    )
    depth_grp.add_argument(
        "--skip-depth", action="store_true",
        help="Skip depth analysis even if --reads provided.",
    )
    depth_grp.add_argument(
        "--depth-window-size", type=int, default=1000,
        help="Window size for mosdepth depth calculation [1000]",
    )
    depth_grp.add_argument(
        "--depth-target-coverage", type=float, default=0,
        help="Target coverage (e.g. 20) for downsampling reads before alignment (0 to disable) [0]",
    )
    depth_grp.add_argument(
        "--keep-depth-bam",
        action="store_true",
        help="Keep aligned BAM file after depth analysis (default: delete to save space)",
    )

    # =========================================================================
    # Synteny mode selection
    # =========================================================================
    synteny = p.add_argument_group("Synteny mode")
    synteny.add_argument(
        "--synteny-mode",
        choices=["protein", "nucleotide"],
        default="protein",
        help="Synteny evidence source for classification. "
        "protein: Use miniprot protein-anchored synteny (requires --ref-gff3). "
        "nucleotide: Use minimap2 whole-genome nucleotide alignment with permissive chaining "
        "for chromosome-scale structural composition analysis. Creates megabase-scale blocks "
        "by chaining through repetitive regions. Suitable for both within-species and cross-species comparisons. "
        "[protein]",
    )

    # =========================================================================
    # Pipeline phase toggles
    # =========================================================================
    toggles = p.add_argument_group("Pipeline phase toggles")
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
        help="Min alignment span for nucleotide synteny block building [10000]",
    )
    block_thresh.add_argument(
        "--assign-minlen-protein", type=int, default=150,
        help="Min target span for protein-anchor synteny block building [150]",
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
        help="Min number of synteny segments in protein mode (nucleotide mode uses ≥1) [5]",
    )
    prot_thresh.add_argument(
        "--min-span-frac", type=float, default=0.20,
        help="Min span fraction of contig for chromosome assignment (both modes) [0.20]",
    )
    prot_thresh.add_argument(
        "--min-span-bp", type=int, default=50_000,
        help="Min absolute span in bp for chromosome assignment (both modes) [50000]",
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
        "--contaminant-min-score", type=int, default=1000,
        help="Min centrifuger score for contaminant classification (~1kb matching sequence with k=31) [1000]",
    )
    contam_thresh.add_argument(
        "--contaminant-min-coverage", type=float, default=0.50,
        help="Min query coverage for contaminant classification (low coverage may indicate conserved genes, not contamination) [0.50]",
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
    # Minimap2 tuning (for --synteny-mode nucleotide)
    # =========================================================================
    mm2_tune = p.add_argument_group("Minimap2 tuning (for --synteny-mode nucleotide)")
    mm2_tune.add_argument("-x", "--preset", default="asm20", help="minimap2 preset [asm20]")
    mm2_tune.add_argument("-k", "--kmer", type=int, default=None, help="minimap2 k-mer size")
    mm2_tune.add_argument("-w", "--window", type=int, default=None, help="minimap2 minimizer window size")

    # =========================================================================
    # Full-length detection and subgenome inference
    # =========================================================================
    fl_grp = p.add_argument_group("Full-length detection and subgenome inference")
    fl_grp.add_argument(
        "--full-length-ref-coverage", type=float, default=0.70,
        help="Reference span coverage threshold for full-length classification [0.70]",
    )
    fl_grp.add_argument(
        "--disable-telomere-detection", action="store_true",
        help="Disable telomere detection (enabled by default)",
    )
    fl_grp.add_argument(
        "--telomere-motif", type=str, default="TTTAGGG",
        help="Telomere repeat motif [TTTAGGG for plants, TTAGGG for vertebrates]",
    )
    fl_grp.add_argument(
        "--telomere-window", type=int, default=10000,
        help="Window size at contig ends to search for telomeres [10000]",
    )
    fl_grp.add_argument(
        "--telomere-min-repeats", type=int, default=3,
        help="Minimum consecutive telomere repeats to call telomere [3]",
    )
    fl_grp.add_argument(
        "--subgenome-k", type=float, default=1.0,
        help="Multiplier for std_dev-based query subgenome clustering threshold [1.0]",
    )
    fl_grp.add_argument(
        "--rearrangement-threshold", type=float, default=0.10,
        help="Minimum off-target span fraction to flag as rearrangement candidate [0.10]",
    )

    args = p.parse_args()

    # Load config file if provided (CLI args override config)
    if args.config:
        from final_finalizer.utils.config import load_config, merge_config_with_args
        try:
            config = load_config(args.config)
            merge_config_with_args(config, args)
        except Exception as e:
            sys.exit(f"[error] Failed to load config file: {e}")

    # Validate synteny mode requirements
    if args.synteny_mode == "protein" and not args.ref_gff3:
        sys.exit("[error] --ref-gff3 is required when using --synteny-mode protein")

    # Validate contaminant coverage threshold
    if not 0.0 <= args.contaminant_min_coverage <= 1.0:
        sys.exit(f"[error] --contaminant-min-coverage must be between 0.0 and 1.0, got {args.contaminant_min_coverage}")

    # Setup logging
    from final_finalizer.utils.logging_config import setup_logging, get_logger
    setup_logging(verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    logger = get_logger("cli")

    if args.ref_id_pattern:
        set_ref_id_patterns(compile_ref_id_patterns(args.ref_id_pattern))

    # --- Phase 1: Reading reference and query FASTA ---
    logger.phase("Phase 1: Reading reference and query FASTA")
    ref = args.ref  # Already a Path from _validate_input_path
    qry = args.query  # Already a Path from _validate_input_path
    outprefix = Path(args.outprefix)

    # Create output directory if it doesn't exist
    out_dir = outprefix.parent
    if out_dir and not out_dir.exists():
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created output directory: {out_dir}")
        except PermissionError as e:
            sys.exit(f"[error] Cannot create output directory {out_dir}: {e}")
        except OSError as e:
            sys.exit(f"[error] Failed to create output directory {out_dir}: {e}")

    ref_lengths_norm, ref_orig_to_norm, ref_norm_to_orig = read_fasta_lengths_with_map(ref)
    ref_ids_raw = set(read_fasta_lengths(ref).keys())
    qry_lengths = read_fasta_lengths(qry)
    logger.info(f"Query contigs: {len(qry_lengths)}")

    # GC content cache paths
    gc_content_tsv = Path(str(outprefix) + ".gc_content.tsv")
    gc_content_metadata = Path(str(outprefix) + ".gc_content.metadata.json")

    # Compute or load cached GC content for reference and query
    if check_gc_cache(gc_content_tsv, gc_content_metadata, ref, qry):
        logger.info(f"Reusing cached GC content: {gc_content_tsv}")
        ref_gc_all, qry_gc = read_gc_content_tsv(gc_content_tsv)
    else:
        logger.info("Computing GC content for reference and query sequences")
        ref_gc_all = calculate_gc_content_fasta(ref)
        qry_gc = calculate_gc_content_fasta(qry)
        write_gc_content_tsv(gc_content_tsv, gc_content_metadata, ref_gc_all, qry_gc, ref, qry)
        logger.info(f"Cached GC content: {gc_content_tsv}")

    # Compute reference GC baseline (nuclear chromosomes only, exclude organelles)
    # Filter to only include sequences matching chromosome patterns
    from final_finalizer.utils.reference_utils import is_nuclear_chromosome
    ref_gc_nuclear = {
        ref_norm_to_orig.get(k, k): v
        for k, v in ref_gc_all.items()
        if is_nuclear_chromosome(k)
    }
    if ref_gc_nuclear:
        ref_gc_mean, ref_gc_std = calculate_gc_stats(ref_gc_nuclear)
        logger.info(f"Reference nuclear GC: mean={ref_gc_mean:.3f}, std={ref_gc_std:.3f}")
    else:
        ref_gc_mean, ref_gc_std = None, None
        logger.warning("Could not compute reference GC baseline (no nuclear chromosomes found)")

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
        logger.info(f"Writing reference lengths: {ref_lengths_tsv}")
        write_ref_lengths_tsv(ref_lengths_tsv, ref)
    else:
        logger.info(f"Reference lengths TSV exists, reusing: {ref_lengths_tsv}")

    # --- Phase 2: Synteny analysis (protein or nucleotide) ---
    # Process GFF3 if provided (optional in nucleotide mode, required in protein mode)
    ref_gff3 = None
    tx2loc = {}
    tx2gene = {}
    if args.ref_gff3:
        filtered_ref_gff3 = Path(str(outprefix) + ".ref_gff3.filtered.gff3")
        ref_gff3 = filter_gff3_by_ref(
            args.ref_gff3,
            ref_ids_raw,
            filtered_ref_gff3,
            ref_norm_to_orig=ref_norm_to_orig,
        )
        if args.synteny_mode == "protein":
            logger.info(f"Protein-anchor mode: parsing GFF3 for protein extraction (ref GFF3: {ref_gff3})")
        else:
            logger.info(f"GFF3 provided in nucleotide mode: will use for gene count statistics (ref GFF3: {ref_gff3})")

        tx2loc, tx2gene = parse_gff3_transcript_coords(ref_gff3, ref_id_map=ref_orig_to_norm)
        if not (tx2loc and tx2gene):
            if args.synteny_mode == "protein":
                raise RuntimeError(f"No transcript features found in GFF3 (mRNA/transcript with ID=...). File: {ref_gff3}")
            else:
                logger.warning(f"No transcript features found in GFF3 - gene count statistics will not be available")
                ref_gff3 = None

    if args.synteny_mode == "protein":
        logger.phase("Phase 2: Protein-anchor synteny analysis")
        if not ref_gff3:
            raise RuntimeError("--ref-gff3 is required for protein mode but GFF3 processing failed")

        run_gffread_extract_proteins(args.gffread, ref, ref_gff3, proteins_faa, gffread_err)
        run_miniprot(args.miniprot, qry, proteins_faa, miniprot_paf_gz, args.threads, args.miniprot_args, miniprot_err)

        ev = parse_miniprot_synteny_evidence_and_segments(
            miniprot_paf_gz=miniprot_paf_gz,
            contig_lengths=qry_lengths,
            tx2loc=tx2loc,
            tx2gene=tx2gene,
            assign_minlen=args.assign_minlen_protein,
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
    else:  # nucleotide mode
        logger.phase("Phase 2: Nucleotide synteny analysis (whole-genome alignment)")
        logger.info("Running minimap2 with permissive chaining for chromosome-scale structural analysis")
        logger.info("Note: Nucleotide mode is optimized for within-species or closely related genomes. "
                    "For highly divergent species, protein mode may provide more reliable assignments.")
        run_minimap2_synteny(
            ref, qry, paf_gz, args.threads, args.preset, args.kmer, args.window, err_log,
            permissive=True  # Always use permissive chaining in nucleotide mode
        )

        ev = parse_paf_chain_evidence_and_segments(
            paf_gz_path=paf_gz,
            contig_lengths=qry_lengths,
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
            ref_id_map=ref_orig_to_norm,
        )

    # --- Synteny gates: filter weak assignments ---
    seg_count, span_bp = build_segment_support_from_rows(ev.chain_segments_rows)

    # Set min_segments based on synteny mode:
    # - Protein mode: use configurable --miniprot-min-segments (default 5)
    # - Nucleotide mode: always 1 (perfect alignments often produce single continuous matches)
    if args.synteny_mode == "protein":
        min_segments = args.miniprot_min_segments
        min_genes = args.miniprot_min_genes
    else:  # nucleotide mode
        min_segments = 1
        min_genes = 0  # Gene count not required in nucleotide mode

    def _passes_gate(
        q: str,
        ref_id: str,
        clen: int,
        min_seg: int,
        min_span_bp: int,
        min_span_frac: float,
        min_gen: int,
    ) -> tuple[bool, int, int, float, int]:
        key = (q, ref_id)
        nseg = int(seg_count.get(key, 0) or 0)
        spbp = int(span_bp.get(key, 0) or 0)
        spfrac = (spbp / clen) if clen > 0 else 0.0
        ngen = int(ev.qr_gene_count.get(key, 0) or 0)

        # Gate: contig must meet all thresholds to be assigned
        ok = (nseg >= min_seg) and (spbp >= min_span_bp) and (spfrac >= min_span_frac) and (ngen >= min_gen)
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
    n_no_candidates = 0  # Contigs with no alignment candidates
    qr_ref_score = ev.qr_weight_all if args.assign_ref_score == "all" else ev.qr_score_topk

    candidates_by_contig = defaultdict(list)
    for (q, ref_id), sc in qr_ref_score.items():
        if sc > 0:
            candidates_by_contig[q].append((float(sc), ref_id))

    logger.info(f"Assignment score per ref: {args.assign_ref_score}")

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
            n_no_candidates += 1
            continue

        cands.sort(key=lambda x: (x[0], int(ev.qr_union_bp.get((q, x[1]), 0) or 0)), reverse=True)
        orig_best = cands[0][1]

        passing = []
        for sc, ref_id in cands:
            ok, _nseg, _spbp, _spfrac, _ngen = _passes_gate(
                q,
                ref_id,
                clen,
                min_segments,
                args.min_span_bp,
                args.min_span_frac,
                min_genes,
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

    # Log gate filtering statistics
    n_total = len(qry_lengths)
    n_with_candidates = n_total - n_no_candidates
    n_passing = n_with_candidates - n_demoted
    logger.info(f"Gate filtering: {n_total} contigs, {n_with_candidates} with candidates, "
                f"{n_passing} passing gates, {n_demoted} demoted, {n_switched} switched")
    logger.info(f"Gate thresholds: min_segments={min_segments}, min_span_bp={args.min_span_bp}, "
                f"min_span_frac={args.min_span_frac}, min_genes={min_genes}")

    if args.synteny_mode == "protein":
        plot_suffix = "protein-anchor synteny"
        logger.done(f"Proteins:         {proteins_faa}")
        logger.done(f"miniprot PAF.gz:  {miniprot_paf_gz}")
        logger.done(f"miniprot err:     {miniprot_err}")
    else:  # nucleotide mode
        plot_suffix = "nucleotide synteny"
        logger.done(f"minimap2 PAF.gz:  {paf_gz}")
        logger.done(f"minimap2 err:     {err_log}")

    # Write segments + evidence summary TSVs
    write_chain_segments_tsv(segments_tsv, ev.chain_segments_rows, ref_norm_to_orig=ref_norm_to_orig)
    write_chain_summary_tsv(chain_summary_tsv, ev.chain_summary_rows, ref_norm_to_orig=ref_norm_to_orig)
    logger.done(f"Segments TSV:      {segments_tsv}")
    logger.done(f"Evidence TSV:      {chain_summary_tsv}")

    # Write macro blocks TSV
    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=ref_norm_to_orig)
    logger.done(f"Macro blocks TSV:  {macro_blocks_tsv}")

    # --- Compute chr_like_minlen if not provided ---
    ref_lengths = ref_lengths_norm
    if args.chr_like_minlen is None:
        min_nuclear = get_min_nuclear_chrom_length(ref_lengths)
        chr_like_minlen = int(min_nuclear * 0.25) if min_nuclear > 0 else 100_000
        logger.info(f"chr_like_minlen not specified, using 25% of smallest nuclear chrom: {chr_like_minlen}")
    else:
        chr_like_minlen = args.chr_like_minlen

    # --- Identify chromosome contigs early (before BLAST phases) ---
    # This allows creating a filtered FASTA for faster organelle/rDNA/contaminant detection
    chromosome_contigs: Set[str] = set()
    for contig, ref_id in best_ref.items():
        if ref_id and qry_lengths.get(contig, 0) >= chr_like_minlen:
            chromosome_contigs.add(contig)

    # Compute assembly chromosome GC baseline (for non-chrom classification confidence)
    # This is more appropriate than reference GC for divergent genomes
    asm_gc_nuclear = {k: v for k, v in qry_gc.items() if k in chromosome_contigs}
    if asm_gc_nuclear:
        asm_gc_mean, asm_gc_std = calculate_gc_stats(asm_gc_nuclear)
        logger.info(f"Assembly chromosome GC: mean={asm_gc_mean:.3f}, std={asm_gc_std:.3f} (n={len(asm_gc_nuclear)})")
    else:
        # Fall back to reference GC if no chromosome contigs identified
        asm_gc_mean, asm_gc_std = ref_gc_mean, ref_gc_std
        logger.warning("No chromosome contigs for assembly GC baseline, using reference GC")

    # Create work directory for classification outputs
    work_dir = outprefix.parent / f"{outprefix.name}_classification"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Create filtered FASTA with only non-chromosome contigs for BLAST-based detection
    # This significantly speeds up organelle/rDNA/contaminant BLAST searches
    non_chrom_contigs = set(qry_lengths.keys()) - chromosome_contigs
    non_chrom_fasta = work_dir / "non_chromosome_contigs.fa"
    if non_chrom_contigs:
        logger.info(f"Creating filtered FASTA for BLAST ({len(non_chrom_contigs)} non-chromosome contigs)")
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
    organelle_hits: Dict = {}

    if not args.skip_organelles:
        logger.phase("Phase 4: Organelle detection")
        chrC_ref, chrM_ref = prepare_organelle_references(
            ref_fasta=ref,
            chrC_ref_arg=getattr(args, 'chrC_ref', None),
            chrM_ref_arg=getattr(args, 'chrM_ref', None),
            work_dir=work_dir / "organelles",
            ref_norm_to_orig=ref_norm_to_orig,
        )

        if chrC_ref or chrM_ref:
            chrC_contig, chrM_contig, organelle_debris, organelle_hits = detect_organelles(
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
        logger.info("Phase 4: Skipping organelle detection (--skip-organelles)")

    # --- Phase 5: rDNA detection ---
    rdna_contigs: Set[str] = set()
    rdna_hits: Dict = {}
    already_classified = {chrC_contig, chrM_contig} | organelle_debris
    already_classified.discard(None)

    if not args.skip_rdna:
        logger.phase("Phase 5: rDNA detection")
        script_dir = Path(__file__).resolve().parent
        rdna_ref = prepare_rdna_reference(args.rdna_ref, script_dir)

        if rdna_ref:
            rdna_contigs, rdna_hits = detect_rdna_contigs(
                query_fasta=blast_query_fasta,
                query_lengths=blast_query_lengths,
                rdna_ref=rdna_ref,
                work_dir=work_dir / "rdna",
                threads=args.threads,
                min_coverage=args.rdna_min_cov,
                exclude_contigs=already_classified,
            )
    else:
        logger.info("Phase 5: Skipping rDNA detection (--skip-rdna)")

    # --- Phase 6: Chromosome debris detection ---
    logger.phase("Phase 6: Chromosome debris detection")
    already_classified = already_classified | rdna_contigs

    chromosome_debris: Set[str] = set()
    chrom_debris_hits: Dict = {}
    if chromosome_contigs:
        chromosome_debris, chrom_debris_hits = detect_chromosome_debris(
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
        logger.phase("Phase 7: Contaminant detection")
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
        logger.info("Phase 7: Skipping contaminant detection")

    # --- Phase 8: Debris/unclassified classification ---
    logger.phase("Phase 8: Debris/unclassified classification")

    already_classified = already_classified | set(contaminants.keys()) | chromosome_contigs
    remaining_contigs = set(qry_lengths.keys()) - already_classified

    additional_debris, _unclassified, other_debris_hits = classify_debris_and_unclassified(
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
    logger.phase("Phase 9: Gene count statistics")
    if ref_gff3:
        ref_gene_counts = count_genes_per_ref_chrom(ref_gff3, ref_id_map=ref_orig_to_norm)
        compute_mean_gene_proportion(
            qr_gene_count=ev.qr_gene_count,
            chromosome_contigs=chromosome_contigs,
            best_ref=best_ref,
            ref_gene_counts=ref_gene_counts,
        )
    else:
        logger.info("No GFF3 available - skipping gene count statistics")
        ref_gene_counts = {}

    # --- Phase 10: Orientation determination ---
    logger.phase("Phase 10: Orientation determination")
    contig_orientations = determine_contig_orientations(
        macro_block_rows=ev.macro_block_rows,
        best_ref=best_ref,
        chromosome_contigs=chromosome_contigs,
        query_lengths=qry_lengths,
    )

    # --- Phase 10.5: Telomere detection (optional) ---
    telomere_results = None
    if not args.disable_telomere_detection:
        logger.phase("Phase 10.5: Telomere detection")
        from final_finalizer.detection.telomere import detect_telomeres
        telomere_results = detect_telomeres(
            query_fasta=qry,
            contig_names=chromosome_contigs,  # Only analyze chromosome-assigned contigs
            motif=args.telomere_motif,
            window_size=args.telomere_window,
            min_repeats=args.telomere_min_repeats,
        )
    else:
        logger.info("Phase 10.5: Skipping telomere detection (--disable-telomere-detection)")

    # --- Phase 11: Classify all contigs ---
    logger.phase("Phase 11: Classifying all contigs")

    # Filter contaminants by coverage threshold (low coverage hits may be false positives)
    contaminants_filtered = {
        contig: hit for contig, hit in contaminants.items()
        if hit.coverage >= args.contaminant_min_coverage
    }
    if len(contaminants_filtered) < len(contaminants):
        n_filtered = len(contaminants) - len(contaminants_filtered)
        logger.info(
            f"Filtered {n_filtered} low-coverage contaminant hits "
            f"(coverage < {args.contaminant_min_coverage:.0%})"
        )

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
        contaminants=contaminants_filtered,
        chromosome_debris=chromosome_debris,
        other_debris=other_debris,
        add_subgenome_suffix=args.add_subgenome_suffix,
        ref_norm_to_orig=ref_norm_to_orig,
        query_gc=qry_gc,
        ref_gc_mean=ref_gc_mean,
        ref_gc_std=ref_gc_std,
        asm_gc_mean=asm_gc_mean,
        asm_gc_std=asm_gc_std,
        organelle_hits=organelle_hits,
        rdna_hits=rdna_hits,
        chrom_debris_hits=chrom_debris_hits,
        other_debris_hits=other_debris_hits,
        # Full-length and subgenome inference parameters
        ref_lengths=ref_lengths,
        telomere_results=telomere_results,
        full_length_threshold=args.full_length_ref_coverage,
        subgenome_k=args.subgenome_k,
        rearrangement_threshold=args.rearrangement_threshold,
    )

    for clf in classifications:
        clf.reversed = contig_orientations.get(clf.original_name, False)

    # --- Phase 11.5: Read depth analysis (optional) ---
    depth_stats: Dict[str, 'DepthStats'] = {}
    if args.reads and not args.skip_depth:
        logger.phase("Phase 11.5: Read depth analysis")
        from final_finalizer.analysis.read_depth import calculate_depth_metrics
        from final_finalizer.models import DepthStats

        depth_work_dir = work_dir / "read_depth"
        depth_stats = calculate_depth_metrics(
            reads=args.reads,
            assembly=qry,
            contig_lengths=qry_lengths,
            work_dir=depth_work_dir,
            threads=args.threads,
            reads_type=args.reads_type,
            window_size=args.depth_window_size,
            target_coverage=args.depth_target_coverage if args.depth_target_coverage > 0 else None,
            keep_bam=args.keep_depth_bam,
        )

        # Attach depth stats to classifications
        for clf in classifications:
            ds = depth_stats.get(clf.original_name)
            if ds:
                clf.depth_mean = ds.mean_depth
                clf.depth_median = ds.median_depth
                clf.depth_std = ds.std_depth
                clf.depth_breadth_1x = ds.breadth_1x
                clf.depth_breadth_10x = ds.breadth_10x

        if depth_stats:
            logger.done(f"Depth analysis complete for {len(depth_stats)} contigs")
    elif args.reads and args.skip_depth:
        logger.info("Phase 11.5: Skipping depth analysis (--skip-depth)")

    # --- Phase 12: Write FASTA outputs ---
    logger.phase("Phase 12: Writing FASTA outputs")
    write_classified_fastas(
        query_fasta=qry,
        classifications=classifications,
        contig_orientations=contig_orientations,
        output_prefix=outprefix,
    )

    # --- Phase 13: Write enhanced summary TSV ---
    logger.phase("Phase 13: Writing enhanced summary TSV")
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

    logger.done(f"Summary:           {summary_tsv}")
    logger.done(f"Ref lengths:       {ref_lengths_tsv}")

    # Write contaminant summary TSV with taxonomic lineage (for alluvial plot)
    # Use filtered contaminants (coverage >= threshold) for consistency with classification
    contaminants_tsv = Path(str(outprefix) + ".contaminants.tsv")
    if contaminants_filtered:
        write_contaminant_summary_tsv(
            output_path=contaminants_tsv,
            contaminants=contaminants_filtered,
            query_lengths=qry_lengths,
            depth_stats=depth_stats if depth_stats else None,
        )
        logger.done(f"Contaminants:      {contaminants_tsv}")

    clf_counts: Dict[str, int] = defaultdict(int)
    for clf in classifications:
        clf_counts[clf.classification] += 1
    logger.done("Classification summary:")
    for cat, count in sorted(clf_counts.items()):
        logger.info(f"       {cat}: {count}")

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
        # Generate depth overview plot if depth data was computed
        if depth_stats:
            run_depth_plot(
                summary_tsv,
                outprefix,
                plot_suffix,
                args.plot_html,
            )
        # Generate contaminant plots if contaminants were detected
        # Note: Coverage filtering already applied when writing contaminants_tsv
        if contaminants_filtered and contaminants_tsv.exists():
            run_contaminant_plot(
                contaminants_tsv,
                outprefix,
                plot_suffix,
                args.plot_html,
            )
            run_contaminant_table(
                contaminants_tsv,
                outprefix,
                plot_suffix,
            )


if __name__ == "__main__":
    main()
