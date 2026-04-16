#!/usr/bin/env python3
"""
dnadis - Genome assembly comprehension and curation tool.

This package classifies contigs from de novo genome assemblies using
protein-anchored synteny evidence, organelle detection, and taxonomic screening.
"""
from __future__ import annotations

# Package version.  Keep in sync with pyproject.toml.
__version__ = "0.1.0"

# flake8: noqa: F401
# Import models from package module (zero dependencies)
from dnadis.models import (
    Block,
    Chain,
    ChainEvidenceResult,
    ContigClassification,
    BlastHitSummary,
)

# Import I/O utilities from package module
from dnadis.utils.io_utils import (
    open_maybe_gzip,
    have_exe,
    run_to_gzip,
    merge_intervals,
)

# Import sequence utilities from package module
from dnadis.utils.sequence_utils import (
    reverse_complement,
    read_fasta_lengths,
    read_fasta_sequences,
    write_fasta,
    write_filtered_fasta,
    is_hifiasm_circular,
)

# Import reference utilities from package module
from dnadis.utils.reference_utils import (
    normalize_ref_id,
    build_ref_id_maps,
    normalize_ref_lengths,
    read_fasta_lengths_with_map,
    normalize_organelle_id,
    get_ref_id_patterns,
    set_ref_id_patterns,
    compile_ref_id_patterns,
    get_min_nuclear_chrom_length,
    split_chrom_subgenome,
    filter_gff3_by_ref,
    write_ref_lengths_tsv,
    parse_gff3_transcript_coords,
    paf_tag_value,
    is_primary_only,
)

# Import alignment external tools from package module
from dnadis.alignment.external_tools import (
    get_minimap2_exe,
    run_minimap2,
    run_minimap2_synteny,
    run_gffread_extract_proteins,
    run_miniprot,
)

# Import BLAST infrastructure from package module
from dnadis.detection.blast import (
    run_makeblastdb,
    run_blastn_megablast,
    parse_blast_coverage,
)

# Import detection modules from package
from dnadis.detection.organelle import (
    extract_organelle_from_ref,
    prepare_organelle_references,
    detect_organelles,
)
from dnadis.detection.rdna import (
    prepare_rdna_reference,
    detect_rdna_contigs,
)
from dnadis.detection.debris import detect_chromosome_debris
from dnadis.detection.contaminant import detect_contaminants

# Import chain parsing from package module
from dnadis.alignment.chain_parsing import (
    parse_paf_chain_evidence_and_segments,
    parse_miniprot_synteny_evidence_and_segments,
)

# Import classification from package module
from dnadis.classification.classifier import (
    compute_orientation_votes,
    determine_contig_orientations,
    count_genes_per_ref_chrom,
    compute_mean_gene_proportion,
    classify_debris_and_unclassified,
    generate_contig_names,
    classify_all_contigs,
)

# Import output modules from package
from dnadis.output.fasta_output import write_classified_fastas
from dnadis.output.tsv_output import (
    write_contig_summary_tsv,
    write_macro_blocks_tsv,
    write_chain_segments_tsv,
    write_chain_summary_tsv,
    compute_summary,
    build_segment_support_from_rows,
)
from dnadis.output.plotting import (
    have_rscript,
    run_assembly_report,
)

# Import stats functions from alignment module
from dnadis.alignment.stats import (
    parse_paf_primary,
    write_stats_tsv,
)

# Import main entry point from CLI module
from dnadis.cli import main

__all__ = [
    # Utilities
    "open_maybe_gzip",
    "have_exe",
    "run_to_gzip",
    "merge_intervals",
    # Reference utilities
    "normalize_ref_id",
    "build_ref_id_maps",
    "normalize_ref_lengths",
    "read_fasta_lengths_with_map",
    "normalize_organelle_id",
    "get_ref_id_patterns",
    "set_ref_id_patterns",
    "compile_ref_id_patterns",
    "get_min_nuclear_chrom_length",
    "split_chrom_subgenome",
    "filter_gff3_by_ref",
    "write_ref_lengths_tsv",
    "parse_gff3_transcript_coords",
    # Sequence utilities
    "reverse_complement",
    "read_fasta_lengths",
    "read_fasta_sequences",
    "write_fasta",
    "write_filtered_fasta",
    "is_hifiasm_circular",
    # External tools
    "get_minimap2_exe",
    "run_minimap2",
    "run_minimap2_synteny",
    "run_gffread_extract_proteins",
    "run_miniprot",
    # PAF helpers
    "paf_tag_value",
    "is_primary_only",
    # BLAST
    "run_makeblastdb",
    "run_blastn_megablast",
    "parse_blast_coverage",
    # Detection
    "extract_organelle_from_ref",
    "prepare_organelle_references",
    "detect_organelles",
    "prepare_rdna_reference",
    "detect_rdna_contigs",
    "detect_chromosome_debris",
    "detect_contaminants",
    # Classification
    "compute_orientation_votes",
    "determine_contig_orientations",
    "count_genes_per_ref_chrom",
    "compute_mean_gene_proportion",
    "classify_debris_and_unclassified",
    "generate_contig_names",
    "classify_all_contigs",
    # Output
    "write_classified_fastas",
    "write_contig_summary_tsv",
    "write_macro_blocks_tsv",
    "write_chain_segments_tsv",
    "write_chain_summary_tsv",
    "compute_summary",
    "build_segment_support_from_rows",
    "have_rscript",
    "run_assembly_report",
    # Chain parsing
    "parse_paf_primary",
    "write_stats_tsv",
    "parse_paf_chain_evidence_and_segments",
    "parse_miniprot_synteny_evidence_and_segments",
    # Dataclasses
    "Block",
    "Chain",
    "ChainEvidenceResult",
    "ContigClassification",
    "BlastHitSummary",
    # Main
    "main",
]
