#!/usr/bin/env python3
"""
Data models for final_finalizer.

Contains dataclasses for synteny blocks, chains, and classification results.
These have zero dependencies on other modules, making them safe to import anywhere.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple


# ----------------------------
# Synteny-block chaining (shared with both modes)
# ----------------------------
@dataclass
class Block:
    """Represents a single alignment block from PAF or miniprot output."""
    qs: int
    qe: int
    rs: int
    re: int
    matches: int
    aln_len: int
    mapq: int
    strand: str
    gene_id: Optional[str] = None


@dataclass
class Chain:
    """Represents a chain of alignment blocks."""
    last_qe: int
    last_r: int
    diag0: int
    q_intervals: list[tuple[int, int]]
    matches_sum: int
    alnlen_sum: int
    gene_ids: Set[str]


@dataclass
class ChainEvidenceResult:
    """Result container for chain evidence and segment parsing functions."""
    qlens_from_paf: Dict[str, int]
    qr_union_bp: Dict[Tuple[str, str], int]
    qr_matches: Dict[Tuple[str, str], int]
    qr_alnlen: Dict[Tuple[str, str], int]
    qr_gene_count: Dict[Tuple[str, str], int]
    contig_total: Dict[str, int]
    contig_refs: Dict[str, Set[str]]
    qr_score_topk: Dict[Tuple[str, str], float]
    qr_weight_all: Dict[Tuple[str, str], float]
    qr_nchains_kept: Dict[Tuple[str, str], int]
    best_ref: Dict[str, str]
    best_score: Dict[str, float]
    best_bp: Dict[str, int]
    second_ref: Dict[str, str]
    second_score: Dict[str, float]
    second_bp: Dict[str, int]
    chain_segments_rows: List[Tuple]
    macro_block_rows: List[Tuple]
    chain_summary_rows: List[Tuple]


# ----------------------------
# Contig classification data structures
# ----------------------------
@dataclass
class ContigClassification:
    """Classification result for a single contig.

    Contains the classification category, confidence level, and supporting evidence
    for a query assembly contig.

    Classification categories:
    - chrom_assigned: Chromosome with reference assignment via synteny
    - chrom_unassigned: Chromosome-length without reference assignment
    - organelle_complete: Complete organelle genome (chrC or chrM)
    - organelle_debris: Partial organelle sequence
    - rDNA: Ribosomal DNA repeat unit
    - contaminant: Sequence from contaminating organism
    - chrom_debris: High-coverage duplicate of assembled chromosome
    - debris: Assembly fragment with reference homology
    - unclassified: No classification evidence

    Confidence levels (high/medium/low) are determined by:
    - Gene proportion: Fraction of reference genes aligned (for chromosomes)
    - GC deviation: Standard deviations from reference nuclear GC mean
    - Coverage: Alignment coverage fraction
    - Identity: Alignment identity
    - Protein hits: Number of miniprot alignments
    """
    original_name: str
    new_name: str
    classification: str  # chrom_assigned, chrom_unassigned, chrom_debris, organelle_complete, organelle_debris, rDNA, contaminant, debris, unclassified
    reversed: bool
    contaminant_taxid: Optional[int]  # NCBI taxonomy ID for contaminants
    contaminant_sci: Optional[str]  # Scientific name for contaminants
    assigned_ref_id: Optional[str]
    ref_gene_proportion: Optional[float]
    contig_len: int
    # Evidence strength fields (for classification confidence)
    gc_content: Optional[float] = None  # GC content (0.0-1.0)
    gc_deviation: Optional[float] = None  # Deviation from reference mean in std devs
    synteny_score: Optional[float] = None  # Synteny evidence strength (0.0-1.0)
    contam_score: Optional[float] = None  # Contaminant evidence strength (0.0-1.0)
    contam_coverage: Optional[float] = None  # Contaminant alignment coverage (0.0-1.0)
    classification_confidence: Optional[str] = None  # high/medium/low
    # Read depth fields (optional, from --reads analysis)
    depth_mean: Optional[float] = None
    depth_median: Optional[float] = None
    depth_std: Optional[float] = None
    depth_breadth_1x: Optional[float] = None
    depth_breadth_10x: Optional[float] = None


@dataclass
class ContaminantHit:
    """Contaminant detection result for a contig."""
    taxid: int
    sci_name: str
    coverage: float  # Fraction of contig covered (0.0-1.0)
    score: int  # Centrifuger score


@dataclass
class ContaminantHitExtended:
    """Contaminant hit with full taxonomic lineage.

    Extends ContaminantHit with NCBI taxonomy hierarchy from TaxonKit.
    Used for phylogenetic breakdown visualizations (alluvial plots).
    """
    taxid: int
    sci_name: str
    coverage: float  # Fraction of contig covered (0.0-1.0)
    score: int  # Centrifuger score
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    tax_class: Optional[str] = None  # 'class' is a reserved keyword
    order: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None


@dataclass
class OrganelleHit:
    """Organelle detection result for a contig."""
    organelle_type: str  # "chrC" or "chrM"
    coverage: float  # Fraction of contig covered by alignments (0.0-1.0)
    identity: float  # Alignment identity (0.0-1.0)
    length_ratio: float  # Contig length / reference length
    is_complete: bool  # True if classified as complete organelle


@dataclass
class RdnaHit:
    """rDNA detection result for a contig."""
    coverage: float  # Fraction of contig covered by rDNA alignments (0.0-1.0)
    identity: float  # Alignment identity (0.0-1.0)


@dataclass
class DebrisHit:
    """Debris detection result for a contig."""
    coverage: float  # Fraction of contig covered by alignments (0.0-1.0)
    identity: float  # Alignment identity (0.0-1.0)
    protein_hits: int  # Number of protein hits
    source: str  # "chromosome" or "reference" indicating detection source


@dataclass
class BlastHitSummary:
    """Aggregated BLAST hits for a contig."""
    contig_name: str
    total_coverage: float
    best_hit_subject: str
    best_hit_taxid: Optional[int]
    best_hit_evalue: float
    total_aligned_bp: int


@dataclass
class DepthStats:
    """Read depth statistics for a contig.

    Computed from mosdepth output across window bins.

    Attributes:
        mean_depth: Mean read depth across the contig
        median_depth: Median read depth across the contig
        std_depth: Standard deviation of read depth
        min_depth: Minimum depth in any window
        max_depth: Maximum depth in any window
        breadth_1x: Fraction of bases with >= 1x coverage (0.0-1.0)
        breadth_10x: Fraction of bases with >= 10x coverage (0.0-1.0)
    """
    mean_depth: float
    median_depth: float
    std_depth: float
    min_depth: float
    max_depth: float
    breadth_1x: float
    breadth_10x: float
