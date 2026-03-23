#!/usr/bin/env python3
"""
Data models for final_finalizer.

Contains dataclasses for synteny blocks, chains, and classification results.
These have zero dependencies on other modules, making them safe to import anywhere.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ----------------------------
# Shared reference context for multi-assembly support
# ----------------------------
@dataclass
class ReferenceContext:
    """Shared reference state prepared once and reused across assemblies.

    Holds everything derived from the reference genome: chromosome lengths,
    ID normalization maps, GC baseline, GFF3 annotations, extracted proteins,
    organelle/rDNA reference sequences, and computed thresholds.
    """
    ref: Path
    ref_lengths_norm: Dict[str, int]
    ref_orig_to_norm: Dict[str, str]
    ref_norm_to_orig: Dict[str, str]
    ref_ids_raw: Set[str]
    ref_gc_all: Dict[str, float]
    ref_gc_mean: Optional[float]
    ref_gc_std: Optional[float]
    ref_gff3: Optional[Path]          # filtered GFF3 (or None)
    tx2loc: Dict                       # transcript → (chrom, start, end, strand)
    tx2gene: Dict                      # transcript → gene_id
    proteins_faa: Optional[Path]       # extracted proteins (protein mode only)
    chrC_ref: Optional[Path]
    chrM_ref: Optional[Path]
    rdna_ref: Optional[Path]
    chr_like_minlen: int
    ref_lengths_tsv: Path


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
class TelomereResult:
    """Telomere detection result for a contig.

    Telomeres are repetitive sequences at chromosome ends. Their presence
    indicates assembly completeness - a contig with both telomeres represents
    a complete chromosome regardless of reference coverage.

    Attributes:
        has_5p_telomere: Telomere detected at 5' end
        has_3p_telomere: Telomere detected at 3' end
        telomere_5p_count: Number of telomere repeats at 5' end
        telomere_3p_count: Number of telomere repeats at 3' end
    """
    has_5p_telomere: bool
    has_3p_telomere: bool
    telomere_5p_count: int
    telomere_3p_count: int


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
    # Reference span coverage (extent-based, robust to divergence)
    qr_ref_span_bp: Dict[Tuple[str, str], int] = None  # (contig, ref_id) -> reference span bp
    # Best chain identity per (contig, ref_id) for subgenome inference
    qr_best_chain_ident: Dict[Tuple[str, str], float] = None  # (contig, ref_id) -> identity
    # Reference coordinate ranges per (contig, ref_id) for scaffolding
    qr_ref_ranges: Dict[Tuple[str, str], Tuple[int, int]] = None  # (contig, ref_id) -> (ref_min, ref_max)
    # Collinearity score per (contig, ref_id) — fraction of aligned bp in reference order
    qr_collinearity: Dict[Tuple[str, str], float] = None  # (contig, ref_id) -> score [0.0, 1.0]


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
    # Full-length vs fragment classification
    ref_coverage: Optional[float] = None  # Fraction of reference chromosome spanned (0.0-1.0)
    is_full_length: Optional[bool] = None  # True if full-length, False if fragment
    full_length_confidence: Optional[str] = None  # "high"/"medium"/"low"
    has_5p_telomere: Optional[bool] = None  # Telomere at 5' end
    has_3p_telomere: Optional[bool] = None  # Telomere at 3' end
    # Query subgenome inference (when multiple contigs map to same ref)
    query_subgenome: Optional[str] = None  # "B", "C", etc. or None if primary/single
    query_subgenome_grp: Optional[int] = None  # Numeric cluster ID (1, 2, 3...)
    seq_identity_vs_ref: Optional[float] = None  # Sequence identity (0.0-1.0)
    # Collinearity score for assigned reference (fraction of aligned bp in reference order)
    collinearity_score: Optional[float] = None  # 0.0-1.0, None if no chains
    # Rearrangement hypothesis detection
    rearrangement_candidates: Optional[str] = None  # Comma-separated off-target chromosomes (e.g., "chr2A,chr5B")


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
class RdnaSubFeature:
    """A sub-feature within a 45S rDNA repeat unit.

    Attributes:
        name: Feature name (e.g., "18S", "ITS1", "5.8S", "ITS2", "28S")
        start: Start position (0-based) relative to consensus
        end: End position (0-based, exclusive) relative to consensus
    """
    name: str
    start: int
    end: int


@dataclass
class RdnaSubFeatureLocus:
    """A sub-feature with coordinates mapped to a specific contig locus.

    Attributes:
        name: Feature name (e.g., "18S_rRNA", "ITS1", "5.8S_rRNA", "ITS2", "28S_rRNA")
        start: Start position on contig (0-based)
        end: End position on contig (0-based, exclusive)
    """
    name: str
    start: int
    end: int


@dataclass
class RdnaConsensus:
    """Consensus 45S rDNA repeat unit derived from the query assembly.

    Attributes:
        sequence: The consensus/exemplar DNA sequence
        length: Length of the consensus sequence
        n_copies_extracted: Number of individual 45S copies extracted
        n_copies_clustered: Number of copies in the largest cluster
        method: Method used to derive consensus ("cdhit+mafft", "cdhit_rep", "blast_central")
        sub_features: List of annotated sub-features (18S, ITS1, 5.8S, ITS2, 28S)
        source_contig: Contig from which the exemplar was extracted (if applicable)
    """
    sequence: str
    length: int
    n_copies_extracted: int
    n_copies_clustered: int
    method: str
    sub_features: List[RdnaSubFeature]
    source_contig: Optional[str] = None


@dataclass
class RdnaLocus:
    """An rDNA locus annotation on a contig.

    Represents a single rDNA hit found by BLASTing the consensus 45S
    against all assembly contigs.

    Attributes:
        contig: Contig name
        start: Start position on contig (0-based)
        end: End position on contig (0-based, exclusive)
        strand: "+" or "-"
        identity: Alignment identity to consensus (0.0-1.0)
        consensus_coverage: Fraction of consensus covered (0.0-1.0)
        copy_type: "full", "partial", or "fragment"
        sub_feature_loci: List of sub-features with mapped contig coordinates
        array_id: ID of the 45S rDNA array this locus belongs to, or None
    """
    contig: str
    start: int
    end: int
    strand: str
    identity: float
    consensus_coverage: float
    copy_type: str
    sub_feature_loci: List[RdnaSubFeatureLocus]
    array_id: Optional[str] = None

    @property
    def is_nor_candidate(self) -> bool:
        """Backward-compatible property: True if part of an rDNA array."""
        return self.array_id is not None

    @property
    def sub_features(self) -> List[str]:
        """Return list of sub-feature names for backward compatibility."""
        return [sf.name for sf in self.sub_feature_loci]


@dataclass
class RdnaArray:
    """A tandem 45S rDNA repeat array on a contig.

    Represents a cluster of tandem rDNA loci detected by proximity-based
    clustering. Arrays are first-class objects used for GFF3 annotation
    (as repeat_region parent features) and summary reporting.

    Attributes:
        array_id: Unique identifier (e.g., "array_1", "array_2")
        contig: Contig name
        start: Start position on contig (0-based, first locus start)
        end: End position on contig (0-based exclusive, last locus end)
        strand: Majority strand of constituent loci
        loci: Constituent loci, sorted by start position
        n_total: Total number of loci in array
        n_full: Number of full-length copies (copy_type == "full")
        n_partial: Number of partial copies
        n_fragment: Number of fragment copies
        identity_median: Median identity of constituent loci
        identity_min: Minimum identity
        identity_max: Maximum identity
    """
    array_id: str
    contig: str
    start: int
    end: int
    strand: str
    loci: List[RdnaLocus]
    n_total: int
    n_full: int
    n_partial: int
    n_fragment: int
    identity_median: float
    identity_min: float
    identity_max: float

    @property
    def span(self) -> int:
        return self.end - self.start


@dataclass
class CompleasmResult:
    """BUSCO completeness assessment from compleasm.

    Percentages are stored as parsed from compleasm's output rather than
    derived from counts, to preserve compleasm's own rounding.
    """
    lineage: str
    n_total: int          # N (total BUSCOs in lineage)
    n_single: int         # S count
    n_duplicated: int     # D count
    n_fragmented: int     # F count
    n_interspersed: int   # I count
    n_missing: int        # M count
    pct_single: float
    pct_duplicated: float
    pct_fragmented: float
    pct_interspersed: float
    pct_missing: float
    summary_path: Path

    def summary_line(self) -> str:
        """One-line summary string for logging (e.g. 'S:85.3% D:7.0% F:3.9% M:3.1%')."""
        return (
            f"S:{self.pct_single:.1f}% D:{self.pct_duplicated:.1f}% "
            f"F:{self.pct_fragmented:.1f}% M:{self.pct_missing:.1f}%"
        )

    def tsv_fields(self) -> list[str]:
        """Return [N, S%, D%, F%, I%, M%] as formatted strings for TSV output."""
        def _fmt(val: float) -> str:
            return f"{val:.2f}"
        return [
            str(self.n_total),
            _fmt(self.pct_single),
            _fmt(self.pct_duplicated),
            _fmt(self.pct_fragmented),
            _fmt(self.pct_interspersed),
            _fmt(self.pct_missing),
        ]


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


# ----------------------------
# Multi-assembly comparison
# ----------------------------
@dataclass
class ChromRefSummary:
    """Per-reference-chromosome summary for one assembly."""
    ref_id: str
    ref_length: int
    n_contigs: int
    total_assigned_bp: int
    best_ref_coverage: Optional[float]
    is_full_length: bool
    has_both_telomeres: bool
    mean_identity: Optional[float]


@dataclass
class AssemblyResult:
    """Summary metrics from one assembly's finalization."""
    # Identity
    assembly_name: str
    assembly_path: Path
    outprefix: Path

    # Contiguity
    total_contigs: int
    total_bp: int
    n50: int
    l50: int
    largest_contig: int

    # Classification counts/bp
    classification_counts: Dict[str, int]
    classification_bp: Dict[str, int]

    # Chromosome completeness
    n_chrom_assigned: int
    n_chrom_unassigned: int
    n_full_length: int
    n_with_both_telomeres: int
    n_with_any_telomere: int
    mean_ref_coverage: Optional[float]
    n_chimeric: int

    # Quality (chrom_assigned contigs)
    mean_identity: Optional[float]
    mean_collinearity: Optional[float]
    mean_gc_deviation: Optional[float]

    # Contamination
    n_contaminants: int
    total_contaminant_bp: int
    n_unique_contaminant_species: int

    # rDNA
    n_rdna_contigs: int
    total_rdna_bp: int
    n_rdna_arrays: int

    # Organelles
    chrC_found: bool
    chrM_found: bool

    # Read depth (optional)
    mean_chrom_depth: Optional[float]

    # Per-reference-chromosome detail
    chrom_ref_coverage: Dict[str, ChromRefSummary]

    # Full classifications list (for per-contig comparisons in Rmd)
    classifications: List[ContigClassification]

    # Output file paths (for Rmd template)
    summary_tsv: Path
    segments_tsv: Path
    evidence_tsv: Path
    macro_blocks_tsv: Path
    contaminants_tsv: Optional[Path] = None
    rdna_annotations_tsv: Optional[Path] = None
    rdna_arrays_tsv: Optional[Path] = None
    agp_tsv: Optional[Path] = None

    # Per-subgenome chrs.fasta paths (polyploid references)
    per_subgenome_chrs: Dict[str, Path] = field(default_factory=dict)

    # Compleasm (BUSCO) completeness results
    compleasm_chrs: Optional['CompleasmResult'] = None
    compleasm_non_chrs: Optional['CompleasmResult'] = None
