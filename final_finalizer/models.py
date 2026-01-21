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
    """Classification result for a single contig."""
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
    contam_coverage: Optional[float] = None  # Contaminant alignment coverage
    classification_confidence: Optional[str] = None  # high/medium/low


@dataclass
class ContaminantHit:
    """Contaminant detection result for a contig."""
    taxid: int
    sci_name: str
    coverage: float  # Fraction of contig covered (0.0-1.0)
    score: int  # Centrifuger score


@dataclass
class BlastHitSummary:
    """Aggregated BLAST hits for a contig."""
    contig_name: str
    total_coverage: float
    best_hit_subject: str
    best_hit_taxid: Optional[int]
    best_hit_evalue: float
    total_aligned_bp: int
