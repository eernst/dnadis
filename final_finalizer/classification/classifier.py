#!/usr/bin/env python3
"""
Classification functions for final_finalizer.

Contains functions for:
- Contig orientation determination
- Gene counting and statistics
- Reference-based debris classification
- Contig naming
- Classification pipeline
"""
from __future__ import annotations

import math
import statistics

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.alignment.external_tools import get_minimap2_exe, run_minimap2
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.models import (
    ChainEvidenceResult,
    ContaminantHit,
    ContigClassification,
    DebrisHit,
    OrganelleHit,
    RdnaHit,
    TelomereResult,
)
from final_finalizer.utils.io_utils import merge_intervals, open_maybe_gzip
from final_finalizer.utils.reference_utils import normalize_ref_id, split_chrom_subgenome
from final_finalizer.utils.sequence_utils import read_fasta_sequences, write_fasta

logger = get_logger("classifier")


# ----------------------------
# Orientation determination
# ----------------------------
def compute_orientation_votes(
    macro_block_rows: List[Tuple],
    contig: str,
    assigned_ref_id: str,
) -> Tuple[int, int]:
    """Compute forward/reverse orientation votes from macro blocks.

    Each macro block contributes its aligned bp (union_bp) as a weighted vote.
    Only blocks matching the assigned reference are counted (off-target blocks
    to other chromosomes are ignored).

    Returns (fwd_bp, rev_bp) — total aligned bp on each strand.
    """
    fwd_bp = 0
    rev_bp = 0

    for row in macro_block_rows:
        # Row format: (contig, contig_len, ref_id, chrom_id, subgenome, strand,
        #              chain_id, qstart, qend, qspan, union_bp, ...)
        if len(row) < 11:
            continue
        row_contig = row[0]
        row_ref_id = row[2]
        row_strand = row[5]
        union_bp = row[10]

        if row_contig != contig:
            continue
        if row_ref_id != assigned_ref_id:
            continue

        if row_strand == "+":
            fwd_bp += union_bp
        elif row_strand == "-":
            rev_bp += union_bp

    return fwd_bp, rev_bp


def determine_contig_orientations(
    macro_block_rows: List[Tuple],
    best_ref: Dict[str, str],
    chromosome_contigs: Set[str],
    query_lengths: Optional[Dict[str, int]] = None,
) -> Dict[str, bool]:
    """Determine which chromosome contigs need to be reverse-complemented.

    Only chromosome-assigned contigs are subject to reorientation based on synteny
    block strand votes. Non-chromosome contigs (debris, contaminants, unclassified)
    are left in their original orientation since they lack reliable synteny evidence.

    Returns:
        Dict[str, bool]: Mapping of contig name to orientation flag.
            True = contig should be reverse-complemented to match reference orientation.
            False = contig is already in correct orientation.
            Only chromosome contigs are included in the returned dict.
    """
    orientations: Dict[str, bool] = {}

    for contig in chromosome_contigs:
        assigned_ref = best_ref.get(contig, "")
        if not assigned_ref:
            orientations[contig] = False
            continue

        fwd_bp, rev_bp = compute_orientation_votes(macro_block_rows, contig, assigned_ref)

        # Reverse if more aligned bp on reverse strand than forward
        should_reverse = rev_bp > fwd_bp
        orientations[contig] = should_reverse

        if should_reverse:
            contig_len = query_lengths.get(contig, 0) if query_lengths else 0
            logger.info(f"Contig {contig} ({contig_len:,} bp) will be reverse-complemented (fwd={fwd_bp:,} bp, rev={rev_bp:,} bp)")

    return orientations


# ----------------------------
# Gene count statistics
# ----------------------------
def count_genes_per_ref_chrom(
    gff3_path: Path,
    ref_id_map: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """Count genes per reference chromosome from GFF3.

    Counts unique gene IDs per chromosome.

    Returns dict: ref_id -> gene_count
    """
    chrom_genes: Dict[str, Set[str]] = defaultdict(set)

    with open_maybe_gzip(gff3_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _src, ftype, _start, _end, _score, _strand, _phase, attrs = parts

            if ftype.lower() != "gene":
                continue

            # Parse gene ID from attributes
            gene_id = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("ID="):
                    gene_id = kv[3:]
                    break

            if gene_id:
                ref_id = ref_id_map.get(seqid, normalize_ref_id(seqid)) if ref_id_map else seqid
                chrom_genes[ref_id].add(gene_id)

    return {chrom: len(genes) for chrom, genes in chrom_genes.items()}


def compute_mean_gene_proportion(
    qr_gene_count: Dict[Tuple[str, str], int],
    chromosome_contigs: Set[str],
    best_ref: Dict[str, str],
    ref_gene_counts: Dict[str, int],
) -> float:
    """Compute mean gene proportion across chromosome-assigned contigs.

    For each chromosome-assigned contig, computes the ratio:
        (genes mapped to contig) / (total genes in assigned reference chromosome)

    Returns the mean of these ratios across all chromosome-assigned contigs.
    """
    ratios = []

    for contig in chromosome_contigs:
        ref_id = best_ref.get(contig, "")
        if not ref_id:
            continue
        mapped_genes = qr_gene_count.get((contig, ref_id), 0)
        ref_genes = ref_gene_counts.get(ref_id, 0)

        if ref_genes > 0:
            ratio = mapped_genes / ref_genes
            ratios.append(ratio)

    if not ratios:
        return 0.0

    mean_ratio = sum(ratios) / len(ratios)
    logger.info(f"Mean gene proportion (mapped/ref): {mean_ratio:.3f} ({len(ratios)} chromosome contigs)")
    return mean_ratio


# ----------------------------
# Reference-based debris classification
# ----------------------------
def classify_debris_and_unclassified(
    remaining_contigs: Set[str],
    query_fasta: Path,
    query_lengths: Dict[str, int],
    ref_fasta: Path,
    chrs_fasta: Optional[Path],
    qr_gene_count: Dict[Tuple[str, str], int],
    work_dir: Path,
    threads: int,
    min_coverage: float,
    min_protein_hits: int,
) -> Tuple[Set[str], Set[str], Dict[str, DebrisHit]]:
    """Classify remaining contigs as debris or unclassified based on reference homology.

    MOTIVATION: After chromosome debris detection (which catches assembly artifacts that
    are near-identical to assembled chromosomes), some contigs may still have homology to
    the reference genome without being classified. These include:
    - Divergent haplotigs that don't align well at nucleotide level to assembled chromosomes
      but retain protein-coding content (detected via miniprot hits)
    - Sequences from genomic regions not represented in the primary chromosome contigs
      but present in the reference (detected via reference nucleotide alignment)

    This reference-based approach complements chromosome debris detection by catching
    sequences with ancestral/evolutionary homology rather than assembly-specific duplicates.

    Classification criteria:
    - Debris: contigs with 50%+ nucleotide alignment coverage to reference OR 2+ miniprot hits
    - Unclassified: everything else (potential novel sequences, contaminants missed by
      earlier screening, or highly divergent sequences)

    Args:
        remaining_contigs: Set of contigs not yet classified
        query_fasta: Path to query FASTA
        query_lengths: Dict of contig name -> length
        ref_fasta: Path to reference FASTA
        chrs_fasta: Optional path to assembled chromosome FASTA (currently unused)
        qr_gene_count: Dict of (contig, ref_id) -> gene count from miniprot
        work_dir: Working directory
        threads: Number of threads
        min_coverage: Min alignment coverage fraction for debris classification
        min_protein_hits: Min miniprot protein hits for debris classification

    Returns:
        Tuple of (debris_contigs, unclassified_contigs, debris_hits)
        where debris_hits maps contig name to DebrisHit with coverage and identity details
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    debris_contigs: Set[str] = set()
    unclassified_contigs: Set[str] = set()
    debris_hits: Dict[str, DebrisHit] = {}

    # Check protein hits for each remaining contig
    contigs_with_protein_hits: Dict[str, int] = {}  # contig -> protein hit count
    for (contig, _), gene_count in qr_gene_count.items():
        if contig in remaining_contigs and gene_count >= min_protein_hits:
            # Track max hits across all ref chroms
            if contig not in contigs_with_protein_hits or gene_count > contigs_with_protein_hits[contig]:
                contigs_with_protein_hits[contig] = gene_count

    # For contigs without enough protein hits, check nucleotide alignment
    contigs_needing_alignment = remaining_contigs - set(contigs_with_protein_hits.keys())

    # Track coverage/identity from alignments
    contig_coverage: Dict[str, float] = {}
    contig_identity: Dict[str, float] = {}

    if contigs_needing_alignment:
        # Write subset of contigs to align
        subset_fasta = work_dir / "debris_candidates.fa"
        query_seqs = read_fasta_sequences(query_fasta)
        subset_seqs = {k: v for k, v in query_seqs.items() if k in contigs_needing_alignment}
        if subset_seqs:
            write_fasta(subset_seqs, subset_fasta)

            # Align against reference using minimap2/mm2plus if available
            if get_minimap2_exe() and subset_fasta.exists():
                paf_out = work_dir / "debris_vs_ref.paf"

                run_minimap2(
                    ref=ref_fasta,
                    qry=subset_fasta,
                    paf_out=paf_out,
                    threads=threads,
                    preset="asm20",
                    extra_args=["-c"],  # output CIGAR for identity calculation
                )

                # Parse PAF to get coverage and identity
                if paf_out.exists():
                    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
                    query_matches: Dict[str, int] = defaultdict(int)
                    query_alnlen: Dict[str, int] = defaultdict(int)

                    with paf_out.open("r") as fh:
                        for line in fh:
                            fields = line.rstrip("\n").split("\t")
                            if len(fields) < 12:
                                continue
                            qname = fields[0]
                            try:
                                qs = int(fields[2])
                                qe = int(fields[3])
                                matches = int(fields[9])
                                aln_len = int(fields[10])
                            except ValueError:
                                continue
                            if qs > qe:
                                qs, qe = qe, qs
                            query_intervals[qname].append((qs, qe))
                            query_matches[qname] += matches
                            query_alnlen[qname] += aln_len

                    for qname, intervals in query_intervals.items():
                        _, total_bp = merge_intervals(intervals)
                        qlen = query_lengths.get(qname, 0)
                        coverage = (total_bp / qlen) if qlen > 0 else 0.0
                        identity = (query_matches[qname] / query_alnlen[qname]) if query_alnlen[qname] > 0 else 0.0

                        contig_coverage[qname] = coverage
                        contig_identity[qname] = identity

                        if coverage >= min_coverage:
                            debris_contigs.add(qname)
                            debris_hits[qname] = DebrisHit(
                                coverage=coverage,
                                identity=identity,
                                protein_hits=0,
                                source="reference",
                            )

    # Add contigs with protein hits to debris
    for contig, protein_hits in contigs_with_protein_hits.items():
        debris_contigs.add(contig)
        debris_hits[contig] = DebrisHit(
            coverage=contig_coverage.get(contig, 0.0),
            identity=contig_identity.get(contig, 0.0),
            protein_hits=protein_hits,
            source="reference",
        )

    # Everything else is unclassified
    unclassified_contigs = remaining_contigs - debris_contigs

    logger.info(f"Debris contigs: {len(debris_contigs)}")
    logger.info(f"Unclassified contigs: {len(unclassified_contigs)}")

    return debris_contigs, unclassified_contigs, debris_hits


# ----------------------------
# Full-length vs fragment classification
# ----------------------------
def classify_full_length(
    ref_coverage: float,
    has_5p_telomere: bool,
    has_3p_telomere: bool,
    full_length_threshold: float = 0.70,
    length_ratio: Optional[float] = None,
    length_ratio_threshold: float = 0.80,
) -> Tuple[bool, str]:
    """Classify a contig as full-length or fragment based on telomeres, coverage,
    and contig-to-reference length ratio.

    Classification logic (in priority order):
    1. Both telomeres present -> full_length (high confidence)
    2. Length ratio >= threshold + one telomere -> full_length (medium confidence)
    3. Length ratio >= threshold + no telomeres -> full_length (low confidence)
    4. One telomere + coverage >= threshold -> full_length (medium confidence)
    5. One telomere + coverage < threshold -> fragment (medium confidence)
    6. No telomeres + coverage >= threshold -> full_length (low confidence)
    7. No telomeres + coverage < threshold -> fragment (low confidence)

    The length ratio check (steps 2-3) catches contigs that are clearly
    chromosome-scale but have low alignment-based ref_coverage due to
    cross-species divergence reducing alignment coverage.

    Args:
        ref_coverage: Fraction of reference chromosome spanned (0.0-1.0)
        has_5p_telomere: Telomere detected at 5' end
        has_3p_telomere: Telomere detected at 3' end
        full_length_threshold: Reference coverage threshold for full-length (default 0.70)
        length_ratio: Contig length / reference chromosome length (optional)
        length_ratio_threshold: Length ratio above which the contig is considered
            chromosome-scale regardless of alignment coverage (default 0.80)

    Returns:
        Tuple of (is_full_length, confidence)
    """
    n_telomeres = int(has_5p_telomere) + int(has_3p_telomere)

    if n_telomeres == 2:
        # Both telomeres = assembly-complete, high confidence
        return True, "high"

    # Length ratio rescue: contig is chromosome-scale even if alignment
    # coverage is low (e.g. cross-species with ~60% nucleotide identity)
    if length_ratio is not None and length_ratio >= length_ratio_threshold:
        if n_telomeres == 1:
            return True, "medium"
        else:
            return True, "low"

    if n_telomeres == 1:
        # One telomere - coverage determines classification
        if ref_coverage >= full_length_threshold:
            return True, "medium"
        else:
            return False, "medium"
    else:
        # No telomeres - coverage-based with low confidence
        if ref_coverage >= full_length_threshold:
            return True, "low"
        else:
            return False, "low"


# ----------------------------
# Rearrangement hypothesis detection
# ----------------------------
def compute_largest_cluster_metrics_per_ref(
    macro_block_rows: List[Tuple],
    max_gap: int = 500000,
) -> Dict[Tuple[str, str], Tuple[int, int, int]]:
    """Compute the largest cluster metrics for each (contig, ref_id) pair.

    Clusters chains to the same ref_id by proximity (ignoring intervening chains
    to other refs), then returns the span, aligned_bp, and gene count of the
    largest cluster.

    This approach:
    - Groups same-target chains that are close together (gap < max_gap)
    - Ignores what's between them (small chains to other targets don't matter)
    - Separates distant same-target chains into different clusters
    - Returns metrics for the largest cluster only

    Args:
        macro_block_rows: List of tuples from chain parsing:
            (contig, contig_len, ref_id, chrom_id, subgenome, strand, chain_id,
             qstart, qend, qspan_bp, union_bp, ..., gene_count_chain[16], ...)
        max_gap: Maximum gap (bp) between consecutive same-target chains to
            cluster them together. Default 500kb.

    Returns:
        Dict mapping (contig, ref_id) -> (span_bp, aligned_bp, gene_count)
        of largest cluster
    """
    # Group chains by (contig, ref_id)
    chains_by_key: Dict[Tuple[str, str], List[Tuple[int, int, int, int]]] = defaultdict(list)

    for row in macro_block_rows:
        if len(row) < 17:
            continue
        contig = row[0]
        ref_id = row[2]
        qstart = int(row[7])
        qend = int(row[8])
        union_bp = int(row[10])  # union_bp is index 10
        gene_count = int(row[16])  # gene_count_chain is index 16
        chains_by_key[(contig, ref_id)].append((qstart, qend, union_bp, gene_count))

    results: Dict[Tuple[str, str], Tuple[int, int, int]] = {}

    for key, chains in chains_by_key.items():
        if not chains:
            continue

        # Sort chains by qstart
        chains.sort(key=lambda x: x[0])

        # Cluster chains by gap between consecutive same-target chains
        clusters: List[List[Tuple[int, int, int, int]]] = []
        current_cluster = [chains[0]]

        for chain in chains[1:]:
            # Gap from end of current cluster to start of this chain
            cluster_end = max(c[1] for c in current_cluster)
            gap = chain[0] - cluster_end

            if gap <= max_gap:
                current_cluster.append(chain)
            else:
                clusters.append(current_cluster)
                current_cluster = [chain]

        clusters.append(current_cluster)

        # Find largest cluster by span
        best_span = 0
        best_aligned = 0
        best_genes = 0
        for cluster in clusters:
            cluster_start = min(c[0] for c in cluster)
            cluster_end = max(c[1] for c in cluster)
            span = cluster_end - cluster_start
            aligned = sum(c[2] for c in cluster)
            genes = sum(c[3] for c in cluster)
            if span > best_span:
                best_span = span
                best_aligned = aligned
                best_genes = genes

        results[key] = (best_span, best_aligned, best_genes)

    return results


def compute_reference_gene_density(
    qr_cluster_metrics: Dict[Tuple[str, str], Tuple[int, int, int]],
    best_ref: Dict[str, str],
) -> float:
    """Compute median on-target gene density (genes/Mb) across assigned contigs.

    For each contig with a best reference assignment, computes the gene density
    of its on-target cluster. Returns the median across all assigned contigs.

    This provides a genome-calibrated reference for evaluating off-target clusters
    in rearrangement detection.  Gene density varies across genomes (gene-rich vs
    gene-poor species) and within genomes (euchromatic arms vs pericentromeres),
    so an auto-calibrated threshold adapts better than a fixed gene count.

    Args:
        qr_cluster_metrics: Dict mapping (contig, ref_id) ->
            (span_bp, aligned_bp, gene_count) from compute_largest_cluster_metrics_per_ref
        best_ref: Dict of contig -> best reference chromosome assignment

    Returns:
        Median on-target gene density in genes per Mb.
        Returns 0.0 if no assigned contigs have cluster data.
    """
    densities: List[float] = []
    for contig, ref_id in best_ref.items():
        metrics = qr_cluster_metrics.get((contig, ref_id))
        if not metrics:
            continue
        span_bp, _, gene_count = metrics
        if span_bp > 0 and gene_count > 0:
            densities.append(gene_count / (span_bp / 1_000_000))

    if not densities:
        return 0.0

    densities.sort()
    n = len(densities)
    if n % 2 == 1:
        return densities[n // 2]
    return (densities[n // 2 - 1] + densities[n // 2]) / 2


def detect_rearrangement_candidates(
    contig: str,
    assigned_ref_id: str,
    contig_len: int,
    qr_cluster_metrics: Dict[Tuple[str, str], Tuple[int, int, int]],
    contig_refs: Set[str],
    threshold: float = 0.10,
    min_density: float = 0.15,
    synteny_mode: str = "protein",
    ref_gene_density: float = 0.0,
    density_frac: float = 0.10,
    min_genes_floor: int = 5,
) -> Optional[str]:
    """Detect potential chromosome rearrangements based on off-target alignment clusters.

    A contig is flagged for rearrangement if its largest cluster of alignments to an
    off-target chromosome has:
    - Span ≥ threshold fraction of the contig length
    - AND one of:
      - Nucleotide mode: density (aligned_bp / span) ≥ min_density
      - Protein mode: gene density ≥ density_frac × median on-target gene density,
        AND gene count ≥ min_genes_floor

    In protein mode, aligned_bp density is inherently low because it represents
    individual protein alignment footprints.  Instead, gene density (genes/Mb) is
    compared to the genome-wide median on-target density, scaled by density_frac.
    This auto-calibrates to the genome: gene-rich genomes get a higher threshold,
    gene-poor genomes get a lower one.  A low density_frac (default 10%) ensures
    that gene-poor regions (pericentromeres) are still detected, since even gene-poor
    chromosomal content has much higher density than scattered paralogous noise.

    Args:
        contig: Contig name
        assigned_ref_id: The reference chromosome this contig is assigned to
        contig_len: Length of the contig in bp
        qr_cluster_metrics: Dict mapping (contig, ref_id) ->
            (span_bp, aligned_bp, gene_count)
        contig_refs: Set of reference chromosomes this contig has alignments to
        threshold: Minimum span fraction to flag as rearrangement (default 0.10)
        min_density: Minimum alignment density within span for nucleotide mode
            (default 0.15)
        synteny_mode: "protein" or "nucleotide"
        ref_gene_density: Median on-target gene density in genes/Mb, from
            compute_reference_gene_density().  Used in protein mode.
        density_frac: Minimum fraction of ref_gene_density required for
            off-target clusters in protein mode (default 0.10 = 10%).
        min_genes_floor: Absolute minimum gene count to prevent single-hit
            artifacts in protein mode (default 5).

    Returns:
        Comma-separated list of off-target chromosomes meeting criteria,
        ordered by span fraction (highest first).
        Returns None if no rearrangements detected.
    """
    if not assigned_ref_id or contig_len <= 0:
        return None

    candidates = []
    for ref_id in contig_refs:
        if ref_id == assigned_ref_id:
            continue
        metrics = qr_cluster_metrics.get((contig, ref_id))
        if not metrics:
            continue
        span_bp, aligned_bp, gene_count = metrics
        if span_bp <= 0:
            continue

        span_fraction = span_bp / contig_len

        if span_fraction < threshold:
            continue

        if synteny_mode == "protein":
            # Auto-calibrated gene density threshold
            cluster_density = gene_count / (span_bp / 1_000_000)
            min_gene_density = ref_gene_density * density_frac
            if cluster_density >= min_gene_density and gene_count >= min_genes_floor:
                candidates.append((span_fraction, ref_id))
        else:
            # Nucleotide mode: use alignment density (blocks are dense)
            density = aligned_bp / span_bp
            if density >= min_density:
                candidates.append((span_fraction, ref_id))

    if not candidates:
        return None

    candidates.sort(key=lambda x: -x[0])
    return ",".join(ref_id for _, ref_id in candidates)


# ----------------------------
# 1D Gaussian Mixture Model
# ----------------------------
_SQRT_2PI = math.sqrt(2 * math.pi)


def _norm_pdf(x: float, mu: float, sigma: float) -> float:
    """Gaussian PDF with numerical guards for degenerate cases.

    Returns large/small sentinel values instead of inf/0 to keep
    downstream sum/log arithmetic well-behaved.
    """
    if sigma < 1e-12:
        return 1e300 if abs(x - mu) < 1e-12 else 1e-300
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * _SQRT_2PI)


def _gmm_1d_em(
    data: List[float],
    k: int,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> Tuple[List[float], List[float], List[float], float]:
    """Fit a 1D Gaussian Mixture Model using Expectation-Maximisation.

    Args:
        data: List of observed values.
        k: Number of mixture components.
        max_iter: Maximum EM iterations.
        tol: Convergence tolerance on log-likelihood.

    Returns:
        (weights, means, stds, log_likelihood) for each component.
    """
    n = len(data)
    if n == 0 or k < 1:
        return [], [], [], float("-inf")

    sorted_data = sorted(data)

    # Initialise means by quantile-based seeding (deterministic)
    means = [sorted_data[int((i + 0.5) * n / k)] for i in range(k)]
    stds = [max(1e-8, (sorted_data[-1] - sorted_data[0]) / (2 * k))] * k
    weights = [1.0 / k] * k

    prev_ll = float("-inf")

    for _iteration in range(max_iter):
        # E-step: compute responsibilities
        resp = []
        log_likelihood = 0.0
        for x in data:
            probs = [weights[j] * _norm_pdf(x, means[j], stds[j]) for j in range(k)]
            total = sum(probs)
            if total < 1e-300:
                total = 1e-300
            log_likelihood += math.log(total)
            resp.append([p / total for p in probs])

        # Check convergence
        if abs(log_likelihood - prev_ll) < tol:
            break
        prev_ll = log_likelihood

        # M-step: update parameters
        for j in range(k):
            nj = sum(resp[i][j] for i in range(n))
            if nj < 1e-10:
                continue
            weights[j] = nj / n
            means[j] = sum(resp[i][j] * data[i] for i in range(n)) / nj
            var_j = sum(resp[i][j] * (data[i] - means[j]) ** 2 for i in range(n)) / nj
            stds[j] = max(1e-8, math.sqrt(var_j))

    return weights, means, stds, log_likelihood


def _gmm_1d_bic(
    data: List[float],
    max_k: int = 3,
) -> Tuple[int, List[float], List[float], List[float], List[int]]:
    """Select the best number of components for a 1D GMM using BIC.

    Tests k=1..max_k and returns the model with the lowest BIC.

    Args:
        data: List of observed values.
        max_k: Maximum number of components to try.

    Returns:
        (best_k, weights, means, stds, labels) where labels[i] is the
        cluster index (0-based) for data[i].
    """
    n = len(data)
    if n < 2:
        return 1, [1.0], [data[0] if data else 0.0], [1e-8], [0] * n

    best_bic = float("inf")
    best_result: Tuple[int, List[float], List[float], List[float]] = (
        1, [1.0], [0.0], [1.0]
    )

    for k in range(1, min(max_k, n) + 1):
        weights, means, stds, ll = _gmm_1d_em(data, k)
        if not weights:
            continue
        num_params = 3 * k - 1
        bic = -2 * ll + num_params * math.log(n)
        if bic < best_bic:
            best_bic = bic
            best_result = (k, weights, means, stds)

    best_k, weights, means, stds = best_result

    # Assign each data point to the most likely component
    labels = []
    for x in data:
        probs = [weights[j] * _norm_pdf(x, means[j], stds[j]) for j in range(best_k)]
        labels.append(max(range(best_k), key=lambda j: probs[j]))

    return best_k, weights, means, stds, labels


def _select_subgenome_k(
    all_idents: List[float],
    ref_to_clfs: Dict[str, List['ContigClassification']],
    ref_ids: List[str],
    qr_best_chain_ident: Dict[Tuple[str, str], float],
    min_consistency: float = 0.70,
) -> Tuple[int, List[float], List[float]]:
    """Select the number of query subgenomes using GMM + paired validation.

    Strategy:
    1. Fit GMMs for k=1..3 and pick best by BIC as a starting candidate.
    2. If BIC picks k=1, also try k=2: validate the k=2 split using
       per-chromosome consistency (for chromosomes with exactly 2 copies,
       what fraction have their copies landing in different clusters?).
       Accept k=2 if consistency >= min_consistency.

    Paired validation is the sole gatekeeper — no minimum separation
    threshold is imposed.  If 70%+ of chromosomes consistently place
    their copies in different clusters, the split is accepted regardless
    of how small the identity gap is.  This avoids rejecting valid
    segmentations for closely related ancestral subgenomes.

    Args:
        all_idents: All identity values from multi-copy chromosomes.
        ref_to_clfs: Dict mapping ref_id -> list of ContigClassifications.
        ref_ids: Reference chromosome IDs in this subgenome group.
        qr_best_chain_ident: Dict of (contig, ref_id) -> best chain identity.
        min_consistency: Minimum fraction of 2-copy chromosomes that must
            have their copies in different clusters (default 0.70).

    Returns:
        (best_k, sorted_means, sorted_stds) with means sorted descending.
    """
    def _assign_labels(data, means_sorted):
        """Assign each value to the nearest cluster mean."""
        k = len(means_sorted)
        return [min(range(k), key=lambda j: abs(x - means_sorted[j])) for x in data]

    def _validate_k_by_pairing(sorted_means, k):
        """Check if a k-cluster split is consistent across chromosomes.

        For chromosomes with exactly k copies, verify that each copy
        lands in a different cluster (i.e., the clusters correspond to
        distinct subgenomes, not random variation).
        """
        if k < 2:
            return True

        n_testable = 0
        n_consistent = 0
        for ref_id in ref_ids:
            ref_clfs = ref_to_clfs[ref_id]
            if len(ref_clfs) != k:
                continue
            n_testable += 1
            copy_idents = [
                qr_best_chain_ident.get((clf.original_name, ref_id), 0.0)
                for clf in ref_clfs
            ]
            labels = _assign_labels(copy_idents, sorted_means)
            if len(set(labels)) == k:
                n_consistent += 1

        if n_testable < 3:
            return False
        return (n_consistent / n_testable) >= min_consistency

    def _sort_components(means, stds, k):
        """Sort GMM components by mean descending (primary = highest identity)."""
        order = sorted(range(k), key=lambda j: means[j], reverse=True)
        return [means[j] for j in order], [stds[j] for j in order]

    # Fit all candidate GMMs once (avoids redundant k=2 refit in rescue path)
    fits: Dict[int, Tuple[List[float], List[float], List[float], float]] = {}
    for k in range(1, min(4, len(all_idents) + 1)):
        w, m, s, ll = _gmm_1d_em(all_idents, k)
        if w:
            fits[k] = (w, m, s, ll)

    # BIC selection
    best_bic = float("inf")
    bic_k = 1
    for k, (w, m, s, ll) in fits.items():
        num_params = 3 * k - 1
        bic = -2 * ll + num_params * math.log(len(all_idents))
        if bic < best_bic:
            best_bic = bic
            bic_k = k

    # Validate BIC result via paired consistency
    if bic_k > 1 and bic_k in fits:
        _, bic_means, bic_stds, _ = fits[bic_k]
        sorted_means, sorted_stds = _sort_components(bic_means, bic_stds, bic_k)
        if _validate_k_by_pairing(sorted_means, bic_k):
            return bic_k, sorted_means, sorted_stds

    # Rescue: try k=2 with paired validation (BIC may be too conservative)
    if 2 in fits:
        _, means_2, stds_2, _ = fits[2]
        sorted_means, sorted_stds = _sort_components(means_2, stds_2, 2)
        if _validate_k_by_pairing(sorted_means, 2):
            return 2, sorted_means, sorted_stds

    # Fall back to k=1
    mean_val = statistics.mean(all_idents)
    std_val = statistics.stdev(all_idents) if len(all_idents) > 1 else 0.0
    return 1, [mean_val], [std_val]


# ----------------------------
# Query subgenome inference
# ----------------------------
def infer_query_subgenomes(
    classifications: List[ContigClassification],
    qr_best_chain_ident: Dict[Tuple[str, str], float],
    subgenome_k: float = 1.0,
) -> None:
    """Infer query subgenomes when multiple contigs map to the same reference.

    Uses a two-stage approach:
    1. **Global GMM + paired validation**: Collect alignment identities from
       all multi-copy reference chromosomes within each reference subgenome,
       fit a 1-D Gaussian Mixture Model.  Model selection uses BIC as a
       first pass, then validates (or rescues) via per-chromosome pairing:
       for chromosomes with exactly k copies, each copy should land in a
       different cluster.  This leverages the biological structure (paired
       subgenomic copies) to detect subgenomes even when the marginal
       identity distributions overlap substantially.
    2. **Per-chromosome assignment**: Use the global cluster means to assign
       each contig to a query subgenome.  Single-copy chromosomes are left
       unsuffixed (group 1).

    Args:
        classifications: List of ContigClassification objects (modified in place)
        qr_best_chain_ident: Dict of (contig, ref_id) -> best chain identity
        subgenome_k: Not used by GMM (kept for API compatibility)
    """
    # Collect identities for all chromosome-assigned contigs
    chrom_clfs = [c for c in classifications if c.classification == "chrom_assigned" and c.assigned_ref_id]

    if not chrom_clfs:
        return

    # Set seq_identity_vs_ref for all contigs
    for clf in chrom_clfs:
        ident = qr_best_chain_ident.get((clf.original_name, clf.assigned_ref_id), 0.0)
        if ident > 0:
            clf.seq_identity_vs_ref = ident

    # Group contigs by reference chromosome
    ref_to_clfs: Dict[str, List[ContigClassification]] = defaultdict(list)
    for clf in chrom_clfs:
        ref_to_clfs[clf.assigned_ref_id].append(clf)

    # Identify multi-copy reference chromosomes (those needing subgenome splitting)
    multi_copy_refs = {rid for rid, clfs in ref_to_clfs.items() if len(clfs) > 1}

    # Group multi-copy chromosomes by reference subgenome (A, T, P, …)
    ref_sg_to_multi_refs: Dict[Optional[str], List[str]] = defaultdict(list)
    for ref_id in multi_copy_refs:
        _, ref_sg = split_chrom_subgenome(ref_id)
        ref_sg_to_multi_refs[ref_sg].append(ref_id)

    # Per reference subgenome: run global GMM on all multi-copy identities
    # Maps ref_subgenome -> (best_k, ordered_cluster_means)
    # Cluster ordering: sorted by mean identity descending so cluster 0 = primary
    ref_sg_gmm: Dict[Optional[str], Tuple[int, List[float]]] = {}

    for ref_sg, ref_ids in ref_sg_to_multi_refs.items():
        # Collect all identities from multi-copy chromosomes in this ref subgenome
        all_idents = []
        for ref_id in ref_ids:
            for clf in ref_to_clfs[ref_id]:
                ident = qr_best_chain_ident.get((clf.original_name, ref_id), 0.0)
                if ident > 0:
                    all_idents.append(ident)

        if len(all_idents) < 4:
            # Too few data points for meaningful clustering
            ref_sg_gmm[ref_sg] = (1, [])
            logger.info(
                f"Subgenome clustering ({ref_sg or 'default'}): "
                f"n={len(all_idents)}, too few for GMM, skipping"
            )
            continue

        # Fit candidate GMMs (k=1..3)
        best_k, sorted_means, sorted_stds = _select_subgenome_k(
            all_idents, ref_to_clfs, ref_ids, qr_best_chain_ident,
        )

        ref_sg_gmm[ref_sg] = (best_k, sorted_means)

        if best_k == 1:
            logger.info(
                f"Subgenome clustering ({ref_sg or 'default'}): "
                f"n={len(all_idents)}, best_k=1, "
                f"mean={sorted_means[0]:.4f}, std={sorted_stds[0]:.4f} — "
                f"single subgenome"
            )
        else:
            cluster_desc = ", ".join(
                f"k{i+1}: mean={sorted_means[i]:.4f} std={sorted_stds[i]:.4f}"
                for i in range(best_k)
            )
            logger.info(
                f"Subgenome clustering ({ref_sg or 'default'}): "
                f"n={len(all_idents)}, best_k={best_k}, {cluster_desc}"
            )

    # Assign per-chromosome subgenome labels
    for ref_id, ref_clfs in ref_to_clfs.items():
        if len(ref_clfs) <= 1:
            # Single contig — no subgenome suffix
            ref_clfs[0].query_subgenome = None
            ref_clfs[0].query_subgenome_grp = 1
            continue

        _, ref_sg = split_chrom_subgenome(ref_id)
        gmm_k, sorted_means = ref_sg_gmm.get(ref_sg, (1, []))

        if gmm_k <= 1 or not sorted_means:
            # No subgenome split detected — assign all to group 1
            for clf in ref_clfs:
                clf.query_subgenome = None
                clf.query_subgenome_grp = 1
            continue

        # Sort contigs by identity descending
        ref_clfs.sort(
            key=lambda c: qr_best_chain_ident.get((c.original_name, ref_id), 0.0),
            reverse=True,
        )

        # Assign each contig to the nearest GMM cluster (by mean)
        assignments: List[int] = []  # 0-based cluster per contig
        for clf in ref_clfs:
            ident = qr_best_chain_ident.get((clf.original_name, ref_id), 0.0)
            best_cluster = min(
                range(gmm_k),
                key=lambda j: abs(ident - sorted_means[j]),
            )
            assignments.append(best_cluster)

        # Rescue: if this chromosome has exactly k copies but they all
        # landed in the same cluster, rank-assign instead (highest
        # identity → primary, next → B, etc.).  This is safe because
        # the global split is already validated across chromosomes;
        # these are borderline cases where the per-chromosome gap is
        # smaller than the distance between cluster means.
        if len(ref_clfs) == gmm_k and len(set(assignments)) == 1:
            # ref_clfs is sorted by identity descending above
            assignments = list(range(gmm_k))

        for clf, cluster in zip(ref_clfs, assignments):
            clf.query_subgenome_grp = cluster + 1  # 1-based
            if cluster == 0:
                clf.query_subgenome = None  # Primary — no suffix
            else:
                clf.query_subgenome = chr(ord('A') + cluster)  # B, C, …


# ----------------------------
# Contig naming
# ----------------------------
def generate_contig_names(
    classifications: List[ContigClassification],
    query_lengths: Dict[str, int],
    add_subgenome_suffix: Optional[str],
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """Generate new contig names with full-length/fragment distinction.

    Naming scheme: chr<ref>(_<query_subgenome>)?(_c<copy>|_f<frag>)?

    Suffixes compose left-to-right; subgenome comes before copy/fragment:
    - chr1A:           Full-length chr1A, single copy, no resolved subgenome
    - chr1A_B:         Full-length chr1A, query subgenome B
    - chr1A_c1/c2:     Multiple full-length copies of the same chromosome (unusual)
    - chr1A_f1/f2:     Chromosome fragments, ordered by descending length
    - chr1A_B_f1:      Longest fragment of chr1A from query subgenome B

    Non-chromosome contigs are named contig_1, contig_2, … by descending length.

    Args:
        classifications: List of ContigClassification with is_full_length and query_subgenome set
        query_lengths: Dict of contig name -> length
        add_subgenome_suffix: Optional suffix to add if ref lacks subgenome
        ref_norm_to_orig: Optional mapping of normalized -> original ref IDs

    Returns:
        Dict mapping original contig name -> new name
    """
    name_mapping: Dict[str, str] = {}

    # Build lookup by original name
    clf_by_name = {clf.original_name: clf for clf in classifications}

    # Group chromosome contigs by (assigned_ref_id, query_subgenome)
    groups: Dict[Tuple[str, Optional[str]], List[ContigClassification]] = defaultdict(list)
    non_chrom_contigs: List[str] = []

    for clf in classifications:
        if clf.classification == "chrom_assigned" and clf.assigned_ref_id:
            key = (clf.assigned_ref_id, clf.query_subgenome)
            groups[key].append(clf)
        else:
            non_chrom_contigs.append(clf.original_name)

    # Name chromosome contigs
    for (ref_id, query_sub), contigs in groups.items():
        base_name = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
        # Add subgenome suffix if reference doesn't have one and user requested it
        if add_subgenome_suffix:
            chrom_id, sub = split_chrom_subgenome(base_name)
            if sub == "NA":
                base_name = f"{chrom_id}{add_subgenome_suffix}"

        # Add query subgenome suffix if present
        if query_sub:
            base_name = f"{base_name}_{query_sub}"

        # Separate full-length and fragments
        full_length = [c for c in contigs if c.is_full_length]
        fragments = [c for c in contigs if not c.is_full_length]

        # Sort each group by length descending
        full_length.sort(key=lambda c: c.contig_len, reverse=True)
        fragments.sort(key=lambda c: c.contig_len, reverse=True)

        # Name full-length contigs
        if len(full_length) == 1:
            name_mapping[full_length[0].original_name] = base_name
        elif len(full_length) > 1:
            # Multiple full-length copies - use _c1, _c2
            for i, clf in enumerate(full_length, start=1):
                name_mapping[clf.original_name] = f"{base_name}_c{i}"

        # Name fragments - use _f1, _f2
        for i, clf in enumerate(fragments, start=1):
            name_mapping[clf.original_name] = f"{base_name}_f{i}"

    # Name non-chromosome contigs
    # Sort by length descending
    non_chrom_contigs.sort(key=lambda c: query_lengths.get(c, 0), reverse=True)
    for i, contig in enumerate(non_chrom_contigs, start=1):
        name_mapping[contig] = f"contig_{i}"

    return name_mapping


# ----------------------------
# Classification pipeline
# ----------------------------
def classify_all_contigs(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    best_ref: Dict[str, str],
    chr_like_minlen: int,
    ev: ChainEvidenceResult,
    ref_gene_counts: Dict[str, int],
    chrC_contig: Optional[str],
    chrM_contig: Optional[str],
    organelle_debris: Set[str],
    rdna_contigs: Set[str],
    contaminants: Dict[str, ContaminantHit],
    chromosome_debris: Set[str],
    other_debris: Set[str],
    add_subgenome_suffix: Optional[str],
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
    query_gc: Optional[Dict[str, float]] = None,
    ref_gc_mean: Optional[float] = None,
    ref_gc_std: Optional[float] = None,
    asm_gc_mean: Optional[float] = None,
    asm_gc_std: Optional[float] = None,
    organelle_hits: Optional[Dict[str, OrganelleHit]] = None,
    rdna_hits: Optional[Dict[str, RdnaHit]] = None,
    chrom_debris_hits: Optional[Dict[str, DebrisHit]] = None,
    other_debris_hits: Optional[Dict[str, DebrisHit]] = None,
    # New parameters for full-length/fragment and subgenome inference
    ref_lengths: Optional[Dict[str, int]] = None,
    telomere_results: Optional[Dict[str, 'TelomereResult']] = None,
    full_length_threshold: float = 0.70,
    subgenome_k: float = 1.0,
    rearrangement_threshold: float = 0.10,
    rearrangement_density_frac: float = 0.10,
    synteny_mode: str = "protein",
) -> List[ContigClassification]:
    """Classify all contigs and assign confidence levels.

    This function classifies each contig into one of nine categories and assigns
    a confidence level (high/medium/low) based on multiple evidence sources.

    Classification categories:
    - chrom_assigned: Chromosome-length contigs with synteny support
    - chrom_unassigned: Chromosome-length contigs without synteny
    - organelle_complete: Complete organelle genomes (chrC, chrM)
    - organelle_debris: Partial organelle sequences
    - rDNA: Ribosomal DNA repeat units
    - contaminant: Sequences from contaminating organisms
    - chrom_debris: High-coverage duplicates of assembled chromosomes
    - debris: Assembly fragments with reference homology
    - unclassified: Sequences with no classification evidence

    Confidence levels are mode-dependent:

    **Protein mode** (default):
    - Gene proportion: Fraction of reference genes aligned to contig
    - GC deviation: Standard deviations from GC baseline
    - Protein hits: Number of miniprot alignments

    **Nucleotide mode**:
    - Reference coverage: Fraction of reference chromosome spanned by synteny blocks
    - Chain identity: Alignment identity of best chain
    - GC deviation: Same as protein mode

    Common to both modes:
    - GC deviation baselines:
      - chrom_assigned: compared to reference nuclear GC
      - other categories: compared to assembly chromosome GC
    - Coverage: Alignment coverage fraction (0.0-1.0)
    - Identity: Alignment identity (0.0-1.0)

    Args:
        query_fasta: Path to query FASTA
        query_lengths: Dict of contig name -> length
        best_ref: Dict of contig -> best reference assignment
        chr_like_minlen: Minimum length for chromosome classification
        ev: Chain evidence results
        ref_gene_counts: Dict of ref_id -> gene count
        chrC_contig: Chloroplast contig name
        chrM_contig: Mitochondrial contig name
        organelle_debris: Set of organelle debris contigs
        rdna_contigs: Set of rDNA contigs
        contaminants: Dict of contig -> ContaminantHit with taxid, name, coverage, score
        chromosome_debris: Set of chromosome debris contigs (from chr-vs-chr alignment)
        other_debris: Set of other debris contigs (from ref/protein alignment)
        add_subgenome_suffix: Optional subgenome suffix to add
        query_gc: Dict of contig name -> GC content (0.0-1.0)
        ref_gc_mean: Mean GC content of reference nuclear chromosomes
        ref_gc_std: Standard deviation of reference GC content
        asm_gc_mean: Mean GC content of assembly chromosome contigs
        asm_gc_std: Standard deviation of assembly chromosome GC content
        organelle_hits: Dict of contig -> OrganelleHit with detection details
        rdna_hits: Dict of contig -> RdnaHit with detection details
        chrom_debris_hits: Dict of contig -> DebrisHit from chromosome debris detection
        other_debris_hits: Dict of contig -> DebrisHit from reference/protein detection
        synteny_mode: "protein" or "nucleotide" — controls confidence criteria

    Returns:
        List of ContigClassification objects with assigned names, categories, and confidence levels
    """

    classifications: List[ContigClassification] = []
    classified_contigs: Set[str] = set()

    # Compute cluster metrics per (contig, ref_id) for rearrangement detection
    qr_cluster_metrics = compute_largest_cluster_metrics_per_ref(ev.macro_block_rows)

    # Compute genome-wide reference gene density for auto-calibrated
    # rearrangement detection in protein mode
    ref_gene_density = compute_reference_gene_density(qr_cluster_metrics, best_ref)

    # Helper to compute GC deviation vs REFERENCE (for chrom_assigned validation)
    def _gc_deviation_vs_ref(contig: str) -> Optional[float]:
        if query_gc is None or ref_gc_mean is None or ref_gc_std is None:
            return None
        gc = query_gc.get(contig)
        if gc is None:
            return None
        if ref_gc_std == 0 or ref_gc_std < 0.001:
            return None  # Cannot compute meaningful deviation with zero/near-zero std
        return abs(gc - ref_gc_mean) / ref_gc_std

    # Helper to compute GC deviation vs ASSEMBLY chromosomes (for non-chrom classifications)
    # This is more appropriate for divergent genomes where the assembly may have
    # different GC than the reference
    def _gc_deviation_vs_asm(contig: str) -> Optional[float]:
        if query_gc is None or asm_gc_mean is None or asm_gc_std is None:
            return None
        gc = query_gc.get(contig)
        if gc is None:
            return None
        if asm_gc_std == 0 or asm_gc_std < 0.001:
            return None  # Cannot compute meaningful deviation with zero/near-zero std
        return abs(gc - asm_gc_mean) / asm_gc_std

    # Helper to get contig GC content
    def _get_gc(contig: str) -> Optional[float]:
        return query_gc.get(contig) if query_gc else None

    # 1. Chromosome-length contigs with reference assignment
    for contig, ref_id in best_ref.items():
        if not ref_id:
            continue
        contig_len = query_lengths.get(contig, 0)
        if contig_len < chr_like_minlen:
            continue

        # Check if it passed synteny gates
        gene_count = ev.qr_gene_count.get((contig, ref_id), 0)
        ref_total_genes = ref_gene_counts.get(ref_id, 0)
        gene_proportion = (gene_count / ref_total_genes) if ref_total_genes > 0 else None
        gc_dev = _gc_deviation_vs_ref(contig)  # Use reference GC for chrom_assigned validation

        # Compute synteny score (0-1): gene proportion in protein mode, ref coverage in nucleotide mode
        if synteny_mode == "nucleotide":
            # Compute ref_cov early so we can use it for synteny_score
            _ref_cov_for_score = None
            if ref_lengths and ev.qr_ref_span_bp:
                _ref_span = ev.qr_ref_span_bp.get((contig, ref_id), 0)
                _ref_len = ref_lengths.get(ref_id, 0)
                if _ref_len > 0:
                    _ref_cov_for_score = _ref_span / _ref_len
            synteny_score = min(1.0, _ref_cov_for_score) if _ref_cov_for_score else 0.0
        else:
            synteny_score = min(1.0, gene_proportion * 2) if gene_proportion else 0.0

        # Compute reference coverage (span-based, robust to divergence)
        ref_cov = None
        if ref_lengths and ev.qr_ref_span_bp:
            ref_span = ev.qr_ref_span_bp.get((contig, ref_id), 0)
            ref_len = ref_lengths.get(ref_id, 0)
            if ref_len > 0:
                ref_cov = ref_span / ref_len

        # Get telomere detection results
        has_5p = False
        has_3p = False
        if telomere_results and contig in telomere_results:
            telo = telomere_results[contig]
            has_5p = telo.has_5p_telomere
            has_3p = telo.has_3p_telomere

        # Classify as full-length or fragment
        is_full = False
        fl_confidence = "low"
        if ref_cov is not None:
            # Compute length ratio for cross-species divergence rescue
            _ref_len = ref_lengths.get(ref_id, 0) if ref_lengths else 0
            _length_ratio = (contig_len / _ref_len) if _ref_len > 0 else None
            is_full, fl_confidence = classify_full_length(
                ref_coverage=ref_cov,
                has_5p_telomere=has_5p,
                has_3p_telomere=has_3p,
                full_length_threshold=full_length_threshold,
                length_ratio=_length_ratio,
            )
        else:
            # No ref_lengths available - default to fragment with low confidence
            is_full = False
            fl_confidence = "low"

        # Compute confidence (mode-dependent criteria)
        if synteny_mode == "nucleotide":
            # Nucleotide mode: use ref_coverage and chain identity instead of gene_proportion
            confidence = "high"
            if ref_cov is None or ref_cov < 0.1:
                confidence = "low"
            elif ref_cov < 0.3:
                confidence = "medium"
            # Identity penalty: low identity suggests spurious assignment
            best_ident = ev.qr_best_chain_ident.get((contig, ref_id), 0.0)
            if best_ident > 0 and best_ident < 0.5:
                confidence = "low"
        else:
            # Protein mode: use gene proportion
            # High: good synteny (>20% genes) AND GC similar to reference (<2 std)
            # Medium: moderate synteny OR borderline GC
            # Low: weak synteny AND/OR very different GC
            confidence = "high"
            if gene_proportion is None or gene_proportion < 0.1:
                confidence = "low"
            elif gene_proportion < 0.2:
                confidence = "medium"
        # GC deviation penalty (both modes)
        if gc_dev is not None and gc_dev > 3.0:
            confidence = "low"
        elif gc_dev is not None and gc_dev > 2.0 and confidence == "high":
            confidence = "medium"

        # Collinearity-based confidence adjustment
        collin_score = None
        if ev.qr_collinearity:
            collin_score = ev.qr_collinearity.get((contig, ref_id))
        n_chains = ev.qr_nchains_kept.get((contig, ref_id), 0)
        if collin_score is not None and n_chains >= 3:
            if collin_score >= 0.9 and n_chains >= 5 and confidence == "medium":
                confidence = "high"
            elif collin_score < 0.3 and n_chains >= 10:
                confidence = "low"
            elif collin_score < 0.5 and n_chains >= 5:
                if confidence == "high":
                    confidence = "medium"
                elif confidence == "medium":
                    confidence = "low"

        # Detect potential rearrangements (largest off-target cluster ≥ threshold)
        rearrangement_candidates = detect_rearrangement_candidates(
            contig=contig,
            assigned_ref_id=ref_id,
            contig_len=contig_len,
            qr_cluster_metrics=qr_cluster_metrics,
            contig_refs=ev.contig_refs.get(contig, set()),
            threshold=rearrangement_threshold,
            synteny_mode=synteny_mode,
            ref_gene_density=ref_gene_density,
            density_frac=rearrangement_density_frac,
        )

        classifications.append(ContigClassification(
            original_name=contig,
            new_name="",  # Will be filled later
            classification="chrom_assigned",
            reversed=False,  # Will be filled later
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id=ref_id,
            ref_gene_proportion=gene_proportion,
            contig_len=contig_len,
            gc_content=_get_gc(contig),
            gc_deviation=gc_dev,
            synteny_score=synteny_score,
            classification_confidence=confidence,
            # Full-length vs fragment fields
            ref_coverage=ref_cov,
            is_full_length=is_full,
            full_length_confidence=fl_confidence,
            has_5p_telomere=has_5p,
            has_3p_telomere=has_3p,
            # Collinearity score
            collinearity_score=collin_score,
            # Rearrangement hypothesis
            rearrangement_candidates=rearrangement_candidates,
        ))
        classified_contigs.add(contig)

    # 1b. Chromosome-length contigs WITHOUT reference assignment
    # (but not contaminants - those are handled separately)
    for contig, contig_len in query_lengths.items():
        if contig in classified_contigs:
            continue
        if contig_len < chr_like_minlen:
            continue
        # Skip known contaminants - they should be classified as contaminants, not chrom_unassigned
        if contig in contaminants:
            continue

        gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline

        # Confidence is inherently low/medium for unassigned contigs
        # Medium: GC similar to assembly chromosomes (might be a real chromosome from a divergent region)
        # Low: GC different from assembly (might be contaminant or unusual sequence)
        confidence = "medium"
        if gc_dev is not None and gc_dev > 2.0:
            confidence = "low"

        # This contig is chromosome-length but has no synteny assignment
        classifications.append(ContigClassification(
            original_name=contig,
            new_name="",  # Will be filled later
            classification="chrom_unassigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id=None,
            ref_gene_proportion=None,
            contig_len=contig_len,
            gc_content=_get_gc(contig),
            gc_deviation=gc_dev,
            synteny_score=0.0,  # No synteny evidence
            classification_confidence=confidence,
        ))
        classified_contigs.add(contig)

    # 2. Organelle contigs
    # organelle_complete: High confidence - passed strict detection criteria
    # (coverage >80%, length within tolerance)
    if chrC_contig and chrC_contig not in classified_contigs:
        hit = organelle_hits.get(chrC_contig) if organelle_hits else None
        # High confidence if coverage >= 90%, medium if >= 80% (our threshold)
        confidence = "high"
        if hit and hit.coverage < 0.90:
            confidence = "medium"
        classifications.append(ContigClassification(
            original_name=chrC_contig,
            new_name="chrC",
            classification="organelle_complete",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chrC",
            ref_gene_proportion=None,
            contig_len=query_lengths.get(chrC_contig, 0),
            gc_content=_get_gc(chrC_contig),
            gc_deviation=None,  # Organelles have different GC, don't compare to nuclear
            classification_confidence=confidence,
        ))
        classified_contigs.add(chrC_contig)

    if chrM_contig and chrM_contig not in classified_contigs:
        hit = organelle_hits.get(chrM_contig) if organelle_hits else None
        # High confidence if coverage >= 90%, medium if >= 80% (our threshold)
        confidence = "high"
        if hit and hit.coverage < 0.90:
            confidence = "medium"
        classifications.append(ContigClassification(
            original_name=chrM_contig,
            new_name="chrM",
            classification="organelle_complete",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chrM",
            ref_gene_proportion=None,
            contig_len=query_lengths.get(chrM_contig, 0),
            gc_content=_get_gc(chrM_contig),
            gc_deviation=None,  # Organelles have different GC, don't compare to nuclear
            classification_confidence=confidence,
        ))
        classified_contigs.add(chrM_contig)

    # Organelle debris: Medium confidence - partial organelle match
    for contig in organelle_debris:
        if contig not in classified_contigs:
            hit = organelle_hits.get(contig) if organelle_hits else None
            # Medium confidence by default (partial match)
            # Low confidence if coverage is very low (<60%)
            confidence = "medium"
            if hit and hit.coverage < 0.60:
                confidence = "low"
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="organelle_debris",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=hit.organelle_type if hit else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=None,  # Organelles have different GC
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 3. rDNA contigs
    # rDNA: Confidence based on coverage and identity
    for contig in rdna_contigs:
        if contig not in classified_contigs:
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            hit = rdna_hits.get(contig) if rdna_hits else None
            # High confidence if coverage >= 80% and identity >= 95%
            # Medium confidence if coverage >= 50% (our threshold)
            # Low confidence if borderline
            confidence = "medium"
            if hit:
                if hit.coverage >= 0.80 and hit.identity >= 0.95:
                    confidence = "high"
                elif hit.coverage < 0.60:
                    confidence = "low"
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="rDNA",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 4. Contaminants
    for contig, hit in contaminants.items():
        if contig not in classified_contigs:
            contig_len = query_lengths.get(contig, 0)
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline

            # Compute confidence based on coverage and GC deviation
            # High confidence: high coverage (>0.8) AND (GC deviation > 2 std from assembly)
            # Medium confidence: moderate coverage (0.5-0.8) OR conflicting evidence
            # Low confidence: low coverage (<0.5) AND some synteny evidence
            confidence = "high"
            if hit.coverage < 0.5:
                confidence = "low"
            elif hit.coverage < 0.8:
                confidence = "medium"
            # GC deviation supports contaminant classification
            if gc_dev is not None and gc_dev > 2.0:
                confidence = "high"

            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="contaminant",
                reversed=False,
                contaminant_taxid=hit.taxid,
                contaminant_sci=hit.sci_name,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=contig_len,
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                contam_score=min(1.0, hit.score / 1e9) if hit.score else None,  # Normalize to 0-1 range
                contam_coverage=hit.coverage,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 5. Chromosome debris (from chr-vs-chr alignment detection)
    # Confidence based on coverage, identity, and GC deviation
    for contig in chromosome_debris:
        if contig not in classified_contigs:
            ref_id = best_ref.get(contig, "")
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            hit = chrom_debris_hits.get(contig) if chrom_debris_hits else None
            # High confidence: ≥90% coverage AND ≥95% identity
            # Medium confidence: passed strict 80% cov, 90% identity threshold
            # Low confidence: borderline or GC very different from assembly chromosomes
            confidence = "high"
            if hit:
                if hit.coverage >= 0.90 and hit.identity >= 0.95:
                    confidence = "high"
                elif hit.coverage < 0.85 or hit.identity < 0.92:
                    confidence = "medium"
            if gc_dev is not None and gc_dev > 3.0:
                confidence = "medium" if confidence == "high" else "low"
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="chrom_debris",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=ref_id if ref_id else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 6. Other debris (from reference/protein alignment detection)
    # Confidence based on coverage, identity, protein hits, and GC deviation
    for contig in other_debris:
        if contig not in classified_contigs:
            # Check if this contig has synteny support
            ref_id = best_ref.get(contig, "")
            classification = "chrom_debris" if ref_id else "debris"
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            hit = other_debris_hits.get(contig) if other_debris_hits else None

            # Confidence based on detection evidence (mode-dependent)
            confidence = "medium"
            if hit:
                if synteny_mode == "nucleotide":
                    # Nucleotide mode: protein_hits is always 0, use coverage/identity only
                    if hit.coverage >= 0.80:
                        confidence = "high"
                    elif hit.coverage < 0.55:
                        confidence = "low"
                else:
                    # Protein mode: use coverage and protein hits
                    # High: ≥80% coverage OR ≥5 protein hits
                    # Medium: ≥50% coverage OR ≥2 protein hits (our threshold)
                    # Low: borderline
                    if hit.coverage >= 0.80 or hit.protein_hits >= 5:
                        confidence = "high"
                    elif hit.coverage < 0.55 and hit.protein_hits < 3:
                        confidence = "low"

            # GC deviation penalty
            if gc_dev is not None and gc_dev > 3.0 and not ref_id:
                confidence = "low"

            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification=classification,
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=ref_id if ref_id else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence=confidence,
            ))
            classified_contigs.add(contig)

    # 7. Unclassified (everything else)
    # Low confidence by definition - no evidence for any classification
    for contig in query_lengths.keys():
        if contig not in classified_contigs:
            gc_dev = _gc_deviation_vs_asm(contig)  # Use assembly chromosome GC baseline
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="unclassified",
                reversed=False,
                contaminant_taxid=None,
                contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
                gc_content=_get_gc(contig),
                gc_deviation=gc_dev,
                classification_confidence="low",
            ))

    # Infer query subgenomes when multiple contigs map to same reference
    if ev.qr_best_chain_ident:
        infer_query_subgenomes(
            classifications=classifications,
            qr_best_chain_ident=ev.qr_best_chain_ident,
            subgenome_k=subgenome_k,
        )

    # Generate names
    name_mapping = generate_contig_names(
        classifications,
        query_lengths,
        add_subgenome_suffix,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    # Update classifications with names
    for clf in classifications:
        if not clf.new_name:
            clf.new_name = name_mapping.get(clf.original_name, clf.original_name)

    return classifications
