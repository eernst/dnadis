#!/usr/bin/env python3
"""Cross-assembly comparison outputs for multi-assembly mode.

Provides functions to build per-assembly summary metrics and write
comparison TSVs after all assemblies have been processed.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from final_finalizer.models import (
    AssemblyResult,
    ChromRefSummary,
    ChainEvidenceResult,
    CompleasmResult,
    ContigClassification,
    DepthStats,
    RdnaArray,
)


def _compute_n50_l50(lengths: List[int]) -> tuple[int, int]:
    """Compute N50 and L50 from a list of contig lengths.

    N50: length of the shortest contig such that contigs of that length
    or longer cover at least 50% of the total assembly size.
    L50: number of contigs of length >= N50.
    """
    if not lengths:
        return 0, 0
    sorted_lens = sorted(lengths, reverse=True)
    total = sum(sorted_lens)
    half = total / 2.0
    cumulative = 0
    for i, length in enumerate(sorted_lens):
        cumulative += length
        if cumulative >= half:
            return length, i + 1
    return sorted_lens[-1], len(sorted_lens)


def build_assembly_result(
    assembly_name: str,
    assembly_path: Path,
    outprefix: Path,
    classifications: List[ContigClassification],
    qry_lengths: Dict[str, int],
    ref_lengths_norm: Dict[str, int],
    ev: ChainEvidenceResult,
    contaminants_filtered: Dict,
    chrC_contig: Optional[str],
    chrM_contig: Optional[str],
    rdna_arrays: List[RdnaArray],
    depth_stats: Dict[str, DepthStats],
    chimera_primary_frac: float,
    chimera_secondary_frac: float,
    summary_tsv: Path,
    segments_tsv: Path,
    evidence_tsv: Path,
    macro_blocks_tsv: Path,
    contaminants_tsv_path: Optional[Path] = None,
    rdna_annotations_tsv: Optional[Path] = None,
    rdna_arrays_tsv: Optional[Path] = None,
    per_subgenome_chrs: Optional[Dict[str, Path]] = None,
    compleasm_chrs: Optional[CompleasmResult] = None,
    compleasm_non_chrs: Optional[CompleasmResult] = None,
) -> AssemblyResult:
    """Build an AssemblyResult summarizing one assembly's finalization.

    Args:
        assembly_name: Short name for this assembly.
        assembly_path: Path to the query FASTA.
        outprefix: Output prefix for this assembly.
        classifications: List of ContigClassification objects.
        qry_lengths: Dict mapping contig name -> length.
        ref_lengths_norm: Normalized reference chromosome lengths.
        ev: ChainEvidenceResult from synteny analysis.
        contaminants_filtered: Dict of filtered contaminant hits.
        chrC_contig: Name of chrC contig, or None.
        chrM_contig: Name of chrM contig, or None.
        rdna_arrays: List of detected rDNA arrays.
        depth_stats: Dict mapping contig name -> DepthStats.
        chimera_primary_frac: Primary fraction threshold for chimera detection.
        chimera_secondary_frac: Secondary fraction threshold for chimera detection.
        summary_tsv: Path to contig_summary.tsv.
        segments_tsv: Path to segments.tsv.
        evidence_tsv: Path to evidence_summary.tsv.
        macro_blocks_tsv: Path to macro_blocks.tsv.
        contaminants_tsv_path: Path to contaminants.tsv, or None.
        rdna_annotations_tsv: Path to rdna_annotations.tsv, or None.
        rdna_arrays_tsv: Path to rdna_arrays.tsv, or None.

    Returns:
        Populated AssemblyResult.
    """
    # --- Contiguity ---
    all_lengths = list(qry_lengths.values())
    total_contigs = len(all_lengths)
    total_bp = sum(all_lengths)
    n50, l50 = _compute_n50_l50(all_lengths)
    largest_contig = max(all_lengths) if all_lengths else 0

    # --- Classification counts/bp ---
    classification_counts: Dict[str, int] = defaultdict(int)
    classification_bp: Dict[str, int] = defaultdict(int)
    for clf in classifications:
        classification_counts[clf.classification] += 1
        classification_bp[clf.classification] += clf.contig_len

    # --- Chromosome completeness ---
    chrom_assigned = [c for c in classifications if c.classification == "chrom_assigned"]
    chrom_unassigned = [c for c in classifications if c.classification == "chrom_unassigned"]

    n_full_length = sum(1 for c in chrom_assigned if c.is_full_length)
    n_with_both_telomeres = sum(
        1 for c in chrom_assigned
        if c.has_5p_telomere and c.has_3p_telomere
    )
    n_with_any_telomere = sum(
        1 for c in chrom_assigned
        if c.has_5p_telomere or c.has_3p_telomere
    )

    ref_coverages = [c.ref_coverage for c in chrom_assigned if c.ref_coverage is not None]
    mean_ref_coverage = sum(ref_coverages) / len(ref_coverages) if ref_coverages else None

    # --- Chimera detection ---
    n_chimeric = 0
    for clf in chrom_assigned:
        q = clf.original_name
        total_aligned_bp = int(ev.contig_total.get(q, 0) or 0)
        n_ref_hits = len(ev.contig_refs.get(q, set()) or set())
        if total_aligned_bp > 0 and n_ref_hits > 1:
            best_union_bp = int(ev.best_bp.get(q, 0) or 0)
            second_union_bp = int(ev.second_bp.get(q, 0) or 0)
            primary_frac = best_union_bp / total_aligned_bp
            secondary_frac = second_union_bp / total_aligned_bp
            if primary_frac < chimera_primary_frac and secondary_frac >= chimera_secondary_frac:
                n_chimeric += 1

    # --- Quality metrics (chrom_assigned) ---
    identities = [c.seq_identity_vs_ref for c in chrom_assigned if c.seq_identity_vs_ref is not None]
    mean_identity = sum(identities) / len(identities) if identities else None

    collinearities = [c.collinearity_score for c in chrom_assigned if c.collinearity_score is not None]
    mean_collinearity = sum(collinearities) / len(collinearities) if collinearities else None

    gc_devs = [abs(c.gc_deviation) for c in chrom_assigned if c.gc_deviation is not None]
    mean_gc_deviation = sum(gc_devs) / len(gc_devs) if gc_devs else None

    # --- Contamination ---
    contaminant_contigs = [c for c in classifications if c.classification == "contaminant"]
    n_contaminants = len(contaminant_contigs)
    total_contaminant_bp = sum(c.contig_len for c in contaminant_contigs)
    contaminant_species = set()
    for c in contaminant_contigs:
        if c.contaminant_sci:
            contaminant_species.add(c.contaminant_sci)
    n_unique_contaminant_species = len(contaminant_species)

    # --- rDNA ---
    rdna_contigs = [c for c in classifications if c.classification == "rDNA"]
    n_rdna_contigs = len(rdna_contigs)
    total_rdna_bp = sum(c.contig_len for c in rdna_contigs)
    n_rdna_arrays_count = len(rdna_arrays) if rdna_arrays else 0

    # --- Organelles ---
    chrC_found = chrC_contig is not None
    chrM_found = chrM_contig is not None

    # --- Read depth ---
    mean_chrom_depth: Optional[float] = None
    if depth_stats:
        chrom_depths = []
        for c in chrom_assigned:
            ds = depth_stats.get(c.original_name)
            if ds:
                chrom_depths.append(ds.mean_depth)
        if chrom_depths:
            mean_chrom_depth = sum(chrom_depths) / len(chrom_depths)

    # --- Per-reference-chromosome detail ---
    chrom_ref_coverage: Dict[str, ChromRefSummary] = {}
    # Group chrom_assigned by assigned_ref_id
    ref_groups: Dict[str, List[ContigClassification]] = defaultdict(list)
    for c in chrom_assigned:
        if c.assigned_ref_id:
            ref_groups[c.assigned_ref_id].append(c)

    for ref_id, ref_len in ref_lengths_norm.items():
        contigs_for_ref = ref_groups.get(ref_id, [])
        n_contigs_ref = len(contigs_for_ref)
        total_assigned_bp_ref = sum(c.contig_len for c in contigs_for_ref)

        best_cov = None
        any_full_length = False
        any_both_telo = False
        ident_vals = []

        for c in contigs_for_ref:
            if c.ref_coverage is not None:
                if best_cov is None or c.ref_coverage > best_cov:
                    best_cov = c.ref_coverage
            if c.is_full_length:
                any_full_length = True
            if c.has_5p_telomere and c.has_3p_telomere:
                any_both_telo = True
            if c.seq_identity_vs_ref is not None:
                ident_vals.append(c.seq_identity_vs_ref)

        chrom_ref_coverage[ref_id] = ChromRefSummary(
            ref_id=ref_id,
            ref_length=ref_len,
            n_contigs=n_contigs_ref,
            total_assigned_bp=total_assigned_bp_ref,
            best_ref_coverage=best_cov,
            is_full_length=any_full_length,
            has_both_telomeres=any_both_telo,
            mean_identity=sum(ident_vals) / len(ident_vals) if ident_vals else None,
        )

    return AssemblyResult(
        assembly_name=assembly_name,
        assembly_path=assembly_path,
        outprefix=outprefix,
        total_contigs=total_contigs,
        total_bp=total_bp,
        n50=n50,
        l50=l50,
        largest_contig=largest_contig,
        classification_counts=dict(classification_counts),
        classification_bp=dict(classification_bp),
        n_chrom_assigned=len(chrom_assigned),
        n_chrom_unassigned=len(chrom_unassigned),
        n_full_length=n_full_length,
        n_with_both_telomeres=n_with_both_telomeres,
        n_with_any_telomere=n_with_any_telomere,
        mean_ref_coverage=mean_ref_coverage,
        n_chimeric=n_chimeric,
        mean_identity=mean_identity,
        mean_collinearity=mean_collinearity,
        mean_gc_deviation=mean_gc_deviation,
        n_contaminants=n_contaminants,
        total_contaminant_bp=total_contaminant_bp,
        n_unique_contaminant_species=n_unique_contaminant_species,
        n_rdna_contigs=n_rdna_contigs,
        total_rdna_bp=total_rdna_bp,
        n_rdna_arrays=n_rdna_arrays_count,
        chrC_found=chrC_found,
        chrM_found=chrM_found,
        mean_chrom_depth=mean_chrom_depth,
        chrom_ref_coverage=chrom_ref_coverage,
        classifications=classifications,
        summary_tsv=summary_tsv,
        segments_tsv=segments_tsv,
        evidence_tsv=evidence_tsv,
        macro_blocks_tsv=macro_blocks_tsv,
        contaminants_tsv=contaminants_tsv_path,
        rdna_annotations_tsv=rdna_annotations_tsv,
        rdna_arrays_tsv=rdna_arrays_tsv,
        per_subgenome_chrs=per_subgenome_chrs or {},
        compleasm_chrs=compleasm_chrs,
        compleasm_non_chrs=compleasm_non_chrs,
    )


def write_comparison_summary_tsv(
    path: Path,
    results: List[AssemblyResult],
) -> None:
    """Write cross-assembly comparison summary TSV.

    One row per assembly with key metrics for side-by-side comparison.

    Args:
        path: Output TSV path.
        results: List of AssemblyResult objects.
    """
    columns = [
        "assembly",
        "total_contigs",
        "total_bp",
        "n50",
        "l50",
        "largest_contig",
        "n_chrom_assigned",
        "n_chrom_unassigned",
        "chrom_assigned_bp",
        "n_full_length",
        "n_both_telomeres",
        "n_any_telomere",
        "mean_ref_coverage",
        "n_chimeric",
        "mean_identity",
        "mean_collinearity",
        "mean_gc_deviation",
        "n_contaminants",
        "contaminant_bp",
        "n_contaminant_species",
        "n_rdna",
        "rdna_bp",
        "n_rdna_arrays",
        "chrC_found",
        "chrM_found",
        "organelle_bp",
        "debris_bp",
        "unclassified_bp",
        "mean_chrom_depth",
        "compleasm_lineage",
        "compleasm_chrs_N",
        "compleasm_chrs_S",
        "compleasm_chrs_D",
        "compleasm_chrs_F",
        "compleasm_chrs_I",
        "compleasm_chrs_M",
        "compleasm_non_chrs_N",
        "compleasm_non_chrs_S",
        "compleasm_non_chrs_D",
        "compleasm_non_chrs_F",
        "compleasm_non_chrs_I",
        "compleasm_non_chrs_M",
    ]

    def _fmt(val, fmt_str=".4f"):
        if val is None:
            return ""
        return f"{val:{fmt_str}}"

    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")

        for r in results:
            organelle_bp = r.classification_bp.get("organelle_complete", 0) + r.classification_bp.get("organelle_debris", 0)
            debris_bp = r.classification_bp.get("debris", 0) + r.classification_bp.get("chrom_debris", 0)
            unclassified_bp = r.classification_bp.get("unclassified", 0)

            row = [
                r.assembly_name,
                str(r.total_contigs),
                str(r.total_bp),
                str(r.n50),
                str(r.l50),
                str(r.largest_contig),
                str(r.n_chrom_assigned),
                str(r.n_chrom_unassigned),
                str(r.classification_bp.get("chrom_assigned", 0)),
                str(r.n_full_length),
                str(r.n_with_both_telomeres),
                str(r.n_with_any_telomere),
                _fmt(r.mean_ref_coverage),
                str(r.n_chimeric),
                _fmt(r.mean_identity),
                _fmt(r.mean_collinearity),
                _fmt(r.mean_gc_deviation),
                str(r.n_contaminants),
                str(r.total_contaminant_bp),
                str(r.n_unique_contaminant_species),
                str(r.n_rdna_contigs),
                str(r.total_rdna_bp),
                str(r.n_rdna_arrays),
                "yes" if r.chrC_found else "no",
                "yes" if r.chrM_found else "no",
                str(organelle_bp),
                str(debris_bp),
                str(unclassified_bp),
                _fmt(r.mean_chrom_depth, ".2f"),
                # Compleasm columns
                r.compleasm_chrs.lineage if r.compleasm_chrs else "",
                *(r.compleasm_chrs.tsv_fields() if r.compleasm_chrs else [""] * 6),
                *(r.compleasm_non_chrs.tsv_fields() if r.compleasm_non_chrs else [""] * 6),
            ]
            fh.write("\t".join(row) + "\n")


def write_chromosome_completeness_tsv(
    path: Path,
    results: List[AssemblyResult],
    ref_lengths_norm: Dict[str, int],
) -> None:
    """Write chromosome completeness matrix in long format.

    One row per (ref_id, assembly) pair. Includes rows for assemblies
    that have no contigs for a given reference chromosome (n_contigs=0).

    Args:
        path: Output TSV path.
        results: List of AssemblyResult objects.
        ref_lengths_norm: Normalized reference chromosome lengths.
    """
    columns = [
        "ref_id",
        "ref_length",
        "assembly",
        "n_contigs",
        "total_assigned_bp",
        "best_ref_coverage",
        "is_full_length",
        "has_both_telomeres",
        "mean_identity",
    ]

    def _fmt(val, fmt_str=".4f"):
        if val is None:
            return ""
        return f"{val:{fmt_str}}"

    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")

        for ref_id in sorted(ref_lengths_norm.keys()):
            ref_len = ref_lengths_norm[ref_id]
            for r in results:
                crs = r.chrom_ref_coverage.get(ref_id)
                if crs:
                    row = [
                        ref_id,
                        str(ref_len),
                        r.assembly_name,
                        str(crs.n_contigs),
                        str(crs.total_assigned_bp),
                        _fmt(crs.best_ref_coverage),
                        "yes" if crs.is_full_length else "no",
                        "yes" if crs.has_both_telomeres else "no",
                        _fmt(crs.mean_identity),
                    ]
                else:
                    row = [
                        ref_id,
                        str(ref_len),
                        r.assembly_name,
                        "0",
                        "0",
                        "",
                        "no",
                        "no",
                        "",
                    ]
                fh.write("\t".join(row) + "\n")
