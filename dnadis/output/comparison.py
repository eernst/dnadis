#!/usr/bin/env python3
"""Cross-assembly comparison outputs for multi-assembly mode.

Provides functions to build per-assembly summary metrics and write
comparison TSVs after all assemblies have been processed.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from dnadis.utils.reference_utils import is_nuclear_chromosome
from dnadis.models import (
    AssemblyResult,
    ChromRefSummary,
    ChainEvidenceResult,
    CompleasmResult,
    ContigClassification,
    DepthStats,
    RdnaArray,
)


# Per-reference-chromosome assembly states (best → missing).  Mirrors the R
# helper compute_chrom_contiguity() in reports/common.R; keep the two in sync.
#
# State is derived from *anchors* (chrom_assigned) only, grouped into copies of
# the chromosome by ``ref_copy`` (overlapping anchors are separate copies,
# tiling anchors are pieces of one copy).  Each copy is T2T (one gap-free
# contig with telomeres at both ends), gapped (one contig with telomeres at
# both ends and scaffold gaps inside), single (one contig) or split (two or
# more contigs);
# the chromosome takes its least contiguous copy's state, and the number of
# copies is reported separately.  Anchors without a ref_copy (no alignment to
# the reference, e.g. placed by sister anchoring) cannot be placed in a copy
# and are left out, unless no anchor has one: then they are taken as one copy.  Fragment contigs (chrom_fragment) do not
# change the state: they are surfaced as a separate count / badge, since they
# are frequently divergent duplicates or eccDNA rather than evidence that the
# chromosome itself is broken.
CHROM_STATE_LEVELS = [
    "T2T", "gapped", "single", "split", "fragments_only", "absent",
]


def _chrom_contiguity_state(
    anchors: List[ContigClassification], n_fragment: int
) -> tuple[str, int]:
    """Return one reference chromosome's (state, number of copies).

    ``anchors`` are its chrom_assigned contigs; ``n_fragment`` counts its
    chrom_fragment contigs.
    """
    if not anchors:
        return ("fragments_only" if n_fragment > 0 else "absent"), 0
    copies: Dict[Optional[int], List[ContigClassification]] = defaultdict(list)
    placed = [c for c in anchors if c.ref_copy is not None]
    for c in placed or anchors:
        copies[c.ref_copy].append(c)
    order = {s: i for i, s in enumerate(CHROM_STATE_LEVELS)}
    state = "T2T"
    for members in copies.values():
        if len(members) > 1:
            copy_state = "split"
        elif members[0].has_5p_telomere and members[0].has_3p_telomere:
            copy_state = "gapped" if members[0].n_gaps else "T2T"
        else:
            copy_state = "single"
        state = max(state, copy_state, key=order.__getitem__)
    return state, len(copies)


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


def _gap_free_pieces(
    lengths: Dict[str, int], gaps: Dict[str, List[Tuple[int, int]]]
) -> List[int]:
    """Lengths of the pieces left after cutting every sequence at its gaps."""
    pieces: List[int] = []
    for name, length in lengths.items():
        pos = 0
        for start, end in gaps.get(name, ()):
            if start > pos:
                pieces.append(start - pos)
            pos = end
        if length > pos:
            pieces.append(length - pos)
    return pieces


def build_assembly_result(
    assembly_name: str,
    assembly_path: Path,
    outprefix: Path,
    classifications: List[ContigClassification],
    qry_lengths: Dict[str, int],
    ref_lengths_norm: Dict[str, int],
    ev: ChainEvidenceResult,
    cobionts_filtered: Dict,
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
    cobionts_tsv_path: Optional[Path] = None,
    rdna_annotations_tsv: Optional[Path] = None,
    rdna_arrays_tsv: Optional[Path] = None,
    agp_tsv: Optional[Path] = None,
    chrs_fasta: Optional[Path] = None,
    per_subgenome_chrs: Optional[Dict[str, Path]] = None,
    compleasm_full: Optional[CompleasmResult] = None,
    compleasm_chrs: Optional[CompleasmResult] = None,
    compleasm_non_chrs: Optional[CompleasmResult] = None,
    synteny_mode: str = "nucleotide",
    qry_gaps: Optional[Dict[str, List[Tuple[int, int]]]] = None,
    organelle_ref_cov: Optional[Dict[str, float]] = None,
) -> AssemblyResult:
    """Build an AssemblyResult summarizing one assembly's analysis.

    Args:
        assembly_name: Short name for this assembly.
        assembly_path: Path to the query FASTA.
        outprefix: Output prefix for this assembly.
        classifications: List of ContigClassification objects.
        qry_lengths: Dict mapping contig name -> length.
        ref_lengths_norm: Normalized reference chromosome lengths.
        ev: ChainEvidenceResult from synteny analysis.
        cobionts_filtered: Dict of filtered cobiont hits.
        chrC_contig: Name of chrC contig, or None.
        chrM_contig: Name of chrM contig, or None.
        organelle_ref_cov: Organelle -> reference coverage of its complete
            contig or partial set (OrganelleCalls.set_ref_coverage).
        rdna_arrays: List of detected rDNA arrays.
        depth_stats: Dict mapping contig name -> DepthStats.
        chimera_primary_frac: Primary fraction threshold for chimera detection.
        chimera_secondary_frac: Secondary fraction threshold for chimera detection.
        summary_tsv: Path to contig_summary.tsv.
        segments_tsv: Path to segments.tsv.
        evidence_tsv: Path to evidence_summary.tsv.
        macro_blocks_tsv: Path to macro_blocks.tsv.
        cobionts_tsv_path: Path to cobionts.tsv, or None.
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
    qry_gaps = qry_gaps or {}
    n_gaps = sum(len(runs) for runs in qry_gaps.values())
    gap_bp = sum(e - s for runs in qry_gaps.values() for s, e in runs)
    contig_n50, contig_l50 = _compute_n50_l50(
        _gap_free_pieces(qry_lengths, qry_gaps)) if qry_gaps else (n50, l50)

    # --- Classification counts/bp ---
    classification_counts: Dict[str, int] = defaultdict(int)
    classification_bp: Dict[str, int] = defaultdict(int)
    for clf in classifications:
        classification_counts[clf.classification] += 1
        if clf.organelle_unit is not None:
            # One unit is written as the organelle and the whole concatemer to
            # the debris FASTA; count the unit as organelle and the rest as
            # debris so class totals still sum to the assembly size
            start, end, _ = clf.organelle_unit
            classification_bp[clf.classification] += end - start
            classification_bp["debris"] += clf.contig_len - (end - start)
        else:
            classification_bp[clf.classification] += clf.contig_len

    # --- Chromosome completeness ---
    chrom_assigned = [c for c in classifications if c.classification == "chrom_assigned"]
    chrom_unassigned = [c for c in classifications if c.classification == "chrom_unassigned"]

    n_chromosome_scale = sum(1 for c in chrom_assigned if c.is_chromosome_scale)
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

    # Aligned-bp-weighted identity over chrom_assigned contigs:
    #   sum(seq_identity_vs_ref * best_ref_union_bp) / sum(contig_len)
    # Unmapped query bp contribute 0, so the metric reflects both alignment
    # breadth and identity rather than identity within aligned regions only.
    #
    # Only defined in nucleotide mode.  In protein mode, best_ref_union_bp is
    # the bp under protein hits (coding exons), so best_bp/contig_len collapses
    # to gene density and the metric tracks inverse genome size rather than
    # alignment quality.  Left as None there; mean_identity carries the
    # protein-mode quality signal and the comparison sort falls back to it.
    weighted_identity: Optional[float] = None
    if synteny_mode != "protein":
        weighted_ident_num = 0.0
        weighted_ident_den = 0
        for c in chrom_assigned:
            if c.seq_identity_vs_ref is None:
                continue
            union_bp = int(ev.best_bp.get(c.original_name, 0) or 0)
            weighted_ident_num += c.seq_identity_vs_ref * union_bp
            weighted_ident_den += c.contig_len
        weighted_identity = (
            weighted_ident_num / weighted_ident_den if weighted_ident_den > 0 else None
        )

    collinearities = [c.collinearity_score for c in chrom_assigned if c.collinearity_score is not None]
    mean_collinearity = sum(collinearities) / len(collinearities) if collinearities else None

    gc_devs = [abs(c.gc_deviation) for c in chrom_assigned if c.gc_deviation is not None]
    mean_gc_deviation = sum(gc_devs) / len(gc_devs) if gc_devs else None

    # --- Cobionts ---
    cobiont_contigs = [c for c in classifications if c.classification == "cobiont"]
    n_cobionts = len(cobiont_contigs)
    total_cobiont_bp = sum(c.contig_len for c in cobiont_contigs)
    cobiont_species = set()
    for c in cobiont_contigs:
        if c.cobiont_sci:
            cobiont_species.add(c.cobiont_sci)
    n_unique_cobiont_species = len(cobiont_species)

    # --- rDNA ---
    rdna_contigs = [c for c in classifications if c.classification == "rDNA"]
    n_rdna_contigs = len(rdna_contigs)
    total_rdna_bp = sum(c.contig_len for c in rdna_contigs)
    n_rdna_arrays_count = len(rdna_arrays) if rdna_arrays else 0

    # --- Organelles ---
    organelle_found: Dict[str, str] = {}
    organelle_n_contigs: Dict[str, int] = {}
    for organelle, contig in (("chrC", chrC_contig), ("chrM", chrM_contig)):
        n_partial = sum(
            1 for c in classifications
            if c.classification == "organelle_partial" and c.assigned_ref_id == organelle
        )
        if contig is not None:
            organelle_found[organelle], organelle_n_contigs[organelle] = "yes", 1
        elif n_partial:
            organelle_found[organelle], organelle_n_contigs[organelle] = "partial", n_partial
        else:
            organelle_found[organelle], organelle_n_contigs[organelle] = "no", 0
    organelle_ref_cov = organelle_ref_cov or {}

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
    # Group chrom_assigned by assigned_ref_id (drives the Overview pill heatmap:
    # n_contigs is the copy-number indicator, quality fields the identity signal).
    ref_groups: Dict[str, List[ContigClassification]] = defaultdict(list)
    for c in chrom_assigned:
        if c.assigned_ref_id:
            ref_groups[c.assigned_ref_id].append(c)

    # Group chrom_fragment by assigned_ref_id (for the contiguity / state view).
    frag_groups: Dict[str, List[ContigClassification]] = defaultdict(list)
    for c in classifications:
        if c.classification == "chrom_fragment" and c.assigned_ref_id:
            frag_groups[c.assigned_ref_id].append(c)

    for ref_id, ref_len in ref_lengths_norm.items():
        # Per-chromosome states (and the n_chr_* counts) cover nuclear chromosomes
        # only; unplaced scaffolds and organelles are not chromosomes to complete.
        if not is_nuclear_chromosome(ref_id):
            continue
        contigs_for_ref = ref_groups.get(ref_id, [])
        n_contigs_ref = len(contigs_for_ref)
        total_assigned_bp_ref = sum(c.contig_len for c in contigs_for_ref)

        best_cov = None
        any_chromosome_scale = False
        any_both_telo = False
        ident_vals = []

        for c in contigs_for_ref:
            if c.ref_coverage is not None:
                if best_cov is None or c.ref_coverage > best_cov:
                    best_cov = c.ref_coverage
            if c.is_chromosome_scale:
                any_chromosome_scale = True
            if c.has_5p_telomere and c.has_3p_telomere:
                any_both_telo = True
            if c.seq_identity_vs_ref is not None:
                ident_vals.append(c.seq_identity_vs_ref)

        # Contiguity state from anchors (chrom_assigned) vs fragments.
        n_fragment = len(frag_groups.get(ref_id, []))
        state, n_copies = _chrom_contiguity_state(contigs_for_ref, n_fragment)

        chrom_ref_coverage[ref_id] = ChromRefSummary(
            ref_id=ref_id,
            ref_length=ref_len,
            n_contigs=n_contigs_ref,
            total_assigned_bp=total_assigned_bp_ref,
            best_ref_coverage=best_cov,
            is_chromosome_scale=any_chromosome_scale,
            has_both_telomeres=any_both_telo,
            mean_identity=sum(ident_vals) / len(ident_vals) if ident_vals else None,
            n_fragment=n_fragment,
            n_copies=n_copies,
            state=state,
            n_gaps=sum(c.n_gaps for c in contigs_for_ref),
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
        n_chromosome_scale=n_chromosome_scale,
        n_with_both_telomeres=n_with_both_telomeres,
        n_with_any_telomere=n_with_any_telomere,
        mean_ref_coverage=mean_ref_coverage,
        n_chimeric=n_chimeric,
        mean_identity=mean_identity,
        weighted_identity=weighted_identity,
        mean_collinearity=mean_collinearity,
        mean_gc_deviation=mean_gc_deviation,
        n_cobionts=n_cobionts,
        total_cobiont_bp=total_cobiont_bp,
        n_unique_cobiont_species=n_unique_cobiont_species,
        n_rdna_contigs=n_rdna_contigs,
        total_rdna_bp=total_rdna_bp,
        n_rdna_arrays=n_rdna_arrays_count,
        chrC_found=organelle_found["chrC"],
        chrM_found=organelle_found["chrM"],
        chrC_n_contigs=organelle_n_contigs["chrC"],
        chrM_n_contigs=organelle_n_contigs["chrM"],
        chrC_ref_cov=organelle_ref_cov.get("chrC") if organelle_found["chrC"] != "no" else None,
        chrM_ref_cov=organelle_ref_cov.get("chrM") if organelle_found["chrM"] != "no" else None,
        mean_chrom_depth=mean_chrom_depth,
        chrom_ref_coverage=chrom_ref_coverage,
        classifications=classifications,
        summary_tsv=summary_tsv,
        segments_tsv=segments_tsv,
        evidence_tsv=evidence_tsv,
        macro_blocks_tsv=macro_blocks_tsv,
        cobionts_tsv=cobionts_tsv_path,
        rdna_annotations_tsv=rdna_annotations_tsv,
        rdna_arrays_tsv=rdna_arrays_tsv,
        agp_tsv=agp_tsv,
        chrs_fasta=chrs_fasta,
        per_subgenome_chrs=per_subgenome_chrs or {},
        compleasm_full=compleasm_full,
        compleasm_chrs=compleasm_chrs,
        compleasm_non_chrs=compleasm_non_chrs,
        n_gaps=n_gaps,
        gap_bp=gap_bp,
        contig_n50=contig_n50,
        contig_l50=contig_l50,
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
        "contig_n50",
        "contig_l50",
        "n_gaps",
        "gap_bp",
        "largest_contig",
        "n_chrom_assigned",
        "n_chrom_unassigned",
        "chrom_assigned_bp",
        "n_chrom_fragment",
        "chrom_fragment_bp",
        "n_circular_element",
        "circular_element_bp",
        "n_chr_t2t",
        "n_chr_gapped",
        "n_chr_single",
        "n_chr_split",
        "n_chr_multi_copy",
        "n_chr_with_fragments",
        "n_chr_fragments_only",
        "n_chr_absent",
        "n_chromosome_scale",
        "n_both_telomeres",
        "n_any_telomere",
        "mean_ref_coverage",
        "n_chimeric",
        "mean_identity",
        "weighted_identity",
        "mean_collinearity",
        "mean_gc_deviation",
        "n_cobionts",
        "cobiont_bp",
        "n_cobiont_species",
        "n_rdna",
        "rdna_bp",
        "n_rdna_arrays",
        "chrC_found",
        "chrM_found",
        "chrC_n_contigs",
        "chrM_n_contigs",
        "chrC_ref_cov",
        "chrM_ref_cov",
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
            organelle_bp = sum(
                r.classification_bp.get(c, 0)
                for c in ("organelle_complete", "organelle_partial", "organelle_debris")
            )
            debris_bp = r.classification_bp.get("debris", 0) + r.classification_bp.get("chrom_debris", 0)
            unclassified_bp = r.classification_bp.get("unclassified", 0)

            # Per-chromosome state tally (excludes absent refs below chr scale is
            # already handled upstream; here we count every reference chromosome).
            state_counts = defaultdict(int)
            n_chr_with_fragments = 0
            n_chr_multi_copy = 0
            for crs in r.chrom_ref_coverage.values():
                state_counts[crs.state] += 1
                if crs.n_fragment > 0:
                    n_chr_with_fragments += 1
                if crs.n_copies > 1:
                    n_chr_multi_copy += 1

            row = [
                r.assembly_name,
                str(r.total_contigs),
                str(r.total_bp),
                str(r.n50),
                str(r.l50),
                str(r.contig_n50),
                str(r.contig_l50),
                str(r.n_gaps),
                str(r.gap_bp),
                str(r.largest_contig),
                str(r.n_chrom_assigned),
                str(r.n_chrom_unassigned),
                str(r.classification_bp.get("chrom_assigned", 0)),
                str(r.classification_counts.get("chrom_fragment", 0)),
                str(r.classification_bp.get("chrom_fragment", 0)),
                str(r.classification_counts.get("circular_element", 0)),
                str(r.classification_bp.get("circular_element", 0)),
                str(state_counts.get("T2T", 0)),
                str(state_counts.get("gapped", 0)),
                str(state_counts.get("single", 0)),
                str(state_counts.get("split", 0)),
                str(n_chr_multi_copy),
                str(n_chr_with_fragments),
                str(state_counts.get("fragments_only", 0)),
                str(state_counts.get("absent", 0)),
                str(r.n_chromosome_scale),
                str(r.n_with_both_telomeres),
                str(r.n_with_any_telomere),
                _fmt(r.mean_ref_coverage),
                str(r.n_chimeric),
                _fmt(r.mean_identity),
                _fmt(r.weighted_identity),
                _fmt(r.mean_collinearity),
                _fmt(r.mean_gc_deviation),
                str(r.n_cobionts),
                str(r.total_cobiont_bp),
                str(r.n_unique_cobiont_species),
                str(r.n_rdna_contigs),
                str(r.total_rdna_bp),
                str(r.n_rdna_arrays),
                r.chrC_found,
                r.chrM_found,
                str(r.chrC_n_contigs),
                str(r.chrM_n_contigs),
                _fmt(r.chrC_ref_cov),
                _fmt(r.chrM_ref_cov),
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
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> None:
    """Write chromosome completeness matrix in long format.

    One row per (ref_id, assembly) pair. Includes rows for assemblies
    that have no contigs for a given reference chromosome (n_contigs=0).

    Args:
        path: Output TSV path.
        results: List of AssemblyResult objects.
        ref_lengths_norm: Normalized reference chromosome lengths.
        ref_norm_to_orig: Map from normalized to original reference IDs. When
            given, ``ref_id`` is written with the reference's own spelling (e.g.
            ``Chr01``), matching ref_lengths.tsv and the contig summaries.
    """
    norm_to_orig = ref_norm_to_orig or {}
    columns = [
        "ref_id",
        "ref_length",
        "assembly",
        "n_contigs",
        "total_assigned_bp",
        "best_ref_coverage",
        "is_chromosome_scale",
        "has_both_telomeres",
        "mean_identity",
        "n_fragment",
        "n_copies",
        "state",
        "n_gaps",
    ]

    def _fmt(val, fmt_str=".4f"):
        if val is None:
            return ""
        return f"{val:{fmt_str}}"

    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")

        for ref_id in sorted(ref_lengths_norm.keys()):
            if not is_nuclear_chromosome(ref_id):
                continue  # skip unplaced scaffolds / organelles
            ref_len = ref_lengths_norm[ref_id]
            for r in results:
                crs = r.chrom_ref_coverage.get(ref_id)
                if crs:
                    row = [
                        norm_to_orig.get(ref_id, ref_id),
                        str(ref_len),
                        r.assembly_name,
                        str(crs.n_contigs),
                        str(crs.total_assigned_bp),
                        _fmt(crs.best_ref_coverage),
                        "yes" if crs.is_chromosome_scale else "no",
                        "yes" if crs.has_both_telomeres else "no",
                        _fmt(crs.mean_identity),
                        str(crs.n_fragment),
                        str(crs.n_copies),
                        crs.state,
                        str(crs.n_gaps),
                    ]
                else:
                    row = [
                        norm_to_orig.get(ref_id, ref_id),
                        str(ref_len),
                        r.assembly_name,
                        "0",
                        "0",
                        "",
                        "no",
                        "no",
                        "",
                        "0",
                        "0",
                        "absent",
                        "0",
                    ]
                fh.write("\t".join(row) + "\n")
