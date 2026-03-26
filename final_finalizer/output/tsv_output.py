#!/usr/bin/env python3
"""TSV output functions for contig summaries and chain data."""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.models import ContigClassification, ContaminantHitExtended, DepthStats, RearrangementCall, RdnaArray, RdnaLocus
from final_finalizer.utils.reference_utils import split_chrom_subgenome

# Use a large but finite value instead of infinity to avoid parsing issues
MAX_SCORE_RATIO = 1e9


def _denormalize_ref_id(ref_id: str, ref_norm_to_orig: Optional[Dict[str, str]]) -> str:
    if not ref_norm_to_orig:
        return ref_id
    return ref_norm_to_orig.get(ref_id, ref_id)


def _ref_fields_for_output(
    ref_id: str,
    ref_norm_to_orig: Optional[Dict[str, str]],
) -> tuple[str, str, str]:
    ref_id_out = _denormalize_ref_id(ref_id, ref_norm_to_orig)
    chrom_id_out, sub_out = split_chrom_subgenome(ref_id_out)
    return ref_id_out, chrom_id_out, sub_out


def _miniprot_ref_from_segment_row(chrom_id: str, sub: str) -> str:
    """Reconstruct reference ID from chromosome and subgenome components.

    Args:
        chrom_id: Chromosome identifier (e.g., "chr1", "Chr5")
        sub: Subgenome identifier (e.g., "A", "B", "At", "Dt") or "NA" if none

    Returns:
        Full reference ID. If sub is "NA", returns just chrom_id.
        Otherwise returns chrom_id + sub (e.g., "chr1A", "chr5At").

    This supports both single-letter subgenomes (Arabidopsis: A/B) and
    multi-character subgenomes (cotton: At/Dt).
    """
    if not chrom_id:
        return ""
    sub = str(sub)
    if sub == "NA":
        return chrom_id
    return f"{chrom_id}{sub}"


def write_contig_summary_tsv(
    output_path: Path,
    classifications: List[ContigClassification],
    contig_orientations: Dict[str, bool],
    all_contig_lengths: Dict[str, int],
    contig_total: Dict[str, int],
    contig_refs: Dict[str, Set[str]],
    qlens_from_paf: Dict[str, int],
    qr_union_bp: Dict[Tuple[str, str], int],
    qr_matches: Dict[Tuple[str, str], int],
    qr_alnlen: Dict[Tuple[str, str], int],
    qr_gene_count: Dict[Tuple[str, str], int],
    best_ref: Dict[str, str],
    best_score: Dict[str, float],
    best_bp: Dict[str, int],
    second_ref: Dict[str, str],
    second_score: Dict[str, float],
    second_bp: Dict[str, int],
    chr_like_minlen: int,
    chimera_primary_frac: float,
    chimera_secondary_frac: float,
    assign_min_frac: float,
    assign_min_ratio: float,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
    scaffold_confidences: Optional[Dict[str, Tuple[float, float, float]]] = None,
) -> None:
    """Write contig_summary.tsv with classification columns.

    Includes classification columns:
    - original_name
    - classification (chrom_assigned, chrom_unassigned, contaminant, etc.)
    - classification_confidence (high/medium/low)
    - reversed
    - contaminant_taxid (NCBI taxonomy ID for contaminants)
    - contaminant_sci (scientific name for contaminants)
    - length (after status)
    - ref_gene_proportion
    - genes_per_Mbp (for chromosome-assigned contigs)

    Evidence strength columns (for classification confidence):
    - gc_content: GC content of contig (0.0-1.0)
    - gc_deviation: Deviation from reference mean in std devs
    - synteny_score: Synteny evidence strength (0.0-1.0)
    - contam_score: Contaminant evidence strength (0.0-1.0)
    - contam_coverage: Contaminant alignment coverage (0.0-1.0)
    """
    # Build lookup by original name
    clf_lookup = {clf.original_name: clf for clf in classifications}

    header = [
        "contig",
        "original_name",
        "classification",
        "classification_confidence",
        "reversed",
        "contaminant_taxid",
        "contaminant_sci",
        "assigned_subgenome",
        "assigned_ref_id",
        "assigned_chrom_id",
        "status",
        "length",
        # Full-length vs fragment fields
        "ref_coverage",
        "is_full_length",
        "full_length_confidence",
        "query_subgenome",
        "seq_identity_vs_ref",
        "has_5p_telomere",
        "has_3p_telomere",
        "rearrangement_candidates",
        # Evidence fields
        "gc_content",
        "gc_deviation",
        "synteny_score",
        "collinearity_score",
        "contam_score",
        "contam_coverage",
        "depth_mean",
        "depth_median",
        "depth_std",
        "depth_breadth_1x",
        "depth_breadth_10x",
        "best_score",
        "second_score",
        "score_ratio",
        "best_ref_union_bp",
        "best_ref_union_frac",
        "ref_gene_proportion",
        "genes_per_Mbp",
        "n_ref_hits",
        "total_aligned_bp",
        "chimeric",
        "chimera_reason",
        "low_coverage",
        "best_matches",
        "best_aln_len",
        "best_identity",
        "best_distance",
        "second_ref_id",
        "second_ref_union_bp",
        "second_matches",
        "second_aln_len",
        "second_identity",
        "second_distance",
        "scaffold_grouping_confidence",
        "scaffold_location_confidence",
        "scaffold_orientation_confidence",
    ]

    def fetch_metrics(q: str, ref_id: str):
        if not ref_id:
            return 0, 0, 0, None, None
        key = (q, ref_id)
        ubp = int(qr_union_bp.get(key, 0) or 0)
        m = int(qr_matches.get(key, 0) or 0)
        al = int(qr_alnlen.get(key, 0) or 0)
        if al > 0:
            ident = m / al
            dist = 1.0 - ident
        else:
            ident, dist = None, None
        return ubp, m, al, ident, dist

    with output_path.open("w") as out:
        out.write("\t".join(header) + "\n")

        for q in sorted(all_contig_lengths.keys()):
            clf = clf_lookup.get(q)

            contig_len = int(all_contig_lengths.get(q, 0) or 0)
            if contig_len <= 0:
                contig_len = int(qlens_from_paf.get(q, 0) or 0)

            total_aligned_bp = int(contig_total.get(q, 0) or 0)
            ref_hits = contig_refs.get(q, set()) or set()
            n_ref_hits = len(ref_hits)

            low_coverage = "yes" if total_aligned_bp < int(chr_like_minlen) else "no"

            # Get chain evidence for reference assignment
            chain_assigned_ref_id = str(best_ref.get(q, "") or "")
            bs = float(best_score.get(q, 0.0) or 0.0)
            sr = float(second_score.get(q, 0.0) or 0.0)

            # For non-chromosome classifications (organelles, rDNA, contaminants, debris),
            # use the classification's assigned_ref_id instead of chain evidence
            clf_class = clf.classification if clf else "unclassified"
            if clf and clf_class in ("organelle_complete", "organelle_debris", "rDNA",
                                     "contaminant", "debris", "chrom_debris", "unclassified"):
                assigned_ref_id = clf.assigned_ref_id if clf.assigned_ref_id else ""
            else:
                assigned_ref_id = chain_assigned_ref_id

            best_union_bp = int(best_bp.get(q, 0) or 0)
            second_ref_id = str(second_ref.get(q, "") or "")
            second_union_bp = int(second_bp.get(q, 0) or 0)

            _mtmp, m_b, al_b, ident_b, dist_b = fetch_metrics(q, assigned_ref_id)
            _stmp, m_s, al_s, ident_s, dist_s = fetch_metrics(q, second_ref_id)

            # For non-chromosome classifications (organelles, debris, etc.), don't require
            # chain evidence (bs > 0) - they have assigned_ref_id from their detection method
            needs_chain_evidence = clf_class not in ("organelle_complete", "organelle_debris",
                                                      "rDNA", "contaminant", "debris",
                                                      "chrom_debris", "unclassified")
            if not assigned_ref_id or (needs_chain_evidence and bs <= 0.0):
                assigned_chrom_id, assigned_subgenome = ("", "NA")
                status = "NO_HITS"
                best_ref_union_frac = 0.0
                score_ratio = 0.0
                assigned_ref_id_out = "NA"
            else:
                assigned_ref_id_out = _denormalize_ref_id(assigned_ref_id, ref_norm_to_orig)
                assigned_chrom_id, assigned_subgenome = split_chrom_subgenome(assigned_ref_id_out)
                best_ref_union_frac = (best_union_bp / contig_len) if contig_len > 0 else 0.0
                score_ratio = (bs / sr) if sr > 0 else (MAX_SCORE_RATIO if bs > 0 else 0.0)

                status = "OK"
                if best_ref_union_frac < float(assign_min_frac):
                    status = "AMBIG_LOW_FRAC"
                elif score_ratio < float(assign_min_ratio):
                    status = "AMBIG_LOW_RATIO"

            chimeric = "no"
            chimera_reason = ""
            if total_aligned_bp > 0 and n_ref_hits > 1:
                primary_frac = (best_union_bp / total_aligned_bp) if total_aligned_bp > 0 else 0.0
                secondary_frac = (second_union_bp / total_aligned_bp) if total_aligned_bp > 0 else 0.0
                if (primary_frac < float(chimera_primary_frac)) and (secondary_frac >= float(chimera_secondary_frac)):
                    chimeric = "yes"
                    chimera_reason = f"primary_frac={primary_frac:.3f},secondary_frac={secondary_frac:.3f}"

            # Classification columns
            original_name = q
            new_name = clf.new_name if clf else q
            classification = clf.classification if clf else "unclassified"
            classification_confidence = clf.classification_confidence if clf and clf.classification_confidence else ""
            reversed_val = "yes" if contig_orientations.get(q, False) else "no"
            contaminant_taxid = str(clf.contaminant_taxid) if clf and clf.contaminant_taxid is not None else ""
            contaminant_sci = clf.contaminant_sci if clf and clf.contaminant_sci else ""
            ref_gene_proportion = clf.ref_gene_proportion if clf and clf.ref_gene_proportion is not None else ""

            # Full-length vs fragment columns
            ref_cov_str = f"{clf.ref_coverage:.4f}" if clf and clf.ref_coverage is not None else ""
            is_full_length = "yes" if clf and clf.is_full_length else ("no" if clf and clf.is_full_length is not None else "")
            fl_confidence = clf.full_length_confidence if clf and clf.full_length_confidence else ""
            query_subgenome = clf.query_subgenome if clf and clf.query_subgenome else ""
            seq_identity = f"{clf.seq_identity_vs_ref:.6f}" if clf and clf.seq_identity_vs_ref is not None else ""
            has_5p_telo = "yes" if clf and clf.has_5p_telomere else ("no" if clf and clf.has_5p_telomere is not None else "")
            has_3p_telo = "yes" if clf and clf.has_3p_telomere else ("no" if clf and clf.has_3p_telomere is not None else "")
            rearrangement_candidates = clf.rearrangement_candidates if clf and clf.rearrangement_candidates else ""

            # Evidence strength columns
            gc_content = f"{clf.gc_content:.4f}" if clf and clf.gc_content is not None else ""
            gc_deviation = f"{clf.gc_deviation:.2f}" if clf and clf.gc_deviation is not None else ""
            synteny_score = f"{clf.synteny_score:.3f}" if clf and clf.synteny_score is not None else ""
            collinearity_score = f"{clf.collinearity_score:.3f}" if clf and clf.collinearity_score is not None else ""
            contam_score = f"{clf.contam_score:.3f}" if clf and clf.contam_score is not None else ""
            contam_coverage = f"{clf.contam_coverage:.3f}" if clf and clf.contam_coverage is not None else ""

            # Read depth columns
            depth_mean = f"{clf.depth_mean:.2f}" if clf and clf.depth_mean is not None else ""
            depth_median = f"{clf.depth_median:.2f}" if clf and clf.depth_median is not None else ""
            depth_std = f"{clf.depth_std:.2f}" if clf and clf.depth_std is not None else ""
            depth_breadth_1x = f"{clf.depth_breadth_1x:.4f}" if clf and clf.depth_breadth_1x is not None else ""
            depth_breadth_10x = f"{clf.depth_breadth_10x:.4f}" if clf and clf.depth_breadth_10x is not None else ""

            # Compute genes per Mbp for chromosome-assigned contigs
            gene_count = qr_gene_count.get((q, assigned_ref_id), 0) if assigned_ref_id else 0
            if contig_len > 0 and gene_count > 0:
                genes_per_Mbp = gene_count / (contig_len / 1_000_000.0)
            else:
                genes_per_Mbp = None

            out.write(
                "\t".join(
                    [
                        new_name,  # contig (new name)
                        original_name,
                        classification,
                        classification_confidence,
                        reversed_val,
                        contaminant_taxid,
                        contaminant_sci,
                        str(assigned_subgenome),
                        str(assigned_ref_id_out),
                        str(assigned_chrom_id),
                        str(status),
                        str(int(contig_len)),  # length
                        # Full-length vs fragment fields
                        ref_cov_str,
                        is_full_length,
                        fl_confidence,
                        query_subgenome,
                        seq_identity,
                        has_5p_telo,
                        has_3p_telo,
                        rearrangement_candidates,
                        # Evidence fields
                        gc_content,
                        gc_deviation,
                        synteny_score,
                        collinearity_score,
                        contam_score,
                        contam_coverage,
                        depth_mean,
                        depth_median,
                        depth_std,
                        depth_breadth_1x,
                        depth_breadth_10x,
                        f"{bs:.3f}",
                        f"{sr:.3f}",
                        f"{score_ratio:.3f}",
                        str(int(best_union_bp)),
                        f"{best_ref_union_frac:.4f}",
                        (f"{ref_gene_proportion:.4f}" if isinstance(ref_gene_proportion, float) else str(ref_gene_proportion)),
                        (f"{genes_per_Mbp:.2f}" if genes_per_Mbp is not None else ""),
                        str(int(n_ref_hits)),
                        str(int(total_aligned_bp)),
                        str(chimeric),
                        str(chimera_reason),
                        str(low_coverage),
                        str(int(m_b)),
                        str(int(al_b)),
                        (f"{ident_b:.6f}" if ident_b is not None else ""),
                        (f"{dist_b:.6f}" if dist_b is not None else ""),
                        (_denormalize_ref_id(second_ref_id, ref_norm_to_orig) if second_ref_id else "NA"),
                        str(int(second_union_bp)),
                        str(int(m_s)),
                        str(int(al_s)),
                        (f"{ident_s:.6f}" if ident_s is not None else ""),
                        (f"{dist_s:.6f}" if dist_s is not None else ""),
                        # Scaffold confidence columns
                        *(lambda sc: (
                            f"{sc[0]:.3f}", f"{sc[1]:.3f}", f"{sc[2]:.3f}"
                        ) if sc else ("", "", ""))(
                            scaffold_confidences.get(q) if scaffold_confidences else None
                        ),
                    ]
                )
                + "\n"
            )


def write_macro_blocks_tsv(
    out_tsv: Path,
    rows,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> None:
    with out_tsv.open("w") as out:
        out.write(
            "\t".join(
                [
                    "contig",
                    "contig_len",
                    "ref_id",
                    "chrom_id",
                    "subgenome",
                    "strand",
                    "chain_id",
                    "qstart",
                    "qend",
                    "qspan_bp",
                    "union_bp",
                    "matches",
                    "aln_len",
                    "identity",
                    "score",
                    "n_segments",
                    "gene_count_chain",
                    "ref_start",
                    "ref_end",
                ]
            )
            + "\n"
        )
        for r in rows:
            # Support both old (17-field) and new (19-field) tuples
            if len(r) >= 19:
                (
                    q, qlen, ref_id, chrom_id, sub, strand, chain_id,
                    qstart, qend, qspan, qbp, msum, alnsum, ident, score,
                    nseg, gene_count_chain, ref_start, ref_end,
                ) = r[:19]
            else:
                (
                    q, qlen, ref_id, chrom_id, sub, strand, chain_id,
                    qstart, qend, qspan, qbp, msum, alnsum, ident, score,
                    nseg, gene_count_chain,
                ) = r[:17]
                ref_start = ""
                ref_end = ""
            ref_id_out, chrom_id_out, sub_out = _ref_fields_for_output(ref_id, ref_norm_to_orig)
            out.write(
                "\t".join(
                    map(
                        str,
                        [
                            q,
                            qlen,
                            ref_id_out,
                            chrom_id_out,
                            sub_out,
                            strand,
                            chain_id,
                            qstart,
                            qend,
                            qspan,
                            qbp,
                            msum,
                            alnsum,
                            ident,
                            score,
                            nseg,
                            gene_count_chain,
                            ref_start,
                            ref_end,
                        ],
                    )
                )
                + "\n"
            )


def write_chain_segments_tsv(
    out_tsv: Path,
    rows,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> None:
    with out_tsv.open("w") as out:
        out.write(
            "\t".join(
                [
                    "contig",
                    "contig_len",
                    "target_chrom_id",
                    "target_subgenome",
                    "strand",
                    "chain_id",
                    "qstart",
                    "qend",
                ]
            )
            + "\n"
        )
        for (q, qlen, chrom_id, sub, strand, chain_id, s, e) in rows:
            ref_id = _miniprot_ref_from_segment_row(str(chrom_id), str(sub))
            _ref_out, chrom_out, sub_out = _ref_fields_for_output(ref_id, ref_norm_to_orig)
            out.write(f"{q}\t{qlen}\t{chrom_out}\t{sub_out}\t{strand}\t{chain_id}\t{s}\t{e}\n")


def write_chain_summary_tsv(
    out_tsv: Path,
    rows,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> None:
    with out_tsv.open("w") as out:
        out.write(
            "\t".join(
                [
                    "contig",
                    "ref_id",
                    "chrom_id",
                    "subgenome",
                    "union_bp",
                    "matches",
                    "aln_len",
                    "identity",
                    "gene_count",
                    "n_chains_kept",
                    "score_topk",
                    "score_all",
                    "best_chain_weight",
                    "best_chain_qbp",
                    "best_chain_identity",
                    "collinearity",
                ]
            )
            + "\n"
        )
        for row in sorted(rows):
            # Support both old (15-field) and new (16-field) tuples
            q, ref_id, chrom_id, sub, ubp, msum, alnsum, ident, gcount, nch, w_topk, w_all, bw, bqbp, bident = row[:15]
            collin = row[15] if len(row) > 15 else None

            ref_id_out, chrom_id_out, sub_out = _ref_fields_for_output(ref_id, ref_norm_to_orig)
            collin_str = f"{float(collin):.3f}" if collin is not None else ""
            out.write(
                "\t".join(
                    [
                        q,
                        ref_id_out,
                        chrom_id_out,
                        sub_out,
                        str(int(ubp)),
                        str(int(msum)),
                        str(int(alnsum)),
                        f"{float(ident):.6f}",
                        str(int(gcount)),
                        str(int(nch)),
                        f"{float(w_topk):.3f}",
                        f"{float(w_all):.3f}",
                        f"{float(bw):.3f}",
                        str(int(bqbp)),
                        f"{float(bident):.6f}",
                        collin_str,
                    ]
                )
                + "\n"
            )


def compute_summary(
    summary_path: Path,
    all_contig_lengths: dict,
    contig_total,
    contig_refs,
    qlens_from_paf,
    qr_union_bp,
    qr_matches,
    qr_alnlen,
    best_ref,
    best_score,
    best_bp,
    second_ref,
    second_score,
    second_bp,
    chr_like_minlen: int,
    chimera_primary_frac: float,
    chimera_secondary_frac: float,
    assign_min_frac: float,
    assign_min_ratio: float,
) -> None:
    """
    Write the canonical per-contig summary TSV.

    Columns (in order):
      contig, assigned_subgenome, assigned_ref_id, assigned_chrom_id, status,
      best_score, second_score, score_ratio, best_ref_union_bp, contig_len, best_ref_union_frac,
      n_ref_hits, total_aligned_bp, chimeric, chimera_reason, low_coverage,
      best_matches, best_aln_len, best_identity, best_distance,
      second_ref_id, second_ref_union_bp, second_matches, second_aln_len, second_identity, second_distance
    """

    header = [
        "contig",
        "assigned_subgenome",
        "assigned_ref_id",
        "assigned_chrom_id",
        "status",
        "best_score",
        "second_score",
        "score_ratio",
        "best_ref_union_bp",
        "contig_len",
        "best_ref_union_frac",
        "n_ref_hits",
        "total_aligned_bp",
        "chimeric",
        "chimera_reason",
        "low_coverage",
        "best_matches",
        "best_aln_len",
        "best_identity",
        "best_distance",
        "second_ref_id",
        "second_ref_union_bp",
        "second_matches",
        "second_aln_len",
        "second_identity",
        "second_distance",
    ]

    def fetch_metrics(q: str, ref_id: str):
        """
        For a given (contig, ref_id), return:
          union_bp, matches, aln_len, identity, distance
        where identity = matches/aln_len, distance = 1-identity, or blanks if aln_len=0.
        """
        if not ref_id:
            return 0, 0, 0, None, None

        key = (q, ref_id)
        ubp = int(qr_union_bp.get(key, 0) or 0)
        m = int(qr_matches.get(key, 0) or 0)
        al = int(qr_alnlen.get(key, 0) or 0)

        if al > 0:
            ident = m / al
            dist = 1.0 - ident
        else:
            ident, dist = None, None

        return ubp, m, al, ident, dist

    with summary_path.open("w") as out:
        out.write("\t".join(header) + "\n")

        for q in sorted(all_contig_lengths.keys()):
            contig_len = int(all_contig_lengths.get(q, 0) or 0)
            if contig_len <= 0:
                contig_len = int(qlens_from_paf.get(q, 0) or 0)

            total_aligned_bp = int(contig_total.get(q, 0) or 0)
            ref_hits = contig_refs.get(q, set()) or set()
            n_ref_hits = len(ref_hits)

            low_coverage = "yes" if total_aligned_bp < int(chr_like_minlen) else "no"

            # Chosen assignment (after any gate-aware reranking)
            assigned_ref_id = str(best_ref.get(q, "") or "")
            bs = float(best_score.get(q, 0.0) or 0.0)
            sr = float(second_score.get(q, 0.0) or 0.0)

            best_union_bp = int(best_bp.get(q, 0) or 0)
            second_ref_id = str(second_ref.get(q, "") or "")
            second_union_bp = int(second_bp.get(q, 0) or 0)

            _mtmp, m_b, al_b, ident_b, dist_b = fetch_metrics(q, assigned_ref_id)
            _stmp, m_s, al_s, ident_s, dist_s = fetch_metrics(q, second_ref_id)

            # Map-style QC
            if not assigned_ref_id or bs <= 0.0:
                assigned_chrom_id, assigned_subgenome = ("", "NA")
                status = "NO_HITS"
                best_ref_union_frac = 0.0
                score_ratio = 0.0
                assigned_ref_id_out = "NA"
            else:
                assigned_chrom_id, assigned_subgenome = split_chrom_subgenome(assigned_ref_id)
                assigned_ref_id_out = assigned_ref_id

                best_ref_union_frac = (best_union_bp / contig_len) if contig_len > 0 else 0.0
                score_ratio = (bs / sr) if sr > 0 else (MAX_SCORE_RATIO if bs > 0 else 0.0)

                status = "OK"
                if best_ref_union_frac < float(assign_min_frac):
                    status = "AMBIG_LOW_FRAC"
                elif score_ratio < float(assign_min_ratio):
                    status = "AMBIG_LOW_RATIO"

            # Chimeric flag based on union_bp fractions across refs (same logic as before)
            chimeric = "no"
            chimera_reason = ""
            if total_aligned_bp > 0 and n_ref_hits > 1:
                primary_frac = (best_union_bp / total_aligned_bp) if total_aligned_bp > 0 else 0.0
                secondary_frac = (second_union_bp / total_aligned_bp) if total_aligned_bp > 0 else 0.0
                if (primary_frac < float(chimera_primary_frac)) and (secondary_frac >= float(chimera_secondary_frac)):
                    chimeric = "yes"
                    chimera_reason = f"primary_frac={primary_frac:.3f},secondary_frac={secondary_frac:.3f}"

            out.write(
                "\t".join(
                    [
                        q,
                        str(assigned_subgenome),
                        str(assigned_ref_id_out),
                        str(assigned_chrom_id),
                        str(status),
                        f"{bs:.3f}",
                        f"{sr:.3f}",
                        f"{score_ratio:.3f}",
                        str(int(best_union_bp)),
                        str(int(contig_len)),
                        f"{best_ref_union_frac:.4f}",
                        str(int(n_ref_hits)),
                        str(int(total_aligned_bp)),
                        str(chimeric),
                        str(chimera_reason),
                        str(low_coverage),
                        str(int(m_b)),
                        str(int(al_b)),
                        (f"{ident_b:.6f}" if ident_b is not None else ""),
                        (f"{dist_b:.6f}" if dist_b is not None else ""),
                        (second_ref_id if second_ref_id else "NA"),
                        str(int(second_union_bp)),
                        str(int(m_s)),
                        str(int(al_s)),
                        (f"{ident_s:.6f}" if ident_s is not None else ""),
                        (f"{dist_s:.6f}" if dist_s is not None else ""),
                    ]
                )
                + "\n"
            )


def build_segment_support_from_rows(segment_rows):
    """
    segment_rows tuples:
      (contig, contig_len, chrom_id, sub, strand, chain_id, s, e)

    Returns:
      seg_count[(q, ref_id)] = number of segments
      span_bp[(q, ref_id)]   = max(qend) - min(qstart)
    """
    seg_count = defaultdict(int)
    min_s: dict[tuple[str, str], int] = {}
    max_e: dict[tuple[str, str], int] = {}

    for (q, _qlen, chrom_id, sub, _strand, _chain_id, s, e) in segment_rows:
        ref_id = _miniprot_ref_from_segment_row(str(chrom_id), str(sub))
        if not ref_id:
            continue

        key = (q, ref_id)
        seg_count[key] += 1

        s = int(s)
        e = int(e)
        if key not in min_s or s < min_s[key]:
            min_s[key] = s
        if key not in max_e or e > max_e[key]:
            max_e[key] = e

    span_bp = {k: max(0, int(max_e[k]) - int(min_s[k])) for k in min_s.keys()}
    return seg_count, span_bp


def write_contaminant_summary_tsv(
    output_path: Path,
    contaminants: Dict[str, ContaminantHitExtended],
    query_lengths: Dict[str, int],
    depth_stats: Optional[Dict[str, DepthStats]] = None,
) -> None:
    """Write detailed contaminant summary TSV with taxonomic lineage.

    This TSV is used for the contamination alluvial plot visualization,
    showing the phylogenetic breakdown of detected contaminants.

    Args:
        output_path: Path to output TSV file
        contaminants: Dict mapping contig name -> ContaminantHitExtended
        query_lengths: Dict mapping contig name -> length
        depth_stats: Optional dict mapping contig name -> DepthStats
    """
    columns = [
        "contig",
        "length",
        "taxid",
        "sci_name",
        "score",
        "coverage",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "depth_mean",
    ]

    with output_path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")

        for contig_name in sorted(contaminants.keys()):
            hit = contaminants[contig_name]
            length = query_lengths.get(contig_name, 0)
            ds = depth_stats.get(contig_name) if depth_stats else None

            row = [
                contig_name,
                str(length),
                str(hit.taxid),
                hit.sci_name or "",
                str(hit.score),
                f"{hit.coverage:.4f}",
                hit.kingdom or "",
                hit.phylum or "",
                hit.tax_class or "",  # 'class' column
                hit.order or "",
                hit.family or "",
                hit.genus or "",
                hit.species or "",
                f"{ds.mean_depth:.2f}" if ds else "",
            ]
            fh.write("\t".join(row) + "\n")


def write_rdna_annotations_tsv(
    output_path: Path,
    loci: List[RdnaLocus],
    classifications: Optional[Dict[str, str]] = None,
) -> None:
    """Write rDNA annotation TSV with per-locus details.

    Args:
        output_path: Path to output TSV file
        loci: List of RdnaLocus annotations
        classifications: Optional dict mapping contig name -> classification string
    """
    columns = [
        "contig",
        "start",
        "end",
        "strand",
        "identity",
        "consensus_coverage",
        "copy_type",
        "sub_features",
        "array_id",
        "contig_classification",
    ]

    with output_path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")

        for locus in sorted(loci, key=lambda l: (l.contig, l.start)):
            clf = classifications.get(locus.contig, "") if classifications else ""
            row = [
                locus.contig,
                str(locus.start),
                str(locus.end),
                locus.strand,
                f"{locus.identity:.4f}",
                f"{locus.consensus_coverage:.4f}",
                locus.copy_type,
                ",".join(locus.sub_features) if locus.sub_features else "",
                locus.array_id if locus.array_id else "",
                clf,
            ]
            fh.write("\t".join(row) + "\n")


def write_rdna_arrays_tsv(
    output_path: Path,
    arrays: List[RdnaArray],
    classifications: Optional[Dict[str, str]] = None,
) -> None:
    """Write rDNA array summary TSV.

    Args:
        output_path: Path to output TSV file
        arrays: List of RdnaArray objects
        classifications: Optional dict mapping contig name -> classification string
    """
    columns = [
        "array_id",
        "contig",
        "start",
        "end",
        "span_kb",
        "strand",
        "n_total",
        "n_full",
        "n_partial",
        "n_fragment",
        "identity_median",
        "identity_min",
        "identity_max",
        "contig_classification",
    ]

    with output_path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")

        for arr in arrays:
            clf = classifications.get(arr.contig, "") if classifications else ""
            row = [
                arr.array_id,
                arr.contig,
                str(arr.start),
                str(arr.end),
                f"{arr.span / 1000:.1f}",
                arr.strand,
                str(arr.n_total),
                str(arr.n_full),
                str(arr.n_partial),
                str(arr.n_fragment),
                f"{arr.identity_median:.4f}",
                f"{arr.identity_min:.4f}",
                f"{arr.identity_max:.4f}",
                clf,
            ]
            fh.write("\t".join(row) + "\n")


# Sequence Ontology term mappings for rDNA sub-features
_RDNA_SO_TERMS = {
    "18S": "rRNA_18S",  # SO:0001000
    "5.8S": "rRNA_5_8S",  # SO:0001001
    "25S": "rRNA_28S",  # SO:0001002 (25S is the plant homolog of 28S)
    "28S": "rRNA_28S",  # SO:0001002
    "ITS1": "internal_transcribed_spacer_region",  # SO:0000372
    "ITS2": "internal_transcribed_spacer_region",  # SO:0000372
}


def _escape_gff3_value(value: str) -> str:
    """Escape special characters in GFF3 attribute values.

    GFF3 requires URL-encoding of: tab, newline, semicolon, equals, percent, ampersand.
    """
    # Characters that need to be percent-encoded in GFF3
    replacements = [
        ("%", "%25"),  # Must be first
        ("\t", "%09"),
        ("\n", "%0A"),
        ("\r", "%0D"),
        (";", "%3B"),
        ("=", "%3D"),
        ("&", "%26"),
        (",", "%2C"),
    ]
    for char, encoded in replacements:
        value = value.replace(char, encoded)
    return value


def write_rdna_annotations_gff3(
    output_path: Path,
    loci: List[RdnaLocus],
    query_lengths: Dict[str, int],
    classifications: Optional[Dict[str, str]] = None,
    arrays: Optional[List[RdnaArray]] = None,
) -> None:
    """Write rDNA annotation GFF3 with hierarchical feature structure.

    Produces a properly formatted GFF3 file with:
    - repeat_region features (SO:0000657) as parent features for rDNA arrays
    - rRNA_gene features (SO:0001637) as parent features for each locus
    - rRNA sub-features (18S, 5.8S, 28S) and ITS regions as children

    GFF3 uses 1-based, fully-closed coordinates [start, end].

    Args:
        output_path: Path to output GFF3 file
        loci: List of RdnaLocus annotations
        query_lengths: Dict mapping contig name -> length (for ##sequence-region)
        classifications: Optional dict mapping contig name -> classification string
        arrays: Optional list of RdnaArray objects for repeat_region parent features
    """
    # Build array lookup for loci that belong to arrays
    array_lookup: Dict[str, RdnaArray] = {}
    if arrays:
        for arr in arrays:
            array_lookup[arr.array_id] = arr

    # Get unique contigs from loci, sorted
    contigs_with_loci = sorted(set(l.contig for l in loci))

    with output_path.open("w") as fh:
        # GFF3 header
        fh.write("##gff-version 3\n")

        # Write sequence-region pragmas for contigs with annotations
        for contig in contigs_with_loci:
            length = query_lengths.get(contig, 0)
            if length > 0:
                fh.write(f"##sequence-region {contig} 1 {length}\n")

        # Write array-level repeat_region features first
        written_arrays: Set[str] = set()
        if arrays:
            for arr in arrays:
                arr_gff_start = arr.start + 1
                arr_gff_end = arr.end
                arr_feature_id = f"rdna_{arr.array_id}"

                arr_attrs = [
                    f"ID={_escape_gff3_value(arr_feature_id)}",
                    f"Name=45S_rDNA_array",
                    f"repeat_type=rDNA_45S",
                    f"n_copies={arr.n_total}",
                    f"n_full={arr.n_full}",
                    f"identity_median={arr.identity_median:.4f}",
                ]

                fh.write("\t".join([
                    arr.contig,
                    "final_finalizer",
                    "repeat_region",  # SO:0000657
                    str(arr_gff_start),
                    str(arr_gff_end),
                    ".",
                    arr.strand,
                    ".",
                    ";".join(arr_attrs),
                ]) + "\n")
                written_arrays.add(arr.array_id)

        # Write features sorted by contig and position
        for locus in sorted(loci, key=lambda l: (l.contig, l.start)):
            # Convert to 1-based GFF3 coordinates
            gff_start = locus.start + 1  # Convert from 0-based to 1-based
            gff_end = locus.end  # Already exclusive, so this is the last included base

            # Build parent feature ID
            parent_id = f"rRNA_gene_{locus.contig}_{gff_start}_{gff_end}"

            # Build attributes for the parent rRNA_gene feature
            attrs = [f"ID={_escape_gff3_value(parent_id)}"]

            # Add Parent attribute if locus belongs to an array
            if locus.array_id and locus.array_id in written_arrays:
                attrs.append(f"Parent=rdna_{locus.array_id}")

            attrs.append(f"copy_type={locus.copy_type}")
            attrs.append(f"identity={locus.identity:.4f}")

            if locus.array_id:
                attrs.append(f"array_id={locus.array_id}")

            if classifications:
                clf = classifications.get(locus.contig, "")
                if clf:
                    attrs.append(f"contig_classification={_escape_gff3_value(clf)}")

            # Write parent feature (rRNA_gene)
            fh.write("\t".join([
                locus.contig,
                "final_finalizer",
                "rRNA_gene",  # SO:0001637
                str(gff_start),
                str(gff_end),
                ".",  # score
                locus.strand,
                ".",  # phase (not applicable)
                ";".join(attrs),
            ]) + "\n")

            # Write sub-features with mapped coordinates
            for sf in locus.sub_feature_loci:
                # Convert sub-feature coordinates to 1-based
                sf_gff_start = sf.start + 1
                sf_gff_end = sf.end  # Already exclusive

                # Get SO term for this feature type
                so_type = _RDNA_SO_TERMS.get(sf.name, "rRNA")

                # Build sub-feature ID
                sf_id = f"{sf.name}_{locus.contig}_{sf_gff_start}_{sf_gff_end}"

                sf_attrs = [
                    f"ID={_escape_gff3_value(sf_id)}",
                    f"Parent={_escape_gff3_value(parent_id)}",
                ]

                # For ITS regions, add Name to distinguish ITS1 from ITS2
                if sf.name in ("ITS1", "ITS2"):
                    sf_attrs.append(f"Name={sf.name}")

                fh.write("\t".join([
                    locus.contig,
                    "final_finalizer",
                    so_type,
                    str(sf_gff_start),
                    str(sf_gff_end),
                    ".",  # score
                    locus.strand,
                    ".",  # phase
                    ";".join(sf_attrs),
                ]) + "\n")


def write_rearrangements_tsv(
    output_path: Path,
    calls: List[RearrangementCall],
) -> None:
    """Write rearrangement calls to a TSV file.

    Args:
        output_path: Output TSV path.
        calls: List of RearrangementCall instances.
    """
    header = [
        "contig",
        "original_name",
        "assigned_ref_id",
        "rearrangement_type",
        "partner_ref_id",
        "query_start",
        "query_end",
        "ref_start",
        "ref_end",
        "span_bp",
        "strand",
        "confidence",
        "evidence",
        "caveat",
    ]
    with output_path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for c in calls:
            fh.write("\t".join(str(v) if v is not None else "" for v in [
                c.contig,
                c.original_name,
                c.assigned_ref_id,
                c.rearrangement_type,
                c.partner_ref_id,
                c.query_start,
                c.query_end,
                c.ref_start,
                c.ref_end,
                c.span_bp,
                c.strand,
                c.confidence,
                c.evidence,
                c.caveat,
            ]) + "\n")
