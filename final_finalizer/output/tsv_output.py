#!/usr/bin/env python3
"""TSV output functions for contig summaries and chain data."""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.models import ContigClassification, ContaminantHitExtended, DepthStats
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
) -> None:
    """Write enhanced contig_summary.tsv with classification columns.

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
        # Evidence fields
        "gc_content",
        "gc_deviation",
        "synteny_score",
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

            # Evidence strength columns
            gc_content = f"{clf.gc_content:.4f}" if clf and clf.gc_content is not None else ""
            gc_deviation = f"{clf.gc_deviation:.2f}" if clf and clf.gc_deviation is not None else ""
            synteny_score = f"{clf.synteny_score:.3f}" if clf and clf.synteny_score is not None else ""
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
                        # Evidence fields
                        gc_content,
                        gc_deviation,
                        synteny_score,
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
                ]
            )
            + "\n"
        )
        for r in rows:
            (
                q,
                qlen,
                ref_id,
                chrom_id,
                sub,
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
            ) = r
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
                ]
            )
            + "\n"
        )
        for (
            q,
            ref_id,
            chrom_id,
            sub,
            ubp,
            msum,
            alnsum,
            ident,
            gcount,
            nch,
            w_topk,
            w_all,
            bw,
            bqbp,
            bident,
        ) in sorted(rows):
            ref_id_out, chrom_id_out, sub_out = _ref_fields_for_output(ref_id, ref_norm_to_orig)
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
