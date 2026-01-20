#!/usr/bin/env python3
"""Legacy primary-only PAF stats functions for QA diagnostics."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Tuple

from final_finalizer.utils.io_utils import merge_intervals, open_maybe_gzip
from final_finalizer.utils.reference_utils import is_primary_only


def parse_paf_primary(paf_gz_path: Path, aln_minlen: int):
    """
    Parse PAF (gz or plain) and aggregate stats per (contig, ref_id), PRIMARY ONLY.
    """
    aln_intervals = defaultdict(list)  # (q, ref_id) -> list[(qs,qe)]
    match_sum = defaultdict(int)       # (q, ref_id) -> matches sum
    alnlen_sum = defaultdict(int)      # (q, ref_id) -> aln_len sum
    contig_refs = defaultdict(set)     # q -> set(ref_ids)
    contig_qlen = {}                   # q -> qlen from PAF (if available)

    with open_maybe_gzip(paf_gz_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            if not is_primary_only(fields):
                continue

            qname = fields[0]
            try:
                qlen = int(fields[1])
                qs = int(fields[2])
                qe = int(fields[3])
            except ValueError:
                continue

            ref_id = fields[5]
            try:
                matches = int(fields[9])
                aln_len = int(fields[10])
            except ValueError:
                continue

            if aln_len < aln_minlen:
                continue

            if qe < qs:
                qs, qe = qe, qs
            if qe <= qs:
                continue

            contig_qlen.setdefault(qname, qlen)

            key = (qname, ref_id)
            aln_intervals[key].append((qs, qe))
            match_sum[key] += matches
            alnlen_sum[key] += aln_len
            contig_refs[qname].add(ref_id)

    aln_sum = defaultdict(int)
    for key, ivs in aln_intervals.items():
        _, total = merge_intervals(ivs)
        aln_sum[key] = total

    contig_total = defaultdict(int)
    for (q, _ref_id), abp in aln_sum.items():
        contig_total[q] += abp

    best_ref_id = defaultdict(str)
    best_aln = defaultdict(int)
    second_best_ref_id = defaultdict(str)
    second_best_aln = defaultdict(int)

    for (q, ref_id), abp in aln_sum.items():
        if abp > best_aln[q]:
            second_best_aln[q] = best_aln[q]
            second_best_ref_id[q] = best_ref_id[q]
            best_aln[q] = abp
            best_ref_id[q] = ref_id
        elif abp > second_best_aln[q]:
            second_best_aln[q] = abp
            second_best_ref_id[q] = ref_id

    return (
        aln_sum,
        match_sum,
        alnlen_sum,
        contig_total,
        contig_refs,
        contig_qlen,
        best_aln,
        best_ref_id,
        second_best_aln,
        second_best_ref_id,
    )


def write_stats_tsv(
    stats_path: Path,
    aln_sum: Dict[Tuple[str, str], int],
    match_sum: Dict[Tuple[str, str], int],
) -> None:
    with stats_path.open("w") as out:
        out.write("contig\tref_id\taligned_bp\tmatches\n")
        for (q, ref_id), aln in sorted(aln_sum.items()):
            m = match_sum[(q, ref_id)]
            out.write(f"{q}\t{ref_id}\t{aln}\t{m}\n")
