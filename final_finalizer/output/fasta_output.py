#!/usr/bin/env python3
"""Output FASTA files for classified contigs."""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from final_finalizer.models import ContigClassification
from final_finalizer.utils.logging_config import get_logger

logger = get_logger("fasta_output")

from final_finalizer.utils.sequence_utils import (
    read_fasta_sequences,
    reverse_complement,
    write_fasta,
)


def sort_key(name: str):
    """Sort key for chromosome/contig names (numeric-aware)."""
    m = re.match(r"[cC]hr(\d+)", name)
    if m:
        return (0, int(m.group(1)), name)
    m = re.match(r"contig_(\d+)", name)
    if m:
        return (1, int(m.group(1)), name)
    return (2, 0, name)


def write_classified_fastas(
    query_fasta: Path,
    classifications: List[ContigClassification],
    contig_orientations: Dict[str, bool],
    output_prefix: Path,
    scaffolded_seqs: Optional[Dict[str, str]] = None,
) -> Dict[str, Path]:
    """Write 6 classified FASTA files.

    When scaffolded_seqs is provided, the chrs.fasta output contains
    scaffolded pseudomolecules instead of individual oriented contigs.

    Returns dict mapping category -> output path:
    - {prefix}.chrs.fasta
    - {prefix}.organelles.fasta
    - {prefix}.rdna.fasta
    - {prefix}.contaminants.fasta
    - {prefix}.debris.fasta
    - {prefix}.unclassified.fasta
    """
    # Read all sequences
    sequences = read_fasta_sequences(query_fasta)

    # Group contigs by classification category
    category_contigs: Dict[str, List[ContigClassification]] = defaultdict(list)

    for clf in classifications:
        # Map classification to output category
        if clf.classification == "chrom_assigned":
            category = "chrs"
        elif clf.classification == "chrom_unassigned":
            category = "unclassified"
        elif clf.classification == "organelle_complete":
            category = "organelles"
        elif clf.classification == "rDNA":
            category = "rdna"
        elif clf.classification == "contaminant":
            category = "contaminants"
        elif clf.classification in ("chrom_debris", "debris", "organelle_debris"):
            category = "debris"
        else:
            category = "unclassified"

        category_contigs[category].append(clf)

    # Write each category
    output_paths: Dict[str, Path] = {}
    categories = ["chrs", "organelles", "rdna", "contaminants", "debris", "unclassified"]

    for category in categories:
        output_path = Path(f"{output_prefix}.{category}.fasta")
        output_paths[category] = output_path

        clfs = category_contigs.get(category, [])

        # For chrs category with scaffolded sequences, write scaffolds instead
        if category == "chrs" and scaffolded_seqs:
            if scaffolded_seqs:
                sorted_seqs = dict(sorted(
                    scaffolded_seqs.items(),
                    key=lambda x: sort_key(x[0]),
                ))
                write_fasta(sorted_seqs, output_path)
                logger.done(f"{category}: {output_path} ({len(sorted_seqs)} scaffolded chromosomes)")
            else:
                output_path.write_text("")
            continue

        if not clfs:
            # Write empty file
            output_path.write_text("")
            continue

        # Build output sequences dict with new names and orientation
        out_seqs: Dict[str, str] = {}
        for clf in clfs:
            orig_seq = sequences.get(clf.original_name, "")
            if not orig_seq:
                continue

            # Reverse complement if needed
            if contig_orientations.get(clf.original_name, False):
                orig_seq = reverse_complement(orig_seq)

            out_seqs[clf.new_name] = orig_seq

        sorted_seqs = dict(sorted(out_seqs.items(), key=lambda x: sort_key(x[0])))
        write_fasta(sorted_seqs, output_path)

        logger.done(f"{category}: {output_path} ({len(sorted_seqs)} contigs)")

    return output_paths
