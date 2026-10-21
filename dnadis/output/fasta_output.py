#!/usr/bin/env python3
"""Output FASTA files for classified contigs."""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from dnadis.models import ContigClassification
from dnadis.utils.io_utils import symlink_if_needed
from dnadis.utils.logging_config import get_logger

logger = get_logger("fasta_output")

from dnadis.utils.reference_utils import split_chrom_subgenome
from dnadis.utils.sequence_utils import (
    read_fasta_sequences,
    reverse_complement,
    write_fasta_if_changed,
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


# Output categories, in write order.  Chromosome contigs (chrom_assigned and
# chrom_fragment) go to "chrs"; with --scaffold the pseudomolecules go to a
# separate "chrs.scaffolds" file.
FASTA_CATEGORIES = ["chrs", "organelles", "rdna", "cobionts", "circles", "debris", "unclassified"]

_CATEGORY_BY_CLASS = {
    "chrom_assigned": "chrs",
    "chrom_fragment": "chrs",
    "chrom_unassigned": "unclassified",
    "organelle_complete": "organelles",
    "rDNA": "rdna",
    "cobiont": "cobionts",
    "circular_element": "circles",
    "chrom_debris": "debris",
    "debris": "debris",
    "organelle_debris": "debris",
}


def fasta_ext(compress: bool) -> str:
    """FASTA output extension: bgzip-compressed unless *compress* is False."""
    return ".fasta.gz" if compress else ".fasta"


def classified_fasta_path(output_prefix: Path, category: str, compress: bool = True) -> Path:
    """Path of a classified FASTA output, e.g. ``{prefix}.chrs.fasta.gz``."""
    return Path(f"{output_prefix}.{category}{fasta_ext(compress)}")


def write_classified_fastas(
    query_fasta: Path,
    classifications: List[ContigClassification],
    contig_orientations: Dict[str, bool],
    output_prefix: Path,
    scaffolded_seqs: Optional[Dict[str, str]] = None,
    compress: bool = True,
    threads: int = 1,
) -> Dict[str, Path]:
    """Write the classified FASTA files (bgzip-compressed by default).

    ``chrs`` always holds the oriented chromosome contigs.  When
    *scaffolded_seqs* is given, the pseudomolecules are written to
    ``{prefix}.chrs.scaffolds.fasta.gz``; if they are identical to the
    contigs (no joins), that file is a symlink to the chrs FASTA.

    Returns dict mapping category -> output path:
    - chrs, organelles, rdna, cobionts, circles, debris, unclassified
    - chrs_scaffolds (only when *scaffolded_seqs* is given)
    """
    sequences = read_fasta_sequences(query_fasta)

    category_contigs: Dict[str, List[ContigClassification]] = defaultdict(list)
    for clf in classifications:
        category = _CATEGORY_BY_CLASS.get(clf.classification, "unclassified")
        category_contigs[category].append(clf)

    output_paths: Dict[str, Path] = {}
    chrs_seqs: Dict[str, str] = {}

    for category in FASTA_CATEGORIES:
        output_path = classified_fasta_path(output_prefix, category, compress)
        output_paths[category] = output_path

        # Build output sequences dict with new names and orientation
        out_seqs: Dict[str, str] = {}
        for clf in category_contigs.get(category, []):
            orig_seq = sequences.get(clf.original_name, "")
            if not orig_seq:
                continue
            if contig_orientations.get(clf.original_name, False):
                orig_seq = reverse_complement(orig_seq)
            out_seqs[clf.new_name] = orig_seq

        sorted_seqs = dict(sorted(out_seqs.items(), key=lambda x: sort_key(x[0])))
        write_fasta_if_changed(sorted_seqs, output_path, threads=threads)
        if sorted_seqs:
            logger.done(f"{category}: {output_path} ({len(sorted_seqs)} contigs)")
        if category == "chrs":
            chrs_seqs = sorted_seqs

    if scaffolded_seqs:
        scaf_path = classified_fasta_path(output_prefix, "chrs.scaffolds", compress)
        output_paths["chrs_scaffolds"] = scaf_path
        sorted_scaf = dict(sorted(scaffolded_seqs.items(), key=lambda x: sort_key(x[0])))
        if list(sorted_scaf.items()) == list(chrs_seqs.items()):
            symlink_if_needed(scaf_path, output_paths["chrs"])
            logger.done(
                f"chrs.scaffolds: {scaf_path} -> {output_paths['chrs'].name} "
                f"(no joins; identical to chrs)"
            )
        else:
            write_fasta_if_changed(sorted_scaf, scaf_path, threads=threads)
            logger.done(f"chrs.scaffolds: {scaf_path} ({len(sorted_scaf)} sequences)")

    return output_paths


def write_per_subgenome_chrs_fastas(
    chrs_fasta: Path,
    output_prefix: Path,
    kind: str = "chrs",
    threads: int = 1,
) -> Dict[str, Path]:
    """Split an oriented chrs FASTA into per-subgenome FASTA files.

    Reads *chrs_fasta* (headers are new_names like chr1A), groups sequences
    by subgenome using split_chrom_subgenome(), and writes one FASTA per
    subgenome: ``{prefix}.{kind}.{sg}.fasta[.gz]`` (compressed like
    *chrs_fasta*).  When every sequence belongs to a single subgenome, that
    file is a symlink to *chrs_fasta* instead of a copy.

    Args:
        chrs_fasta: Path to the oriented chrs (or chrs.scaffolds) FASTA.
        output_prefix: Output prefix for per-subgenome files.
        kind: File-name infix: ``chrs`` or ``chrs.scaffolds``.
        threads: Compression threads.

    Returns:
        Dict mapping subgenome -> Path, e.g. {"A": Path("...chrs.A.fasta.gz")}.
        Empty dict if the FASTA is missing/empty or all subgenomes are "NA".
    """
    if not chrs_fasta.exists():
        return {}

    sequences = read_fasta_sequences(chrs_fasta)
    if not sequences:
        return {}

    sg_seqs: Dict[str, Dict[str, str]] = defaultdict(dict)
    for name, seq in sequences.items():
        _chrom_id, sg = split_chrom_subgenome(name)
        if sg == "NA":
            continue
        sg_seqs[sg][name] = seq

    if not sg_seqs:
        return {}

    ext = ".fasta.gz" if chrs_fasta.name.endswith(".gz") else ".fasta"
    single_set = len(sg_seqs) == 1 and len(next(iter(sg_seqs.values()))) == len(sequences)

    result: Dict[str, Path] = {}
    for sg in sorted(sg_seqs):
        out_path = Path(f"{output_prefix}.{kind}.{sg}{ext}")
        if single_set:
            symlink_if_needed(out_path, chrs_fasta)
            logger.info(f"Per-subgenome FASTA ({sg}): {out_path} -> {chrs_fasta.name}")
        else:
            sorted_seqs = dict(sorted(sg_seqs[sg].items(), key=lambda x: sort_key(x[0])))
            write_fasta_if_changed(sorted_seqs, out_path, threads=threads)
            logger.info(f"Per-subgenome FASTA ({sg}): {out_path} ({len(sorted_seqs)} contigs)")
        result[sg] = out_path

    return result


_OUTPUT_FASTA_RE = (
    r"\.(?:chrs(?:\.[^.]+)*|organelles|rdna|cobionts|circles|circular_elements"
    r"|debris|unclassified|scaffolded)\.fasta(?:\.gz)?"
    r"|\.(?:scaffolded|chrs\.scaffolds)\.agp"
)


def find_stale_outputs(output_prefix: Path, current: Iterable[Path]) -> List[Path]:
    """FASTA/AGP outputs at *output_prefix* that this run did not write.

    Catches files left by earlier runs under other names (uncompressed or
    pre-rename outputs, scaffolds from a run with --scaffold, subgenome
    splits that no longer apply).
    """
    current_names = {p.name for p in current}
    pattern = re.compile(re.escape(output_prefix.name) + "(?:" + _OUTPUT_FASTA_RE + ")$")
    stale = []
    for p in sorted(output_prefix.parent.glob(f"{output_prefix.name}.*")):
        if p.name not in current_names and pattern.match(p.name):
            stale.append(p)
    return stale
