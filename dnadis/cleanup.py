#!/usr/bin/env python3
"""``dnadis cleanup RESULTS_DIR``: remove bulky intermediates after a run.

Removes files that no later step reads and that a rerun recreates cheaply
(copies of the input split by class, BLAST databases, RagTag work
directories, compleasm's raw miniprot/hmmsearch output).  Keeps everything
the reports, comparison, phylogeny and rerun caches rely on: compleasm
summaries, full tables and translated proteins; read-depth regions and their
cache metadata; synteny PAFs; small cached alignment and centrifuger results.

Read-alignment BAMs (kept only with --keep-depth-bam) are removed unless
--keep-bams is given.

Opt-in:
  --sister-anchors   sister-anchoring anchor/candidate FASTAs (large; a rerun
                     with --sister-anchoring redoes the sister alignments)
  --legacy-outputs   FASTA/AGP outputs under pre-bgzip names, only where the
                     current-name replacement exists
"""
from __future__ import annotations

import argparse
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from dnadis.utils.logging_config import get_logger, setup_logging

logger = get_logger("cleanup")

# Inside each compleasm lineage directory: everything except the files
# dnadis reads (summary.txt one level up, full_table_busco_format.tsv,
# translated_protein.fasta).  The .done markers go with the outputs they
# vouch for, so a forced compleasm rerun does not skip those steps.
_COMPLEASM_BULK = [
    "miniprot_output.gff", "hmmer_output", "detected_genes.gff",
    "gene_marker.fasta", "full_table.tsv", "miniprot.done", "hmmsearch.done",
]

# Per-assembly work directory ({asm}/{asm}_classification/), as glob patterns.
_WORKDIR_BULK = [
    "non_chromosome_contigs.fa",
    "residual_for_cobiont_screen.fa",
    "organelles",
    "rdna",
    "rdna_period",
    "scaffolding",
    "chr_debris/*.fa",
    "redundancy/*.fa",
    "debris/*.fa",
    "cobionts/contigs.fq",
    "rdna_consensus/reannotation/query_db.*",
    "read_depth/*.mosdepth.*.dist.txt",
    "read_depth/*.mosdepth.summary.txt",
    "read_depth/downsampled_reads*",
]

# Read alignments behind the depth track (present with --keep-depth-bam).
_READ_BAMS = ["read_depth/*.bam", "read_depth/*.bam.bai", "read_depth/*.bam.csi"]

_LEGACY_RENAMES = {
    # old suffix -> current suffix that must exist before the old one goes
    ".scaffolded.fasta": ".chrs.scaffolds.agp",
    ".scaffolded.agp": ".chrs.scaffolds.agp",
    ".circular_elements.fasta": ".circles.fasta.gz",
}


@dataclass
class CleanupPlan:
    """Paths to remove, grouped by a short label for the summary."""

    groups: Dict[str, List[Path]] = field(default_factory=dict)

    def add(self, label: str, paths: Iterable[Path]) -> None:
        for p in paths:
            if p.is_symlink() or p.exists():
                self.groups.setdefault(label, []).append(p)


def _size(path: Path) -> int:
    """Bytes on disk under *path*, not following symlinks."""
    if path.is_symlink() or path.is_file():
        return path.lstat().st_size
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(root, f)).st_size
            except OSError:
                pass
    return total


def _remove(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TB"


def _assembly_dirs(results_dir: Path) -> List[Path]:
    """Assembly output directories: those holding an ``{name}_classification`` dir."""
    candidates = [results_dir] + sorted(p for p in results_dir.iterdir() if p.is_dir())
    return [d for d in candidates if (d / f"{d.name}_classification").is_dir()]


def _legacy_outputs(asm_dir: Path) -> List[Path]:
    """Pre-bgzip-layout outputs whose current-name replacement exists."""
    prefix = asm_dir / asm_dir.name
    stale: List[Path] = []
    for old_suffix, new_suffix in _LEGACY_RENAMES.items():
        old = Path(f"{prefix}{old_suffix}")
        if old.exists() and Path(f"{prefix}{new_suffix}").exists():
            stale.append(old)
    for old in sorted(asm_dir.glob(f"{asm_dir.name}.*.fasta")):
        if old in stale or old.name.endswith(".rdna_consensus.fasta"):
            continue
        if Path(f"{old}.gz").exists():
            stale.append(old)
    return stale


def build_plan(
    results_dir: Path,
    sister_anchors: bool = False,
    legacy_outputs: bool = False,
    keep_bams: bool = False,
) -> CleanupPlan:
    """Collect the paths ``cleanup`` would remove under *results_dir*."""
    plan = CleanupPlan()
    for asm_dir in _assembly_dirs(results_dir):
        work = asm_dir / f"{asm_dir.name}_classification"
        for pattern in _WORKDIR_BULK:
            plan.add("work-dir intermediates", work.glob(pattern))
        if not keep_bams:
            for pattern in _READ_BAMS:
                plan.add("read-alignment BAMs", work.glob(pattern))
        if legacy_outputs:
            plan.add("legacy-name outputs", _legacy_outputs(asm_dir))

    # compleasm runs: per assembly, outgroups, reference
    for marker in results_dir.rglob("full_table_busco_format.tsv"):
        lineage_dir = marker.parent
        plan.add("compleasm raw output", (lineage_dir / n for n in _COMPLEASM_BULK))

    pairwise = results_dir / "pairwise"
    if pairwise.is_dir():
        plan.add("pairwise PAFs", pairwise.rglob("*.paf.gz"))

    if sister_anchors:
        sister = results_dir / "sister_anchoring"
        plan.add("sister-anchoring FASTAs", [sister / "anchors", sister / "candidates.fa"])
    return plan


def run_cleanup(
    results_dir: Path,
    dry_run: bool = False,
    sister_anchors: bool = False,
    legacy_outputs: bool = False,
    keep_bams: bool = False,
) -> int:
    """Remove (or list, with *dry_run*) bulky intermediates.  Returns bytes freed."""
    plan = build_plan(results_dir, sister_anchors, legacy_outputs, keep_bams)
    total = 0
    verb = "Would remove" if dry_run else "Removed"
    for label, paths in plan.groups.items():
        group_bytes = 0
        for p in paths:
            size = _size(p)
            group_bytes += size
            logger.debug(f"{verb} {p} ({_human(size)})")
            if not dry_run:
                _remove(p)
        total += group_bytes
        logger.info(f"{verb} {len(paths)} {label}: {_human(group_bytes)}")
    if not plan.groups:
        logger.info(f"Nothing to remove under {results_dir}")
    else:
        logger.done(f"{verb} {_human(total)} under {results_dir}")
    return total


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="dnadis cleanup",
        description=(
            "Remove bulky intermediates from a dnadis results folder. Outputs, "
            "report inputs and rerun caches are kept."
        ),
    )
    p.add_argument("results_dir", type=Path, metavar="RESULTS_DIR",
                   help="dnadis output folder (-o), or one assembly's subfolder")
    p.add_argument("-n", "--dry-run", action="store_true",
                   help="List what would be removed and how much space it frees, without removing it")
    p.add_argument("--keep-bams", action="store_true",
                   help="Keep read-alignment BAMs (written with --keep-depth-bam; "
                        "needed for split-read evidence on a rerun)")
    p.add_argument("--sister-anchors", action="store_true",
                   help="Also remove sister-anchoring anchor/candidate FASTAs "
                        "(a rerun with --sister-anchoring redoes the sister alignments)")
    p.add_argument("--legacy-outputs", action="store_true",
                   help="Also remove FASTA/AGP outputs under pre-bgzip names "
                        "(*.fasta, *.scaffolded.*, *.circular_elements.fasta) "
                        "where the current-name file exists")
    p.add_argument("-v", "--verbose", action="store_true", help="List every path")
    args = p.parse_args(argv)

    setup_logging(verbose=args.verbose)
    results_dir = args.results_dir.resolve()
    if not results_dir.is_dir():
        p.error(f"not a directory: {results_dir}")
    if not _assembly_dirs(results_dir):
        p.error(f"no dnadis assembly output (*/<name>_classification/) under {results_dir}")
    run_cleanup(results_dir, args.dry_run, args.sister_anchors, args.legacy_outputs, args.keep_bams)
    return 0
