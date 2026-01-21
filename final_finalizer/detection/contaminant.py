#!/usr/bin/env python3
"""
Contaminant detection for final_finalizer.

Contains functions for identifying contaminated contigs using centrifuger
taxonomic classification.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.models import ContaminantHit
from final_finalizer.utils.io_utils import have_exe, open_maybe_gzip


def _fasta_to_fastq_stream(fasta_path: Path, fastq_path: Path) -> None:
    """Convert FASTA to FASTQ with dummy quality scores.

    Centrifuger requires FASTQ input. For contigs, we use a constant high
    quality score (I = Phred 40) since these are assembled sequences.
    """
    with open_maybe_gzip(fasta_path, "rt") as fh_in, fastq_path.open("w") as fh_out:
        name: Optional[str] = None
        seq_parts: List[str] = []

        for line in fh_in:
            if line.startswith(">"):
                # Write previous record
                if name is not None and seq_parts:
                    seq = "".join(seq_parts)
                    fh_out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

        # Write last record
        if name is not None and seq_parts:
            seq = "".join(seq_parts)
            fh_out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def _validate_centrifuger_index(idx_prefix: str) -> bool:
    """Check if centrifuger index files exist."""
    idx_path = Path(idx_prefix)
    # Centrifuger index consists of .1.cfr, .2.cfr, .3.cfr files
    required_exts = [".1.cfr", ".2.cfr", ".3.cfr"]
    for ext in required_exts:
        if not Path(str(idx_path) + ext).exists():
            return False
    return True


def _get_centrifuger_name_table(idx_prefix: str, work_dir: Path) -> Dict[int, str]:
    """Get taxid -> scientific name mapping from centrifuger index.

    Runs centrifuger-inspect --name-table and caches the result.
    Returns dict mapping taxid (int) -> scientific name (str).
    """
    cache_path = work_dir / "centrifuger_names.tsv"

    # Use cached result if available
    if cache_path.exists():
        name_table: Dict[int, str] = {}
        with cache_path.open("r") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    try:
                        name_table[int(parts[0])] = parts[1]
                    except ValueError:
                        continue
        return name_table

    # Run centrifuger-inspect
    if not have_exe("centrifuger-inspect"):
        print("[warn] centrifuger-inspect not found, cannot look up scientific names", file=sys.stderr)
        return {}

    cmd = ["centrifuger-inspect", "-x", idx_prefix, "--name-table"]
    try:
        ret = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if ret.returncode != 0:
            print(f"[warn] centrifuger-inspect failed: {ret.stderr}", file=sys.stderr)
            return {}
    except Exception as e:
        print(f"[warn] centrifuger-inspect error: {e}", file=sys.stderr)
        return {}

    # Parse output: "taxid | name | scientific name |"
    name_table = {}
    for line in ret.stdout.splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) >= 2:
            try:
                taxid = int(parts[0])
                sci_name = parts[1]
                name_table[taxid] = sci_name
            except ValueError:
                continue

    # Cache result
    work_dir.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w") as fh:
        for taxid, name in sorted(name_table.items()):
            fh.write(f"{taxid}\t{name}\n")

    print(f"[info] Loaded {len(name_table)} taxid -> name mappings from centrifuger index", file=sys.stderr)
    return name_table


def detect_contaminants(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    centrifuger_idx: str,
    work_dir: Path,
    threads: int,
    min_score: int,
    exclude_contigs: Set[str],
) -> Dict[str, ContaminantHit]:
    """Identify contaminant contigs using centrifuger taxonomic classification.

    Centrifuger is much faster than BLAST for taxonomic classification.
    Any contig that gets a significant classification hit is considered
    a potential contaminant.

    Args:
        query_fasta: Query FASTA file
        query_lengths: Dict of contig name -> length
        centrifuger_idx: Path prefix to centrifuger index
        work_dir: Working directory for intermediate files
        threads: Number of threads
        min_score: Minimum centrifuger score to consider a hit significant
        exclude_contigs: Set of contigs to exclude from screening

    Returns:
        Dict mapping contig_name -> ContaminantHit with taxid, name, coverage, score
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Check for centrifuger executable
    if not have_exe("centrifuger"):
        print("[warn] centrifuger not found in PATH, skipping contaminant detection", file=sys.stderr)
        return {}

    # Validate index
    if not _validate_centrifuger_index(centrifuger_idx):
        print(f"[warn] centrifuger index not found: {centrifuger_idx}", file=sys.stderr)
        return {}

    print(f"[info] Using centrifuger index: {centrifuger_idx}", file=sys.stderr)

    # Convert FASTA to FASTQ (centrifuger requires FASTQ)
    fastq_path = work_dir / "contigs.fq"
    if not fastq_path.exists():
        print("[info] Converting FASTA to FASTQ for centrifuger", file=sys.stderr)
        _fasta_to_fastq_stream(query_fasta, fastq_path)

    # Run centrifuger
    output_path = work_dir / "centrifuger_results.tsv"
    err_path = work_dir / "centrifuger.err"

    if not output_path.exists():
        cmd = [
            "centrifuger",
            "-x", centrifuger_idx,
            "-u", str(fastq_path),
            "-t", str(threads),
        ]

        print(f"[info] Running centrifuger -> {output_path}", file=sys.stderr)

        err_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w") as out_fh, err_path.open("w") as err_fh:
            ret = subprocess.run(cmd, stdout=out_fh, stderr=err_fh, check=False)

        if ret.returncode != 0:
            print(f"[warn] centrifuger failed with return code {ret.returncode}", file=sys.stderr)
            return {}

    # Parse results
    contaminants: Dict[str, ContaminantHit] = {}

    if not output_path.exists() or output_path.stat().st_size == 0:
        return contaminants

    # Get taxid -> scientific name mapping
    name_table = _get_centrifuger_name_table(centrifuger_idx, work_dir)

    # Centrifuger output columns:
    # 1: Read ID, 2: Sequence ID, 3: Taxid, 4: Score, 5: Second-best score,
    # 6: Matching bp, 7: Read length, 8: Number of classifications
    with output_path.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue

            contig_name = fields[0]
            if contig_name in exclude_contigs:
                continue

            try:
                taxid = int(fields[2])
                score = int(fields[3])
                matching_bp = int(fields[5])
                read_len = int(fields[6])
            except ValueError:
                continue

            # Skip low-score hits
            if score < min_score:
                continue

            # Calculate coverage
            coverage = matching_bp / read_len if read_len > 0 else 0.0

            # Look up scientific name
            sci_name = name_table.get(taxid, "")

            contaminants[contig_name] = ContaminantHit(
                taxid=taxid,
                sci_name=sci_name,
                coverage=coverage,
                score=score,
            )
            contig_len = query_lengths.get(contig_name, read_len)
            print(
                f"[info] Contaminant: {contig_name} ({contig_len:,} bp, taxid={taxid}, {sci_name}, score={score}, cov={coverage:.2f})",
                file=sys.stderr,
            )

    print(f"[info] Contaminant contigs: {len(contaminants)}", file=sys.stderr)
    return contaminants
