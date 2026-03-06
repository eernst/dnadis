#!/usr/bin/env python3
"""
Consensus 45S rDNA pipeline for final_finalizer.

Builds a consensus nuclear 45S rDNA repeat unit from the query assembly,
annotates sub-features (18S, ITS1, 5.8S, ITS2, 25S/28S), and uses it
as a species-specific probe to:
1. Quantify rDNA load per chromosome-assigned contig
2. Detect NOR locations (nucleolar organizer regions)
3. Annotate sub-feature boundaries within each rDNA locus
4. Improve rDNA contig classification with a better-matched probe

Steps (within CLI phase 13):
  Step 1: Self-alignment → repeat boundary detection → copy extraction
  Step 2: Clustering → exemplar/consensus selection
  Step 3: Sub-feature annotation on the consensus
  Step 4: Re-annotation of all contigs with consensus probe
"""
from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.alignment.external_tools import (
    get_minimap2_exe,
    run_cdhit_est,
    run_mafft,
    run_minimap2,
)
from final_finalizer.detection.blast import (
    run_blastn_megablast,
    run_makeblastdb,
)
from final_finalizer.models import RdnaArray, RdnaConsensus, RdnaLocus, RdnaSubFeature, RdnaSubFeatureLocus
from final_finalizer.utils.io_utils import file_exists_and_valid, have_exe, merge_intervals
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.sequence_utils import (
    read_fasta_sequences,
    write_fasta,
)

logger = get_logger("rdna_consensus")


# ---------------------------------------------------------------------------
# Step 1: Extract rDNA-containing regions, self-align, detect repeat period,
#          extract individual 45S copies
# ---------------------------------------------------------------------------

def _extract_rdna_regions(
    query_fasta: Path,
    rdna_hit_intervals: Dict[str, List[Tuple[int, int]]],
    query_lengths: Dict[str, int],
    output_fasta: Path,
    flank_bp: int = 2000,
    min_region_bp: int = 5000,
) -> Dict[str, Tuple[str, int, int]]:
    """Extract rDNA-containing regions from the query assembly.

    Merges overlapping rDNA BLAST hits per contig, adds flanking context,
    and writes the extracted regions to a FASTA file.

    Args:
        query_fasta: Query assembly FASTA
        rdna_hit_intervals: Per-contig BLAST hit intervals from rDNA detection
        query_lengths: Contig lengths
        output_fasta: Output FASTA path for extracted regions
        flank_bp: Flanking context to add on each side [2000]
        min_region_bp: Minimum region size to extract [5000]

    Returns:
        Dict mapping region name -> (source_contig, region_start, region_end)
    """
    if file_exists_and_valid(output_fasta):
        logger.info(f"rDNA regions FASTA exists, reusing: {output_fasta}")
        # Reconstruct region map from FASTA names
        region_map = {}
        seqs = read_fasta_sequences(output_fasta)
        for name in seqs:
            # Parse name format: contig:start-end
            if ":" in name and "-" in name.split(":")[-1]:
                contig = name.rsplit(":", 1)[0]
                coords = name.rsplit(":", 1)[1]
                s, e = coords.split("-")
                region_map[name] = (contig, int(s), int(e))
        return region_map

    # Merge intervals per contig and add flanking
    regions_to_extract: Dict[str, List[Tuple[int, int]]] = {}
    for contig, intervals in rdna_hit_intervals.items():
        if not intervals:
            continue
        merged, total_bp = merge_intervals(intervals)
        if total_bp < min_region_bp:
            continue
        clen = query_lengths.get(contig, 0)
        # Add flanking and merge again
        flanked = []
        for s, e in merged:
            fs = max(0, s - flank_bp)
            fe = min(clen, e + flank_bp) if clen > 0 else e + flank_bp
            flanked.append((fs, fe))
        flanked_merged, _ = merge_intervals(flanked)
        regions_to_extract[contig] = flanked_merged

    if not regions_to_extract:
        logger.warning("No rDNA regions large enough to extract")
        return {}

    # Read sequences and extract regions
    all_seqs = read_fasta_sequences(query_fasta)
    extracted: Dict[str, str] = {}
    region_map: Dict[str, Tuple[str, int, int]] = {}

    for contig, regions in regions_to_extract.items():
        seq = all_seqs.get(contig, "")
        if not seq:
            continue
        for s, e in regions:
            region_name = f"{contig}:{s}-{e}"
            extracted[region_name] = seq[s:e]
            region_map[region_name] = (contig, s, e)

    if extracted:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(extracted, output_fasta)
        logger.info(f"Extracted {len(extracted)} rDNA regions -> {output_fasta}")

    return region_map


def _detect_repeat_period(
    regions_fasta: Path,
    work_dir: Path,
    threads: int,
    min_period: int = 8000,
    max_period: int = 15000,
) -> Optional[int]:
    """Detect the dominant repeat period via minimap2 self-alignment.

    Aligns the extracted rDNA regions against themselves, then builds
    a histogram of inter-hit distances on the same contig to find the
    dominant repeat unit length (expected ~10-13 kb for 45S).

    Args:
        regions_fasta: FASTA with extracted rDNA regions
        work_dir: Working directory for intermediate files
        threads: Number of threads
        min_period: Minimum expected repeat period [8000]
        max_period: Maximum expected repeat period [15000]

    Returns:
        Detected repeat period in bp, or None if not detected
    """
    paf_out = work_dir / "rdna_self_align.paf"

    mapper = get_minimap2_exe()
    if not mapper:
        logger.warning("minimap2 not available for self-alignment; using default period 10000")
        return 10000

    ok = run_minimap2(
        ref=regions_fasta,
        qry=regions_fasta,
        paf_out=paf_out,
        threads=threads,
        preset="asm5",
        extra_args=["-X"],  # self-mapping mode
        err_path=work_dir / "self_align.err",
    )

    if not ok or not file_exists_and_valid(paf_out):
        logger.warning("Self-alignment failed; using default period 10000")
        return 10000

    # Parse PAF: collect offsets of hits on the same query sequence
    distances: List[int] = []
    hits_per_seq: Dict[str, List[int]] = defaultdict(list)

    with paf_out.open("r") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            tname = fields[5]
            if qname != tname:
                continue
            try:
                qs = int(fields[2])
                qe = int(fields[3])
                ts = int(fields[7])
                te = int(fields[8])
            except ValueError:
                continue
            # Skip self-alignment (diagonal)
            if abs(qs - ts) < 1000:
                continue
            offset = abs(ts - qs)
            if min_period <= offset <= max_period:
                distances.append(offset)

    if not distances:
        logger.info("No repeat period detected from self-alignment; using default 10000")
        return 10000

    # Find dominant period: bin distances into 500 bp bins and take mode
    bin_size = 500
    binned = [d // bin_size * bin_size for d in distances]
    counter = Counter(binned)
    most_common_bin, count = counter.most_common(1)[0]

    # Refine: take median of distances in the most common bin
    in_bin = [d for d in distances if most_common_bin <= d < most_common_bin + bin_size]
    period = int(statistics.median(in_bin))

    logger.info(f"Detected rDNA repeat period: {period} bp (from {len(distances)} measurements)")
    return period


def _extract_individual_copies(
    regions_fasta: Path,
    region_map: Dict[str, Tuple[str, int, int]],
    repeat_period: int,
    output_fasta: Path,
    min_copy_frac: float = 0.7,
) -> int:
    """Extract individual 45S copies from rDNA regions using the detected period.

    Tiles each region into windows of repeat_period length and extracts
    copies that are at least min_copy_frac * repeat_period in length.

    Args:
        regions_fasta: FASTA with extracted rDNA regions
        region_map: Mapping of region name -> (contig, start, end)
        repeat_period: Detected repeat period in bp
        output_fasta: Output FASTA for individual copies
        min_copy_frac: Minimum fraction of period for a copy to be kept [0.7]

    Returns:
        Number of copies extracted
    """
    if file_exists_and_valid(output_fasta):
        seqs = read_fasta_sequences(output_fasta)
        logger.info(f"rDNA copies FASTA exists, reusing: {output_fasta} ({len(seqs)} copies)")
        return len(seqs)

    seqs = read_fasta_sequences(regions_fasta)
    min_len = int(repeat_period * min_copy_frac)
    copies: Dict[str, str] = {}
    copy_idx = 0

    for region_name, seq in seqs.items():
        region_len = len(seq)
        # Tile the region into copies
        pos = 0
        while pos + min_len <= region_len:
            end = min(pos + repeat_period, region_len)
            copy_seq = seq[pos:end]
            if len(copy_seq) >= min_len:
                source_contig = region_map.get(region_name, ("unknown", 0, 0))[0]
                copy_name = f"rdna_copy_{copy_idx}_{source_contig}_{pos}_{end}"
                copies[copy_name] = copy_seq
                copy_idx += 1
            pos += repeat_period

    if copies:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(copies, output_fasta)
        logger.info(f"Extracted {len(copies)} individual rDNA copies -> {output_fasta}")

    return len(copies)


# ---------------------------------------------------------------------------
# Step 2: Cluster copies and select exemplar/consensus
# ---------------------------------------------------------------------------

def _parse_cdhit_clusters(clstr_path: Path) -> Dict[int, List[str]]:
    """Parse cd-hit-est .clstr file to extract cluster membership.

    Returns:
        Dict mapping cluster_id -> list of sequence names
    """
    clusters: Dict[int, List[str]] = defaultdict(list)
    current_cluster = -1

    with clstr_path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[-1])
            elif current_cluster >= 0 and line:
                # Parse sequence name from line like: "0  9378nt, >name... *"
                parts = line.split(">")
                if len(parts) >= 2:
                    name = parts[1].split("...")[0]
                    clusters[current_cluster].append(name)

    return dict(clusters)


def _blast_all_vs_all_central(
    copies_fasta: Path,
    work_dir: Path,
    threads: int,
) -> Optional[str]:
    """Select the most central copy by BLAST all-vs-all identity.

    Fallback when cd-hit-est is not available.

    Returns:
        Name of the most central (highest average identity to all others) copy,
        or None if failed.
    """
    db_path = work_dir / "rdna_copies_db"
    run_makeblastdb(copies_fasta, db_path, err_path=work_dir / "makeblastdb_copies.err")

    blast_out = work_dir / "rdna_allvsall.txt"
    run_blastn_megablast(
        query_fasta=copies_fasta,
        db_paths=[str(db_path)],
        output_path=blast_out,
        threads=threads,
        max_hsps=1,
        max_target_seqs=500,
        err_path=work_dir / "blastn_allvsall.err",
    )

    if not file_exists_and_valid(blast_out):
        return None

    # Parse: accumulate identity per query
    identity_sums: Dict[str, float] = defaultdict(float)
    identity_counts: Dict[str, int] = defaultdict(int)

    with blast_out.open("r") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qseqid = fields[0]
            sseqid = fields[1]
            if qseqid == sseqid:
                continue  # Skip self-hit
            try:
                pident = float(fields[2])
            except ValueError:
                continue
            identity_sums[qseqid] += pident
            identity_counts[qseqid] += 1

    if not identity_sums:
        return None

    # Pick copy with highest mean identity to all others
    best_name = max(
        identity_sums.keys(),
        key=lambda n: identity_sums[n] / identity_counts[n] if identity_counts[n] > 0 else 0,
    )

    mean_ident = identity_sums[best_name] / identity_counts[best_name]
    logger.info(f"BLAST central copy: {best_name} (mean identity={mean_ident:.1f}%)")
    return best_name


def _mafft_consensus(
    aligned_fasta: Path,
    min_support: float = 0.50,
) -> Optional[str]:
    """Build majority-rule consensus from MAFFT-aligned sequences.

    For each alignment column, take the most frequent non-gap base.
    Columns where fewer than *min_support* fraction of sequences have
    a non-gap base are skipped.  This trims low-confidence flanking
    regions that arise when input copies start at different phases of
    the repeat unit.

    Args:
        aligned_fasta: MAFFT-aligned multi-FASTA.
        min_support: Minimum fraction of sequences that must have a
            non-gap base for a column to be included [0.50].

    Returns:
        Consensus sequence string, or None if no alignment.
    """
    seqs = read_fasta_sequences(aligned_fasta)
    if not seqs:
        return None

    sequences = list(seqs.values())
    if not sequences:
        return None

    n_seqs = len(sequences)
    aln_len = max(len(s) for s in sequences)
    min_bases = max(1, int(n_seqs * min_support))
    consensus_parts = []

    for i in range(aln_len):
        bases = []
        for seq in sequences:
            if i < len(seq):
                base = seq[i].upper()
                if base not in ("-", "."):
                    bases.append(base)
        if len(bases) >= min_bases:
            counter = Counter(bases)
            consensus_parts.append(counter.most_common(1)[0][0])
        # Skip low-support and pure-gap columns

    return "".join(consensus_parts)


def cluster_and_select_exemplar(
    copies_fasta: Path,
    work_dir: Path,
    threads: int,
    identity_threshold: float = 0.95,
) -> Tuple[Optional[str], str, int]:
    """Cluster rDNA copies and select exemplar or build consensus.

    Tries three strategies in order of preference:
    1. cd-hit-est + MAFFT: cluster, then MSA of largest cluster -> consensus
    2. cd-hit-est only: use cluster representative
    3. BLAST all-vs-all: pick most central copy

    Args:
        copies_fasta: FASTA with individual rDNA copies
        work_dir: Working directory
        threads: Number of threads
        identity_threshold: Clustering identity threshold [0.95]

    Returns:
        Tuple of (consensus_sequence, method_used, n_copies_in_cluster)
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    all_seqs = read_fasta_sequences(copies_fasta)
    if not all_seqs:
        return None, "none", 0

    if len(all_seqs) == 1:
        name = next(iter(all_seqs))
        return all_seqs[name], "single_copy", 1

    # Try cd-hit-est
    cdhit_prefix = work_dir / "rdna_cdhit"
    cdhit_ok = run_cdhit_est(
        input_fasta=copies_fasta,
        output_prefix=cdhit_prefix,
        identity=identity_threshold,
        threads=threads,
        err_path=work_dir / "cdhit.err",
    )

    if cdhit_ok:
        clstr_path = Path(str(cdhit_prefix) + ".clstr")
        if file_exists_and_valid(clstr_path):
            clusters = _parse_cdhit_clusters(clstr_path)
            if clusters:
                # Find largest cluster
                largest_id = max(clusters, key=lambda k: len(clusters[k]))
                largest_members = clusters[largest_id]
                n_clustered = len(largest_members)
                logger.info(f"Largest CD-HIT-EST cluster: {n_clustered} copies")

                # Try MAFFT consensus on largest cluster
                if n_clustered >= 3:
                    cluster_fasta = work_dir / "largest_cluster.fa"
                    cluster_seqs = {n: all_seqs[n] for n in largest_members if n in all_seqs}
                    if cluster_seqs:
                        write_fasta(cluster_seqs, cluster_fasta)

                        aligned_fasta = work_dir / "largest_cluster_aligned.fa"
                        mafft_ok = run_mafft(
                            input_fasta=cluster_fasta,
                            output_fasta=aligned_fasta,
                            threads=threads,
                            err_path=work_dir / "mafft.err",
                        )
                        if mafft_ok and file_exists_and_valid(aligned_fasta):
                            consensus = _mafft_consensus(aligned_fasta)
                            if consensus and len(consensus) > 1000:
                                logger.info(f"MAFFT consensus: {len(consensus)} bp")
                                return consensus, "cdhit+mafft", n_clustered

                # Fallback: use cd-hit-est representative
                rep_seqs = read_fasta_sequences(cdhit_prefix)
                if rep_seqs:
                    # The representative of the largest cluster should be in the rep file
                    # Try to find a member of the largest cluster in representatives
                    for member in largest_members:
                        if member in rep_seqs:
                            logger.info(f"Using CD-HIT-EST representative: {member}")
                            return rep_seqs[member], "cdhit_rep", n_clustered
                    # Fallback: just use the first representative
                    first_name = next(iter(rep_seqs))
                    return rep_seqs[first_name], "cdhit_rep", n_clustered

    # Fallback: BLAST all-vs-all centrality selection
    logger.info("cd-hit-est unavailable or failed; using BLAST all-vs-all centrality")
    central_name = _blast_all_vs_all_central(copies_fasta, work_dir, threads)
    if central_name and central_name in all_seqs:
        return all_seqs[central_name], "blast_central", len(all_seqs)

    # Last resort: just pick the longest copy
    longest = max(all_seqs.keys(), key=lambda n: len(all_seqs[n]))
    return all_seqs[longest], "longest_copy", len(all_seqs)


# ---------------------------------------------------------------------------
# Step 3: Sub-feature annotation using Infernal/Rfam
# ---------------------------------------------------------------------------

def _find_rfam_cm_database() -> Optional[Path]:
    """Find the bundled Rfam covariance model database.

    Auto-presses the database if index files are missing.

    Returns:
        Path to euk-rrna.cm if found and pressed, else None.
    """
    import subprocess

    # Look relative to this module
    module_dir = Path(__file__).resolve().parent.parent
    cm_path = module_dir / "data" / "rfam" / "euk-rrna.cm"

    if not cm_path.exists():
        return None

    # Check if pressed (index files exist)
    i1m_path = Path(str(cm_path) + ".i1m")
    if i1m_path.exists():
        return cm_path

    # Auto-press the database
    if not have_exe("cmpress"):
        logger.warning(f"Rfam CM found but cmpress not available to create index")
        return None

    logger.info(f"Pressing Rfam CM database: {cm_path}")
    try:
        ret = subprocess.call(
            ["cmpress", str(cm_path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if ret == 0 and i1m_path.exists():
            logger.info("Rfam CM database pressed successfully")
            return cm_path
        else:
            logger.warning(f"cmpress failed with return code {ret}")
            return None
    except Exception as e:
        logger.warning(f"Failed to press CM database: {e}")
        return None


def annotate_sub_features_infernal(
    consensus_seq: str,
    work_dir: Path,
    threads: int,
    cm_database: Optional[Path] = None,
) -> List[RdnaSubFeature]:
    """Annotate rRNA sub-features using Infernal covariance models.

    Uses cmscan against Rfam models for structure-based annotation,
    providing more accurate boundaries than homology-based projection.

    The bundled Rfam database contains:
    - RF00001: 5S_rRNA (not in 45S, but included for completeness)
    - RF00002: 5_8S_rRNA
    - RF01960: SSU_rRNA_eukarya (18S)
    - RF02543: LSU_rRNA_eukarya (25S/28S)

    ITS1 and ITS2 are derived as the gaps between rRNA genes.

    Args:
        consensus_seq: The consensus 45S sequence
        work_dir: Working directory
        threads: Number of threads
        cm_database: Optional path to CM database (auto-detected if None)

    Returns:
        List of RdnaSubFeature with coordinates on the consensus,
        or empty list if Infernal is unavailable or fails.
    """
    import subprocess

    # Check for cmscan
    if not have_exe("cmscan"):
        logger.info("cmscan not available; will use BLAST-based annotation")
        return []

    # Find CM database
    if cm_database is None:
        cm_database = _find_rfam_cm_database()
    if cm_database is None:
        logger.warning("Rfam CM database not found; will use BLAST-based annotation")
        return []

    work_dir.mkdir(parents=True, exist_ok=True)

    # Write consensus sequence
    consensus_fasta = work_dir / "consensus_for_infernal.fa"
    write_fasta({"rdna_consensus": consensus_seq}, consensus_fasta)

    # Run cmscan
    tbl_out = work_dir / "cmscan_rrna.tbl"
    cmscan_out = work_dir / "cmscan_rrna.out"
    err_path = work_dir / "cmscan.err"

    cmd = [
        "cmscan",
        "--cpu", str(threads),
        "--tblout", str(tbl_out),
        "-o", str(cmscan_out),
        "--noali",  # Skip alignment output for speed
        str(cm_database),
        str(consensus_fasta),
    ]

    logger.info(f"Running Infernal cmscan for sub-feature annotation")

    with err_path.open("wb") as err_fh:
        ret = subprocess.call(cmd, stderr=err_fh)

    if ret != 0:
        logger.warning(f"cmscan failed with return code {ret}; falling back to BLAST")
        return []

    if not file_exists_and_valid(tbl_out):
        logger.warning("cmscan produced no output; falling back to BLAST")
        return []

    # Parse tabular output
    # Format: target_name accession query_name accession mdl mdl_from mdl_to seq_from seq_to strand ...
    hits: Dict[str, Tuple[int, int, float]] = {}  # model_name -> (start, end, score)

    # Model name mappings to our feature names
    MODEL_TO_FEATURE = {
        "SSU_rRNA_eukarya": "18S",
        "5_8S_rRNA": "5.8S",
        "LSU_rRNA_eukarya": "25S",
        "5S_rRNA": "5S",  # Not in 45S but may appear
    }

    with tbl_out.open("r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 15:
                continue

            try:
                model_name = fields[0]
                seq_from = int(fields[7])
                seq_to = int(fields[8])
                strand = fields[9]
                score = float(fields[14])
                inc = fields[16]  # '!' for significant, '?' for below threshold
            except (ValueError, IndexError):
                continue

            # Only use significant hits
            if inc != "!":
                continue

            # Skip 5S (not part of 45S)
            if model_name == "5S_rRNA":
                continue

            feature_name = MODEL_TO_FEATURE.get(model_name)
            if feature_name is None:
                continue

            # Normalize coordinates (cmscan outputs 1-based)
            start = min(seq_from, seq_to) - 1  # Convert to 0-based
            end = max(seq_from, seq_to)

            # Keep best hit per model
            prev = hits.get(feature_name)
            if prev is None or score > prev[2]:
                hits[feature_name] = (start, end, score)

    if not hits:
        logger.warning("No significant rRNA hits from cmscan")
        return []

    # Build feature list
    result: List[RdnaSubFeature] = []

    for name in ("18S", "5.8S", "25S"):
        if name in hits:
            start, end, score = hits[name]
            result.append(RdnaSubFeature(name=name, start=start, end=end))
            logger.info(f"  {name:5s}: {start}-{end} ({end - start} bp, score={score:.1f})")

    if not result:
        return []

    # Sort by position
    result.sort(key=lambda f: f.start)

    # Derive ITS regions from gaps between rRNA genes
    rrna_by_name = {f.name: f for f in result}
    its_features: List[RdnaSubFeature] = []

    if "18S" in rrna_by_name and "5.8S" in rrna_by_name:
        f18s = rrna_by_name["18S"]
        f58s = rrna_by_name["5.8S"]
        if f18s.end < f58s.start:
            its_features.append(RdnaSubFeature(
                name="ITS1", start=f18s.end, end=f58s.start,
            ))
            logger.info(f"  ITS1 : {f18s.end}-{f58s.start} ({f58s.start - f18s.end} bp, derived)")
        elif f58s.end < f18s.start:
            # Reverse order on consensus
            its_features.append(RdnaSubFeature(
                name="ITS1", start=f58s.end, end=f18s.start,
            ))
            logger.info(f"  ITS1 : {f58s.end}-{f18s.start} ({f18s.start - f58s.end} bp, derived)")

    if "5.8S" in rrna_by_name and "25S" in rrna_by_name:
        f58s = rrna_by_name["5.8S"]
        f25s = rrna_by_name["25S"]
        if f58s.end < f25s.start:
            its_features.append(RdnaSubFeature(
                name="ITS2", start=f58s.end, end=f25s.start,
            ))
            logger.info(f"  ITS2 : {f58s.end}-{f25s.start} ({f25s.start - f58s.end} bp, derived)")
        elif f25s.end < f58s.start:
            its_features.append(RdnaSubFeature(
                name="ITS2", start=f25s.end, end=f58s.start,
            ))
            logger.info(f"  ITS2 : {f25s.end}-{f58s.start} ({f58s.start - f25s.end} bp, derived)")

    result.extend(its_features)
    result.sort(key=lambda f: f.start)

    logger.info(f"Infernal annotated {len(result)} sub-features on consensus")
    return result


# ---------------------------------------------------------------------------
# Step 4: Re-annotation of all contigs with the consensus probe
# ---------------------------------------------------------------------------

def annotate_contigs_with_consensus(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    consensus: RdnaConsensus,
    work_dir: Path,
    threads: int,
    classifications: Optional[Dict[str, str]] = None,
    min_tandem_copies: int = 3,
    max_tandem_gap: int = 50000,
) -> Tuple[List[RdnaLocus], List[RdnaArray]]:
    """BLAST the consensus 45S against all contigs to annotate rDNA loci.

    Args:
        query_fasta: Query assembly FASTA (all contigs)
        query_lengths: Contig lengths
        consensus: RdnaConsensus object with sequence and sub-features
        work_dir: Working directory
        threads: Number of threads
        classifications: Optional dict of contig_name -> classification string
        min_tandem_copies: Minimum tandem copies for array detection [3]
        max_tandem_gap: Maximum gap between tandem copies [50000]

    Returns:
        Tuple of (list of RdnaLocus annotations, list of RdnaArray objects)
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Write consensus FASTA
    consensus_fasta = work_dir / "rdna_consensus_probe.fa"
    write_fasta({"rdna_consensus": consensus.sequence}, consensus_fasta)

    # Create BLAST db from query assembly
    query_db = work_dir / "query_db"
    run_makeblastdb(query_fasta, query_db, err_path=work_dir / "makeblastdb_query.err")

    # BLAST consensus against all contigs
    blast_out = work_dir / "consensus_vs_contigs.txt"
    run_blastn_megablast(
        query_fasta=consensus_fasta,
        db_paths=[str(query_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=500,
        max_target_seqs=1000,
        err_path=work_dir / "blastn_consensus_vs_contigs.err",
    )

    if not file_exists_and_valid(blast_out):
        logger.warning("BLAST of consensus vs contigs failed")
        return [], []

    # Parse BLAST hits
    # Fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    hits_per_contig: Dict[str, List[dict]] = defaultdict(list)

    with blast_out.open("r") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                sseqid = fields[1]  # Subject = contig
                pident = float(fields[2])
                aln_length = int(fields[3])
                qstart = int(fields[6])
                qend = int(fields[7])
                sstart = int(fields[8])
                send = int(fields[9])
            except ValueError:
                continue

            # Determine strand and normalize coordinates
            if sstart <= send:
                strand = "+"
                s_low, s_high = sstart, send
            else:
                strand = "-"
                s_low, s_high = send, sstart

            if qstart > qend:
                qstart, qend = qend, qstart

            hits_per_contig[sseqid].append({
                "pident": pident,
                "aln_length": aln_length,
                "qstart": qstart,
                "qend": qend,
                "sstart": s_low,
                "send": s_high,
                "strand": strand,
            })

    # Process hits into loci
    loci: List[RdnaLocus] = []

    for contig, hits in hits_per_contig.items():
        # Sort by subject start position
        hits.sort(key=lambda h: h["sstart"])

        # Merge overlapping/adjacent hits into loci
        merged_loci = _merge_hits_into_loci(hits, consensus.length, consensus.sub_features)

        for locus in merged_loci:
            locus_dict = dict(locus)
            loci.append(RdnaLocus(
                contig=contig,
                start=locus_dict["start"],
                end=locus_dict["end"],
                strand=locus_dict["strand"],
                identity=locus_dict["identity"],
                consensus_coverage=locus_dict["consensus_coverage"],
                copy_type=locus_dict["copy_type"],
                sub_feature_loci=locus_dict["sub_feature_loci"],
            ))

    # Detect tandem 45S rDNA arrays
    arrays = _detect_arrays(loci, min_tandem_copies, max_tandem_gap)

    logger.info(f"Annotated {len(loci)} rDNA loci across {len(hits_per_contig)} contigs")
    if arrays:
        array_loci = sum(a.n_total for a in arrays)
        logger.info(f"45S rDNA arrays: {len(arrays)} arrays with {array_loci} loci")

    return loci, arrays


def _merge_hits_into_loci(
    hits: List[dict],
    consensus_length: int,
    sub_features: List[RdnaSubFeature],
    merge_gap: int = 2000,
) -> List[dict]:
    """Merge overlapping BLAST hits into rDNA loci and classify them.

    Args:
        hits: Sorted list of BLAST hit dicts
        consensus_length: Length of the consensus sequence
        sub_features: Sub-feature annotations on the consensus
        merge_gap: Maximum gap between hits to merge [2000]

    Returns:
        List of locus dicts with classification info
    """
    if not hits:
        return []

    # Merge overlapping hits on the subject (contig)
    merged_groups: List[List[dict]] = []
    current_group = [hits[0]]

    for hit in hits[1:]:
        if hit["sstart"] <= current_group[-1]["send"] + merge_gap:
            current_group.append(hit)
        else:
            merged_groups.append(current_group)
            current_group = [hit]
    merged_groups.append(current_group)

    loci = []
    for group in merged_groups:
        start = min(h["sstart"] for h in group)
        end = max(h["send"] for h in group)

        # Determine strand (majority vote)
        strand_votes = Counter(h["strand"] for h in group)
        strand = strand_votes.most_common(1)[0][0]

        # Compute identity (weighted by alignment length)
        total_matches = sum(h["pident"] * h["aln_length"] / 100.0 for h in group)
        total_alnlen = sum(h["aln_length"] for h in group)
        identity = total_matches / total_alnlen if total_alnlen > 0 else 0.0

        # Compute consensus coverage (fraction of consensus covered)
        query_intervals = [(h["qstart"], h["qend"]) for h in group]
        _, covered_bp = merge_intervals(query_intervals)
        consensus_coverage = covered_bp / consensus_length if consensus_length > 0 else 0.0

        # Classify copy type based on coverage of the transcribed region
        # (18S through 25S) rather than the full consensus including IGS.
        # The IGS is the most rapidly evolving part of the 45S unit and
        # often fails to align via BLAST, making functionally complete
        # copies appear as "partial" if scored against full consensus length.
        if sub_features:
            tx_start = min(f.start for f in sub_features)
            tx_end = max(f.end for f in sub_features)
            tx_len = tx_end - tx_start
            if tx_len > 0:
                # Compute coverage restricted to the transcribed region
                tx_intervals = []
                for qs, qe in query_intervals:
                    # Clip to transcribed region
                    cs = max(qs, tx_start)
                    ce = min(qe, tx_end)
                    if ce > cs:
                        tx_intervals.append((cs, ce))
                _, tx_covered = merge_intervals(tx_intervals) if tx_intervals else ([], 0)
                tx_coverage = tx_covered / tx_len
            else:
                tx_coverage = consensus_coverage
        else:
            tx_coverage = consensus_coverage

        if tx_coverage >= 0.90:
            copy_type = "full"
        elif tx_coverage >= 0.50:
            copy_type = "partial"
        else:
            copy_type = "fragment"

        # Determine which sub-features are present and map their coordinates
        # to the contig using the BLAST hit alignment
        sub_feature_loci: List[RdnaSubFeatureLocus] = []
        for feat in sub_features:
            # Collect contig coordinate intervals for this feature across all hits
            contig_intervals = []
            for h in group:
                if h["qstart"] <= feat.end and h["qend"] >= feat.start:
                    # Clip feature coords to this hit's consensus span
                    clip_start = max(feat.start, h["qstart"])
                    clip_end = min(feat.end, h["qend"])
                    if clip_end <= clip_start:
                        continue

                    # Map consensus coordinates to contig coordinates
                    # The mapping depends on strand orientation
                    if h["strand"] == "+":
                        # Forward strand: contig coords increase with consensus coords
                        contig_start = h["sstart"] + (clip_start - h["qstart"])
                        contig_end = h["sstart"] + (clip_end - h["qstart"])
                    else:
                        # Reverse strand: contig coords decrease with consensus coords
                        # For "-" strand, sstart < send but maps to reversed consensus
                        contig_end = h["send"] - (clip_start - h["qstart"])
                        contig_start = h["send"] - (clip_end - h["qstart"])

                    if contig_end > contig_start:
                        contig_intervals.append((contig_start, contig_end))

            if contig_intervals:
                # Merge overlapping intervals
                merged, total_bp = merge_intervals(contig_intervals)
                feat_len = feat.end - feat.start

                # Check if we have at least 50% coverage of the feature
                if feat_len > 0 and total_bp / feat_len >= 0.5:
                    # Use the full merged extent as the feature bounds
                    merged_start = min(iv[0] for iv in merged)
                    merged_end = max(iv[1] for iv in merged)
                    sub_feature_loci.append(RdnaSubFeatureLocus(
                        name=feat.name,
                        start=merged_start,
                        end=merged_end,
                    ))

        # Sort sub-features by start position
        sub_feature_loci.sort(key=lambda sf: sf.start)

        loci.append({
            "start": start,
            "end": end,
            "strand": strand,
            "identity": identity,
            "consensus_coverage": consensus_coverage,
            "copy_type": copy_type,
            "sub_feature_loci": sub_feature_loci,
        })

    return loci


def _build_array(
    array_id: str,
    contig: str,
    cluster: List[RdnaLocus],
) -> RdnaArray:
    """Construct an RdnaArray from a cluster of loci.

    Computes copy type counts, identity statistics, and strand majority vote.
    Sets array_id on each constituent locus.
    """
    for locus in cluster:
        locus.array_id = array_id

    n_full = sum(1 for l in cluster if l.copy_type == "full")
    n_partial = sum(1 for l in cluster if l.copy_type == "partial")
    n_fragment = sum(1 for l in cluster if l.copy_type == "fragment")

    identities = [l.identity for l in cluster]
    identity_median = statistics.median(identities)
    identity_min = min(identities)
    identity_max = max(identities)

    # Strand majority vote
    strand_votes = Counter(l.strand for l in cluster)
    strand = strand_votes.most_common(1)[0][0]

    return RdnaArray(
        array_id=array_id,
        contig=contig,
        start=cluster[0].start,
        end=cluster[-1].end,
        strand=strand,
        loci=list(cluster),
        n_total=len(cluster),
        n_full=n_full,
        n_partial=n_partial,
        n_fragment=n_fragment,
        identity_median=identity_median,
        identity_min=identity_min,
        identity_max=identity_max,
    )


def _detect_arrays(
    loci: List[RdnaLocus],
    min_tandem_copies: int,
    max_tandem_gap: int,
) -> List[RdnaArray]:
    """Detect tandem 45S rDNA repeat arrays from loci.

    Groups loci by contig, sorts by start position, clusters by gap,
    and builds RdnaArray objects for clusters with enough copies.
    Sets array_id on constituent loci.

    Args:
        loci: List of RdnaLocus annotations
        min_tandem_copies: Minimum copies for an array
        max_tandem_gap: Maximum gap between adjacent loci in an array

    Returns:
        List of RdnaArray objects, numbered sequentially.
    """
    arrays: List[RdnaArray] = []
    array_counter = 0

    # Group by contig
    by_contig: Dict[str, List[RdnaLocus]] = defaultdict(list)
    for locus in loci:
        by_contig[locus.contig].append(locus)

    # Process contigs in sorted order for deterministic array numbering
    for contig in sorted(by_contig.keys()):
        contig_loci = by_contig[contig]
        contig_loci.sort(key=lambda l: l.start)

        # Find tandem clusters
        cluster: List[RdnaLocus] = [contig_loci[0]]
        for locus in contig_loci[1:]:
            if locus.start - cluster[-1].end <= max_tandem_gap:
                cluster.append(locus)
            else:
                if len(cluster) >= min_tandem_copies:
                    array_counter += 1
                    arrays.append(_build_array(
                        f"array_{array_counter}", contig, cluster,
                    ))
                cluster = [locus]
        # Check last cluster
        if len(cluster) >= min_tandem_copies:
            array_counter += 1
            arrays.append(_build_array(
                f"array_{array_counter}", contig, cluster,
            ))

    return arrays


# ---------------------------------------------------------------------------
# Step 4b: Reclassify contigs using consensus-based coverage
# ---------------------------------------------------------------------------

def identify_rdna_contigs_from_loci(
    loci: List[RdnaLocus],
    query_lengths: Dict[str, int],
    min_coverage: float = 0.50,
    exclude_contigs: Optional[Set[str]] = None,
) -> Tuple[Set[str], Dict[str, float]]:
    """Identify contigs with significant rDNA content using consensus-based loci.

    Computes per-contig rDNA coverage by merging all locus intervals and
    comparing to contig length. This uses the species-specific consensus
    probe, so it's more sensitive than the initial seed-based detection
    for divergent species.

    Args:
        loci: List of RdnaLocus annotations from consensus re-annotation
        query_lengths: Contig lengths
        min_coverage: Minimum rDNA coverage fraction for reclassification [0.50]
        exclude_contigs: Contigs to exclude (e.g., already classified as chromosomes)

    Returns:
        Tuple of:
        - Set of contig names passing coverage threshold
        - Dict mapping contig name -> rDNA coverage fraction (for all contigs with hits)
    """
    exclude = exclude_contigs or set()

    # Collect intervals per contig
    contig_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for locus in loci:
        contig_intervals[locus.contig].append((locus.start, locus.end))

    rdna_contigs: Set[str] = set()
    coverage_map: Dict[str, float] = {}

    for contig, intervals in contig_intervals.items():
        if contig in exclude:
            continue
        clen = query_lengths.get(contig, 0)
        if clen <= 0:
            continue
        _, total_bp = merge_intervals(intervals)
        cov = total_bp / clen
        coverage_map[contig] = cov
        if cov >= min_coverage:
            rdna_contigs.add(contig)

    return rdna_contigs, coverage_map


# ---------------------------------------------------------------------------
# Main orchestrator: build_rdna_consensus
# ---------------------------------------------------------------------------

def build_rdna_consensus(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    rdna_hit_intervals: Dict[str, List[Tuple[int, int]]],
    seed_ref_path: Path,
    work_dir: Path,
    threads: int,
    rdna_ref_features_path: Optional[Path] = None,
    classifications: Optional[Dict[str, str]] = None,
    min_tandem_copies: int = 3,
) -> Tuple[Optional[RdnaConsensus], List[RdnaLocus], List[RdnaArray]]:
    """Build consensus 45S rDNA and annotate all contigs.

    This is the main entry point that orchestrates Phases 2-5.

    Args:
        query_fasta: Query assembly FASTA
        query_lengths: Contig lengths
        rdna_hit_intervals: Per-contig BLAST hit intervals from rDNA detection
        seed_ref_path: Seed rDNA reference FASTA (Arabidopsis or user-provided)
        work_dir: Working directory for all intermediate files
        threads: Number of threads
        rdna_ref_features_path: Optional path to sub-feature TSV for seed reference
        classifications: Optional dict of contig_name -> classification
        min_tandem_copies: Minimum tandem copies for array detection [3]

    Returns:
        Tuple of (RdnaConsensus or None, list of RdnaLocus annotations, list of RdnaArray objects)
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    logger.phase("rDNA consensus: Step 1 - Extract rDNA regions and copies")

    if not rdna_hit_intervals:
        logger.warning("No rDNA hit intervals available; skipping consensus building")
        return None, [], []

    # Step 1a: Extract rDNA-containing regions
    regions_fasta = work_dir / "rdna_regions.fa"
    region_map = _extract_rdna_regions(
        query_fasta=query_fasta,
        rdna_hit_intervals=rdna_hit_intervals,
        query_lengths=query_lengths,
        output_fasta=regions_fasta,
    )

    if not region_map:
        logger.warning("No rDNA regions extracted; skipping consensus building")
        return None, [], []

    # Step 1b: Detect repeat period via self-alignment
    repeat_period = _detect_repeat_period(
        regions_fasta=regions_fasta,
        work_dir=work_dir / "self_align",
        threads=threads,
    )

    if not repeat_period:
        logger.warning("Could not determine repeat period; skipping consensus building")
        return None, [], []

    # Step 1c: Extract individual copies
    copies_fasta = work_dir / "rdna_copies.fa"
    n_copies = _extract_individual_copies(
        regions_fasta=regions_fasta,
        region_map=region_map,
        repeat_period=repeat_period,
        output_fasta=copies_fasta,
    )

    if n_copies == 0:
        logger.warning("No rDNA copies extracted; skipping consensus building")
        return None, [], []

    # Step 2: Cluster and select exemplar/consensus
    logger.phase("rDNA consensus: Step 2 - Cluster copies and build consensus")
    consensus_seq, method, n_clustered = cluster_and_select_exemplar(
        copies_fasta=copies_fasta,
        work_dir=work_dir / "clustering",
        threads=threads,
    )

    if not consensus_seq:
        logger.warning("Failed to build consensus; skipping rDNA annotation")
        return None, [], []

    logger.info(f"Consensus: {len(consensus_seq)} bp, method={method}, "
                f"{n_copies} copies extracted, {n_clustered} in cluster")

    # Step 3: Sub-feature annotation using Infernal/Rfam
    logger.phase("rDNA consensus: Step 3 - Annotate sub-features (Infernal)")
    sub_features = annotate_sub_features_infernal(
        consensus_seq=consensus_seq,
        work_dir=work_dir / "annotation",
        threads=threads,
    )

    if not sub_features:
        logger.warning(
            "Sub-feature annotation unavailable (requires Infernal: "
            "conda install -c bioconda infernal)"
        )

    # Build RdnaConsensus object
    consensus = RdnaConsensus(
        sequence=consensus_seq,
        length=len(consensus_seq),
        n_copies_extracted=n_copies,
        n_copies_clustered=n_clustered,
        method=method,
        sub_features=sub_features,
    )

    # Write consensus FASTA
    consensus_out = work_dir / "rdna_consensus.fa"
    write_fasta({"rdna_consensus": consensus_seq}, consensus_out)

    # Step 4: Re-annotate all contigs
    logger.phase("rDNA consensus: Step 4 - Annotate all contigs")
    loci, arrays = annotate_contigs_with_consensus(
        query_fasta=query_fasta,
        query_lengths=query_lengths,
        consensus=consensus,
        work_dir=work_dir / "reannotation",
        threads=threads,
        classifications=classifications,
        min_tandem_copies=min_tandem_copies,
    )

    return consensus, loci, arrays
