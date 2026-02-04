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

Pipeline phases:
  Phase 2: Self-alignment → repeat boundary detection → copy extraction
  Phase 3: Clustering → exemplar/consensus selection
  Phase 4: Sub-feature annotation on the consensus
  Phase 5: Re-annotation of all contigs with consensus probe
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
from final_finalizer.models import RdnaConsensus, RdnaLocus, RdnaSubFeature
from final_finalizer.utils.io_utils import file_exists_and_valid, have_exe, merge_intervals
from final_finalizer.utils.logging_config import get_logger
from final_finalizer.utils.sequence_utils import (
    read_fasta_sequences,
    write_fasta,
)

logger = get_logger("rdna_consensus")


# ---------------------------------------------------------------------------
# Phase 2: Extract rDNA-containing regions, self-align, detect repeat period,
#           extract individual 45S copies
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
        rdna_hit_intervals: Per-contig BLAST hit intervals from Phase 1
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
# Phase 3: Cluster copies and select exemplar/consensus
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
# Phase 4: Sub-feature annotation on the consensus
# ---------------------------------------------------------------------------

def _load_features_tsv(features_path: Path) -> List[RdnaSubFeature]:
    """Load sub-feature coordinates from a TSV file.

    Format: feature_name<tab>start<tab>end (1-based inclusive)

    Returns:
        List of RdnaSubFeature objects (converted to 0-based half-open)
    """
    features = []
    with features_path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            try:
                start = int(parts[1]) - 1  # Convert to 0-based
                end = int(parts[2])  # Already exclusive in 0-based
            except ValueError:
                continue
            features.append(RdnaSubFeature(name=name, start=start, end=end))
    return features


def _find_default_features_tsv(script_dir: Path) -> Optional[Path]:
    """Find the bundled Arabidopsis 45S features TSV."""
    search_paths = [
        script_dir / "data" / "athal-45s-features.tsv",
        script_dir.parent / "data" / "athal-45s-features.tsv",
    ]
    for path in search_paths:
        if path.exists():
            return path
    return None


def annotate_sub_features(
    consensus_seq: str,
    seed_ref_path: Path,
    seed_features: List[RdnaSubFeature],
    work_dir: Path,
    threads: int,
) -> List[RdnaSubFeature]:
    """Map sub-features from the seed reference onto the consensus.

    Extracts the three rRNA subunit sequences (18S, 5.8S, 25S) from
    the seed reference and BLASTs each individually against the
    consensus.  This gives direct, per-subunit coordinate placement
    that works across divergent species (unlike the previous approach
    of linear projection from a single alignment).

    ITS1 and ITS2 are defined as the gaps between the three subunits:
      ITS1 = 18S end  →  5.8S start
      ITS2 = 5.8S end →  25S start

    Uses ``blastn -task blastn`` (word size 11) instead of megablast
    so that even the short 5.8S (~164 bp) aligns reliably across
    kingdoms.

    Args:
        consensus_seq: The consensus 45S sequence
        seed_ref_path: Path to the seed rDNA reference FASTA
        seed_features: Sub-feature coordinates on the seed reference
        work_dir: Working directory
        threads: Number of threads

    Returns:
        List of RdnaSubFeature with coordinates on the consensus
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Identify the three rRNA subunit features from the seed features
    rRNA_NAMES = ("18S", "5.8S", "25S")
    subunit_feats = {f.name: f for f in seed_features if f.name in rRNA_NAMES}

    if len(subunit_feats) < 3:
        logger.warning(
            f"Seed features must include 18S, 5.8S, and 25S; "
            f"found {list(subunit_feats.keys())}"
        )
        return []

    # Read seed reference sequence
    seed_seqs = read_fasta_sequences(seed_ref_path)
    if not seed_seqs:
        logger.warning("Could not read seed reference FASTA")
        return []
    seed_seq = next(iter(seed_seqs.values()))

    # Extract each subunit sequence from the seed and write as a
    # multi-FASTA query
    subunit_fasta = work_dir / "subunit_queries.fa"
    subunit_seqs = {}
    for name in rRNA_NAMES:
        feat = subunit_feats[name]
        subunit_seqs[name] = seed_seq[feat.start : feat.end]
    write_fasta(subunit_seqs, subunit_fasta)

    # Write consensus as BLAST database
    consensus_fasta = work_dir / "consensus_for_annotation.fa"
    write_fasta({"rdna_consensus": consensus_seq}, consensus_fasta)

    db_path = work_dir / "consensus_db"
    run_makeblastdb(
        consensus_fasta, db_path,
        err_path=work_dir / "makeblastdb_consensus.err",
    )

    # BLAST subunits against consensus using blastn (not megablast)
    # for cross-species sensitivity
    blast_out = work_dir / "subunits_vs_consensus.txt"
    _run_blastn_subunits(
        query_fasta=subunit_fasta,
        db_path=db_path,
        output_path=blast_out,
        threads=threads,
        err_path=work_dir / "blastn_subunits.err",
    )

    if not file_exists_and_valid(blast_out):
        logger.warning("BLAST subunit annotation failed")
        return []

    # Parse best hit per subunit (highest bitscore)
    best_hits: Dict[str, Tuple[int, int, float]] = {}  # name -> (cons_start, cons_end, bitscore)

    with blast_out.open("r") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                qname = fields[0]
                sstart = int(fields[8])
                send = int(fields[9])
                bitscore = float(fields[11])
            except (ValueError, IndexError):
                continue

            if qname not in rRNA_NAMES:
                continue

            s_lo = min(sstart, send) - 1   # convert to 0-based
            s_hi = max(sstart, send)

            prev = best_hits.get(qname)
            if prev is None or bitscore > prev[2]:
                best_hits[qname] = (s_lo, s_hi, bitscore)

    mapped = [name for name in rRNA_NAMES if name in best_hits]
    if not mapped:
        logger.warning("No subunit aligned to consensus")
        return []

    # Build the feature list: 18S, ITS1, 5.8S, ITS2, 25S
    result: List[RdnaSubFeature] = []

    for name in rRNA_NAMES:
        if name in best_hits:
            s, e, _ = best_hits[name]
            result.append(RdnaSubFeature(name=name, start=s, end=e))

    # Sort by start position so we can derive ITS regions
    result.sort(key=lambda f: f.start)

    # Insert ITS1 and ITS2 between adjacent subunits
    # Expected order on the consensus: 18S ... 5.8S ... 25S
    # (or reverse, depending on consensus phase)
    rrna_by_name = {f.name: f for f in result}
    its_features: List[RdnaSubFeature] = []

    if "18S" in rrna_by_name and "5.8S" in rrna_by_name:
        f18s = rrna_by_name["18S"]
        f58s = rrna_by_name["5.8S"]
        if f18s.end <= f58s.start:
            its_features.append(RdnaSubFeature(
                name="ITS1", start=f18s.end, end=f58s.start,
            ))
        elif f58s.end <= f18s.start:
            # Reverse order on consensus
            its_features.append(RdnaSubFeature(
                name="ITS1", start=f58s.end, end=f18s.start,
            ))

    if "5.8S" in rrna_by_name and "25S" in rrna_by_name:
        f58s = rrna_by_name["5.8S"]
        f25s = rrna_by_name["25S"]
        if f58s.end <= f25s.start:
            its_features.append(RdnaSubFeature(
                name="ITS2", start=f58s.end, end=f25s.start,
            ))
        elif f25s.end <= f58s.start:
            its_features.append(RdnaSubFeature(
                name="ITS2", start=f25s.end, end=f58s.start,
            ))

    result.extend(its_features)
    result.sort(key=lambda f: f.start)

    for f in result:
        logger.info(f"  {f.name:5s}: {f.start}-{f.end} ({f.end - f.start} bp)")
    logger.info(f"Mapped {len(result)} sub-features to consensus")
    return result


def _run_blastn_subunits(
    query_fasta: Path,
    db_path: Path,
    output_path: Path,
    threads: int,
    err_path: Optional[Path] = None,
) -> None:
    """Run blastn (not megablast) for cross-species rRNA subunit search.

    Uses word_size=11 and evalue=1e-5 for sensitivity on short, conserved
    sequences like the 5.8S rRNA (~164 bp).
    """
    import subprocess

    if output_path.exists():
        logger.info(f"BLAST subunit output exists, reusing: {output_path}")
        return

    cmd = [
        "blastn",
        "-num_threads", str(threads),
        "-task", "blastn",
        "-db", str(db_path),
        "-query", str(query_fasta),
        "-evalue", "1e-5",
        "-max_hsps", "5",
        "-max_target_seqs", "1",
        "-outfmt", "6 std staxids",
        "-out", str(output_path),
    ]

    logger.info(f"Running blastn subunit annotation -> {output_path}")

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)
        with err_path.open("wb") as err_fh:
            ret = subprocess.call(cmd, stderr=err_fh)
    else:
        ret = subprocess.call(cmd, stderr=subprocess.DEVNULL)

    if ret != 0:
        raise RuntimeError(f"blastn subunit annotation failed with return code {ret}")


# ---------------------------------------------------------------------------
# Phase 5: Re-annotation of all contigs with the consensus probe
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
) -> List[RdnaLocus]:
    """BLAST the consensus 45S against all contigs to annotate rDNA loci.

    Args:
        query_fasta: Query assembly FASTA (all contigs)
        query_lengths: Contig lengths
        consensus: RdnaConsensus object with sequence and sub-features
        work_dir: Working directory
        threads: Number of threads
        classifications: Optional dict of contig_name -> classification string
        min_tandem_copies: Minimum tandem copies for NOR candidate [3]
        max_tandem_gap: Maximum gap between tandem copies [50000]

    Returns:
        List of RdnaLocus annotations
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
        return []

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
                sub_features=locus_dict["sub_features"],
                is_nor_candidate=False,  # Set below
            ))

    # Detect NOR candidates: contigs with >= min_tandem_copies in a cluster
    _mark_nor_candidates(loci, min_tandem_copies, max_tandem_gap)

    logger.info(f"Annotated {len(loci)} rDNA loci across {len(hits_per_contig)} contigs")
    nor_count = sum(1 for l in loci if l.is_nor_candidate)
    if nor_count:
        nor_contigs = set(l.contig for l in loci if l.is_nor_candidate)
        logger.info(f"NOR candidates: {len(nor_contigs)} contigs with {nor_count} loci")

    return loci


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

        # Determine which sub-features are present
        present_features = []
        for feat in sub_features:
            # Check if any hit overlaps this sub-feature on the consensus
            for h in group:
                if h["qstart"] <= feat.end and h["qend"] >= feat.start:
                    # Compute overlap fraction
                    overlap_start = max(h["qstart"], feat.start)
                    overlap_end = min(h["qend"], feat.end)
                    feat_len = feat.end - feat.start
                    if feat_len > 0 and (overlap_end - overlap_start) / feat_len >= 0.5:
                        present_features.append(feat.name)
                        break

        loci.append({
            "start": start,
            "end": end,
            "strand": strand,
            "identity": identity,
            "consensus_coverage": consensus_coverage,
            "copy_type": copy_type,
            "sub_features": present_features,
        })

    return loci


def _mark_nor_candidates(
    loci: List[RdnaLocus],
    min_tandem_copies: int,
    max_tandem_gap: int,
) -> None:
    """Mark loci that are part of tandem clusters as NOR candidates.

    Modifies loci in place.
    """
    # Group by contig
    by_contig: Dict[str, List[RdnaLocus]] = defaultdict(list)
    for locus in loci:
        by_contig[locus.contig].append(locus)

    for contig, contig_loci in by_contig.items():
        # Sort by start position
        contig_loci.sort(key=lambda l: l.start)

        # Find tandem clusters
        cluster: List[RdnaLocus] = [contig_loci[0]]
        for locus in contig_loci[1:]:
            if locus.start - cluster[-1].end <= max_tandem_gap:
                cluster.append(locus)
            else:
                if len(cluster) >= min_tandem_copies:
                    for l in cluster:
                        l.is_nor_candidate = True
                cluster = [locus]
        # Check last cluster
        if len(cluster) >= min_tandem_copies:
            for l in cluster:
                l.is_nor_candidate = True


# ---------------------------------------------------------------------------
# Phase 5b: Reclassify contigs using consensus-based coverage
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
) -> Tuple[Optional[RdnaConsensus], List[RdnaLocus]]:
    """Build consensus 45S rDNA and annotate all contigs.

    This is the main entry point that orchestrates Phases 2-5.

    Args:
        query_fasta: Query assembly FASTA
        query_lengths: Contig lengths
        rdna_hit_intervals: Per-contig BLAST hit intervals from Phase 1
        seed_ref_path: Seed rDNA reference FASTA (Arabidopsis or user-provided)
        work_dir: Working directory for all intermediate files
        threads: Number of threads
        rdna_ref_features_path: Optional path to sub-feature TSV for seed reference
        classifications: Optional dict of contig_name -> classification
        min_tandem_copies: Minimum tandem copies for NOR candidate [3]

    Returns:
        Tuple of (RdnaConsensus or None, list of RdnaLocus annotations)
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    logger.phase("rDNA consensus: Phase 2 - Extract rDNA regions and copies")

    if not rdna_hit_intervals:
        logger.warning("No rDNA hit intervals available; skipping consensus building")
        return None, []

    # Phase 2a: Extract rDNA-containing regions
    regions_fasta = work_dir / "rdna_regions.fa"
    region_map = _extract_rdna_regions(
        query_fasta=query_fasta,
        rdna_hit_intervals=rdna_hit_intervals,
        query_lengths=query_lengths,
        output_fasta=regions_fasta,
    )

    if not region_map:
        logger.warning("No rDNA regions extracted; skipping consensus building")
        return None, []

    # Phase 2b: Detect repeat period via self-alignment
    repeat_period = _detect_repeat_period(
        regions_fasta=regions_fasta,
        work_dir=work_dir / "self_align",
        threads=threads,
    )

    if not repeat_period:
        logger.warning("Could not determine repeat period; skipping consensus building")
        return None, []

    # Phase 2c: Extract individual copies
    copies_fasta = work_dir / "rdna_copies.fa"
    n_copies = _extract_individual_copies(
        regions_fasta=regions_fasta,
        region_map=region_map,
        repeat_period=repeat_period,
        output_fasta=copies_fasta,
    )

    if n_copies == 0:
        logger.warning("No rDNA copies extracted; skipping consensus building")
        return None, []

    # Phase 3: Cluster and select exemplar/consensus
    logger.phase("rDNA consensus: Phase 3 - Cluster copies and build consensus")
    consensus_seq, method, n_clustered = cluster_and_select_exemplar(
        copies_fasta=copies_fasta,
        work_dir=work_dir / "clustering",
        threads=threads,
    )

    if not consensus_seq:
        logger.warning("Failed to build consensus; skipping rDNA annotation")
        return None, []

    logger.info(f"Consensus: {len(consensus_seq)} bp, method={method}, "
                f"{n_copies} copies extracted, {n_clustered} in cluster")

    # Phase 4: Sub-feature annotation
    logger.phase("rDNA consensus: Phase 4 - Annotate sub-features")
    sub_features: List[RdnaSubFeature] = []

    # Load seed features
    if rdna_ref_features_path and rdna_ref_features_path.exists():
        seed_features = _load_features_tsv(rdna_ref_features_path)
    else:
        # Try bundled default
        script_dir = Path(__file__).resolve().parent.parent
        default_features = _find_default_features_tsv(script_dir)
        if default_features:
            seed_features = _load_features_tsv(default_features)
            logger.info(f"Using default features: {default_features}")
        else:
            seed_features = []
            logger.warning("No sub-feature annotation available")

    if seed_features:
        sub_features = annotate_sub_features(
            consensus_seq=consensus_seq,
            seed_ref_path=seed_ref_path,
            seed_features=seed_features,
            work_dir=work_dir / "annotation",
            threads=threads,
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

    # Phase 5: Re-annotate all contigs
    logger.phase("rDNA consensus: Phase 5 - Annotate all contigs")
    loci = annotate_contigs_with_consensus(
        query_fasta=query_fasta,
        query_lengths=query_lengths,
        consensus=consensus,
        work_dir=work_dir / "reannotation",
        threads=threads,
        classifications=classifications,
        min_tandem_copies=min_tandem_copies,
    )

    return consensus, loci
