#!/usr/bin/env python3
"""
Organelle (chrC/chrM) detection for dnadis.

Contains functions for identifying chloroplast and mitochondrial contigs
using BLAST-based alignment to reference organelle sequences.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

from dnadis.detection.blast import run_makeblastdb, run_blastn_megablast
from dnadis.models import OrganelleCalls, OrganelleHit
from dnadis.utils.io_utils import file_exists_and_valid, merge_intervals
from dnadis.utils.logging_config import get_logger
from dnadis.utils.reference_utils import normalize_organelle_id
from dnadis.utils.sequence_utils import (
    write_filtered_fasta,
    is_hifiasm_circular,
    read_fasta_sequences,
    write_fasta,
)

logger = get_logger("organelle")


def extract_organelle_from_ref(
    ref_fasta: Path,
    organelle_id: str,
    output_path: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> bool:
    """Extract a specific sequence (chrC or chrM) from reference FASTA.

    Returns True if found and extracted, False otherwise.
    """
    if file_exists_and_valid(output_path):
        logger.info(f"Reusing cached {organelle_id} reference: {output_path}")
        return True

    sequences = read_fasta_sequences(ref_fasta)
    if ref_norm_to_orig:
        candidate = ref_norm_to_orig.get(organelle_id)
        if candidate and candidate in sequences:
            write_fasta({candidate: sequences[candidate]}, output_path)
            return True

    for name, seq in sequences.items():
        if normalize_organelle_id(name) == organelle_id:
            write_fasta({name: seq}, output_path)
            return True
    return False


def prepare_organelle_references(
    ref_fasta: Path,
    chrC_ref_arg: Optional[str],
    chrM_ref_arg: Optional[str],
    work_dir: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Prepare organelle reference FASTAs.

    Priority:
    1. User-provided path (if it exists as a file)
    2. chrC/chrM from reference FASTA if present
    3. None (skip this organelle)

    Returns (chrC_fasta_path, chrM_fasta_path)
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    chrC_path: Optional[Path] = None
    chrM_path: Optional[Path] = None

    # Handle chloroplast reference
    if chrC_ref_arg:
        chrC_arg_path = Path(chrC_ref_arg)
        if chrC_arg_path.exists():
            chrC_path = chrC_arg_path
            logger.info(f"Using user-provided chrC reference: {chrC_path}")
        else:
            logger.warning(f"chrC reference not found: {chrC_ref_arg}")
    else:
        # Try to extract from reference
        chrC_extracted = work_dir / "chrC_ref.fa"
        if extract_organelle_from_ref(ref_fasta, "chrC", chrC_extracted, ref_norm_to_orig):
            chrC_path = chrC_extracted
            logger.info(f"Extracted chrC from reference: {chrC_path}")
        else:
            logger.info("No chrC found in reference FASTA")

    # Handle mitochondrial reference
    if chrM_ref_arg:
        chrM_arg_path = Path(chrM_ref_arg)
        if chrM_arg_path.exists():
            chrM_path = chrM_arg_path
            logger.info(f"Using user-provided chrM reference: {chrM_path}")
        else:
            logger.warning(f"chrM reference not found: {chrM_ref_arg}")
    else:
        # Try to extract from reference
        chrM_extracted = work_dir / "chrM_ref.fa"
        if extract_organelle_from_ref(ref_fasta, "chrM", chrM_extracted, ref_norm_to_orig):
            chrM_path = chrM_extracted
            logger.info(f"Extracted chrM from reference: {chrM_path}")
        else:
            logger.info("No chrM found in reference FASTA")

    return chrC_path, chrM_path


@dataclass(frozen=True)
class OrganelleThresholds:
    """Criteria for organelle calls (see ``call_organelles``).

    A contig is ``complete`` when it covers ``min_ref_cov`` of the reference,
    at least ``min_query_cov`` of the contig aligns, and its length is within
    the organelle's (min, max) ratio to the reference length.  When no contig
    passes, the fewest contigs (each at ``min_query_cov``) that together cover
    ``set_min_ref_cov`` of the reference are ``partial``.  Other contigs with
    at least ``debris_min_query_cov`` aligned are ``debris``.
    """
    min_ref_cov: float = 0.95
    min_query_cov: float = 0.90
    set_min_ref_cov: float = 0.90
    debris_min_query_cov: float = 0.50
    same_genome_identity: float = 0.99
    len_ratio: Tuple[Tuple[str, float, float], ...] = (
        ("chrC", 0.90, 1.50),
        ("chrM", 0.50, 3.00),
    )

    def len_bounds(self, organelle: str) -> Tuple[float, float]:
        for name, lo, hi in self.len_ratio:
            if name == organelle:
                return lo, hi
        return 0.0, float("inf")


# A set member must add at least this fraction of the reference.
_SET_MIN_GAIN = 0.01

# Two contigs belong to the same organelle genome at or above
# OrganelleThresholds.same_genome_identity with at least this fraction of the
# query aligned.  Copies of one plastome in an assembly are near-identical;
# plastomes of related species are typically 94-99% identical.
_SAME_GENOME_MIN_COV = 0.50

# Without read depth, the genome with the most organelle sequence is used only
# when it has at least this many times the sequence of the next one.
_GROUP_DOMINANCE = 2.0

# (organelle, complete candidates, all contigs assigned to that organelle)
#   -> (candidate -> genome group, group -> bp of assigned contigs in it)
GenomeGrouper = Callable[[str, List[str], List[str]], Tuple[Dict[str, int], Dict[int, int]]]


@dataclass
class _PairEvidence:
    query_intervals: List[Tuple[int, int]] = field(default_factory=list)
    ref_intervals: List[Tuple[int, int]] = field(default_factory=list)
    matches: float = 0.0
    aln_len: int = 0
    # (query start, query end, identity, alignment length) per HSP
    hsps: List[Tuple[int, int, float, int]] = field(default_factory=list)


def parse_organelle_blast(blast_out: Path) -> Dict[Tuple[str, str], _PairEvidence]:
    """Collect query/reference intervals and matches per (contig, organelle)."""
    pairs: Dict[Tuple[str, str], _PairEvidence] = defaultdict(_PairEvidence)
    if not (blast_out.exists() and blast_out.stat().st_size > 0):
        return pairs
    with blast_out.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                pident = float(fields[2])
                aln_length = int(fields[3])
                qstart, qend = sorted((int(fields[6]), int(fields[7])))
                sstart, send = sorted((int(fields[8]), int(fields[9])))
            except ValueError:
                continue
            ev = pairs[(fields[0], fields[1])]
            # BLAST coordinates are 1-based inclusive; store half-open intervals
            ev.query_intervals.append((qstart - 1, qend))
            ev.ref_intervals.append((sstart - 1, send))
            ev.matches += pident * aln_length / 100.0
            ev.aln_len += aln_length
            ev.hsps.append((qstart - 1, qend, pident / 100.0, aln_length))
    return pairs


def _covered_bp(intervals: List[Tuple[int, int]]) -> int:
    return merge_intervals(intervals)[1] if intervals else 0


def call_organelles(
    pairs: Dict[Tuple[str, str], _PairEvidence],
    query_lengths: Dict[str, int],
    ref_lengths: Dict[str, int],
    thresholds: OrganelleThresholds = OrganelleThresholds(),
    circular_contigs: Optional[Set[str]] = None,
    genome_groups: Optional[GenomeGrouper] = None,
) -> OrganelleCalls:
    """Call complete, partial and debris organelle contigs from BLAST evidence.

    Each contig is assigned to the organelle it matches best (aligned fraction
    of the contig x identity): plastid-derived stretches of the mitogenome
    otherwise label chloroplast fragments as mitochondrial.

    Complete candidates are ranked circular first, then by length closest to
    the reference.  A candidate that is not selected stays organelle debris.

    A sample can hold the organelle genome of another species (a second
    plant in the culture, contamination), and its contigs can pass the
    complete tests.  With ``genome_groups``, candidates are grouped into
    genomes by their identity to each other, and only the genome holding the
    most organelle sequence in the assembly is ranked; the others are
    reported in ``OrganelleCalls.other_genomes``.
    """
    circular = circular_contigs if circular_contigs is not None else {
        c for c, _ in pairs if is_hifiasm_circular(c)
    }

    # Best organelle per contig
    best: Dict[str, Tuple[str, OrganelleHit, _PairEvidence]] = {}
    for (contig, organelle), ev in pairs.items():
        ref_len = ref_lengths.get(organelle, 0)
        qlen = query_lengths.get(contig, 0)
        if ref_len <= 0 or qlen <= 0:
            continue
        hit = OrganelleHit(
            organelle_type=organelle,
            coverage=_covered_bp(ev.query_intervals) / qlen,
            identity=(ev.matches / ev.aln_len) if ev.aln_len else 0.0,
            length_ratio=qlen / ref_len,
            is_complete=False,
            ref_coverage=_covered_bp(ev.ref_intervals) / ref_len,
        )
        key = (hit.coverage * hit.identity, hit.ref_coverage)
        prev = best.get(contig)
        if prev is None or key > (prev[1].coverage * prev[1].identity, prev[1].ref_coverage):
            best[contig] = (organelle, hit, ev)

    calls = OrganelleCalls()
    for organelle in sorted(ref_lengths):
        ref_len = ref_lengths[organelle]
        members = {c: (h, ev) for c, (o, h, ev) in best.items() if o == organelle}
        lo, hi = thresholds.len_bounds(organelle)

        candidates = [
            c for c, (h, _) in members.items()
            if h.ref_coverage >= thresholds.min_ref_cov
            and h.coverage >= thresholds.min_query_cov
            and lo <= h.length_ratio <= hi
        ]
        candidates.sort(key=lambda c: (
            c not in circular,
            abs(query_lengths[c] - ref_len),
            -members[c][0].ref_coverage,
            c,
        ))

        if genome_groups is not None and len(candidates) > 1:
            group_of, group_bp = genome_groups(
                organelle, candidates,
                sorted(c for c, (h, _) in members.items()
                       if h.coverage >= thresholds.debris_min_query_cov),
            )
            groups = {group_of[c] for c in candidates}
            if len(groups) > 1:
                calls.genome_groups[organelle] = dict(group_of)
                calls.complete_candidates[organelle] = list(candidates)
                first_rank = {g: min(i for i, c in enumerate(candidates) if group_of[c] == g)
                              for g in groups}
                by_bp = sorted(groups, key=lambda g: (-group_bp.get(g, 0), first_rank[g]))
                main, runner_up = by_bp[0], by_bp[1]
                if group_bp.get(main, 0) >= _GROUP_DOMINANCE * group_bp.get(runner_up, 0):
                    others = [c for c in candidates if group_of[c] != main]
                    candidates = [c for c in candidates if group_of[c] == main]
                    calls.other_genomes[organelle] = others
                    logger.warning(
                        f"{organelle}: candidates come from {len(groups)} different genomes "
                        f"(< {thresholds.same_genome_identity:.1%} identity to each other). Using "
                        f"the one with the most {organelle} sequence in the assembly "
                        f"({group_bp.get(main, 0):,} vs {group_bp.get(runner_up, 0):,} bp); "
                        f"kept as organelle_debris: " + ", ".join(others)
                    )
                else:
                    logger.warning(
                        f"{organelle}: candidates come from {len(groups)} different genomes "
                        f"(< {thresholds.same_genome_identity:.1%} identity to each other) with "
                        f"similar amounts of sequence ({group_bp.get(main, 0):,} vs "
                        f"{group_bp.get(runner_up, 0):,} bp); ranking all of them. With read "
                        f"depth, the genome with the most read bases is used."
                    )

        chosen: Optional[str] = candidates[0] if candidates else None
        partial: List[str] = []
        if chosen is None:
            partial, set_cov = _greedy_organelle_set(
                {c: ev for c, (h, ev) in members.items()
                 if h.coverage >= thresholds.min_query_cov},
                ref_len, thresholds.set_min_ref_cov, circular, query_lengths,
            )
            if partial:
                calls.partial[organelle] = partial
                calls.set_ref_coverage[organelle] = set_cov
        else:
            calls.set_ref_coverage[organelle] = members[chosen][0].ref_coverage

        for contig, (hit, _) in members.items():
            if contig == chosen:
                hit.status, hit.is_complete = "complete", True
            elif contig in partial:
                hit.status = "partial"
            elif hit.coverage >= thresholds.debris_min_query_cov:
                calls.debris.add(contig)
            else:
                continue
            calls.hits[contig] = hit

        if organelle == "chrC":
            calls.chrC_contig = chosen
        elif organelle == "chrM":
            calls.chrM_contig = chosen

        if chosen:
            h = members[chosen][0]
            logger.info(
                f"Selected {organelle}: {chosen} ({query_lengths[chosen]:,} bp, "
                f"{h.length_ratio:.2f}x reference, ref cov={h.ref_coverage:.3f}, "
                f"query cov={h.coverage:.3f}, ident={h.identity:.4f})"
            )
        elif partial:
            logger.info(
                f"No single {organelle} contig passes; {len(partial)} contigs together "
                f"cover {calls.set_ref_coverage[organelle]:.3f} of the reference: "
                + ", ".join(partial)
            )
        elif members:
            top = max(members, key=lambda c: members[c][0].ref_coverage)
            h = members[top][0]
            logger.info(
                f"No {organelle} call; best contig {top} covers {h.ref_coverage:.3f} "
                f"of the reference ({h.length_ratio:.2f}x reference length, "
                f"query cov={h.coverage:.3f})"
            )
    return calls


def _greedy_organelle_set(
    eligible: Dict[str, _PairEvidence],
    ref_len: int,
    min_ref_cov: float,
    circular: Set[str],
    query_lengths: Dict[str, int],
) -> Tuple[List[str], float]:
    """Fewest contigs that together cover ``min_ref_cov`` of the reference.

    Greedy set cover: each step adds the contig that covers the most reference
    bp not yet covered (ties: circular, then shorter).  Returns ([], coverage)
    when the threshold is not reached.
    """
    chosen: List[str] = []
    covered: List[Tuple[int, int]] = []
    covered_bp = 0
    remaining = dict(eligible)
    while remaining and covered_bp < min_ref_cov * ref_len:
        gains = {
            c: _covered_bp(covered + ev.ref_intervals) - covered_bp
            for c, ev in remaining.items()
        }
        pick = min(gains, key=lambda c: (-gains[c], c not in circular, query_lengths.get(c, 0), c))
        if gains[pick] < _SET_MIN_GAIN * ref_len:
            break
        chosen.append(pick)
        covered = merge_intervals(covered + remaining.pop(pick).ref_intervals)[0]
        covered_bp = _covered_bp(covered)
    frac = covered_bp / ref_len if ref_len else 0.0
    return (chosen, frac) if frac >= min_ref_cov else ([], frac)


def detect_organelles(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    chrC_ref: Optional[Path],
    chrM_ref: Optional[Path],
    work_dir: Path,
    threads: int,
    thresholds: OrganelleThresholds = OrganelleThresholds(),
    circular_contigs: Optional[Set[str]] = None,
) -> OrganelleCalls:
    """BLAST contigs against chrC/chrM and call organelle contigs.

    See ``call_organelles`` for the criteria.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    if not chrC_ref and not chrM_ref:
        logger.info("No organelle references available, skipping organelle detection")
        return OrganelleCalls()

    # Rename sequences to "chrC" and "chrM" so BLAST subject IDs match expected names
    combined_seqs: Dict[str, str] = {}
    for name, ref in (("chrC", chrC_ref), ("chrM", chrM_ref)):
        if ref:
            seqs = read_fasta_sequences(ref)
            if seqs:
                combined_seqs[name] = next(iter(seqs.values()))
    ref_lengths = {name: len(seq) for name, seq in combined_seqs.items()}
    combined_ref = work_dir / "organelle_refs.fa"
    write_fasta(combined_seqs, combined_ref)

    organelle_db = work_dir / "organelle_refs"
    run_makeblastdb(combined_ref, organelle_db, err_path=work_dir / "makeblastdb.err")

    # Organelle genomes (chloroplast inverted repeats, mitochondrial repeats)
    # produce many HSPs per query-subject pair.
    blast_out = work_dir / "organelle_blast.txt"
    run_blastn_megablast(
        query_fasta=query_fasta,
        db_paths=[str(organelle_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=1000,
        err_path=work_dir / "blastn.err",
    )

    return call_organelles(
        parse_organelle_blast(blast_out),
        query_lengths,
        ref_lengths,
        thresholds=thresholds,
        circular_contigs=circular_contigs,
        genome_groups=blast_genome_grouper(
            query_fasta, query_lengths, work_dir, threads,
            min_identity=thresholds.same_genome_identity,
        ),
    )


def blast_genome_grouper(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    work_dir: Path,
    threads: int,
    min_identity: float = 0.99,
) -> GenomeGrouper:
    """Group complete candidates into genomes with megablast (see call_organelles).

    Candidates at >= ``min_identity`` over half their length are linked
    (single linkage).  Every contig assigned to the organelle counts toward the group
    of the candidate it matches best, if at that identity.
    """
    def grouper(organelle: str, candidates: List[str], assigned: List[str]):
        gdir = work_dir / f"genome_groups_{organelle}"
        gdir.mkdir(parents=True, exist_ok=True)
        cand_fa, query_fa = gdir / "candidates.fa", gdir / "assigned.fa"
        write_filtered_fasta(query_fasta, cand_fa, set(candidates))
        write_filtered_fasta(query_fasta, query_fa, set(assigned) | set(candidates))
        db = gdir / "candidates"
        run_makeblastdb(cand_fa, db, err_path=gdir / "makeblastdb.err")
        out = gdir / "assigned_vs_candidates.txt"
        run_blastn_megablast(
            query_fasta=query_fa, db_paths=[str(db)], output_path=out, threads=threads,
            max_hsps=100, max_target_seqs=len(candidates), err_path=gdir / "blastn.err",
        )
        return group_by_identity(parse_organelle_blast(out), candidates, assigned,
                                 query_lengths, min_identity)

    return grouper


def group_by_identity(
    pairs: Dict[Tuple[str, str], _PairEvidence],
    candidates: List[str],
    assigned: List[str],
    query_lengths: Dict[str, int],
    min_identity: float = 0.99,
) -> Tuple[Dict[str, int], Dict[int, int]]:
    """Genome groups and the bp of assigned contigs in each.

    Returns a group for every candidate and every assigned contig that
    matches one, and the summed contig length per group.
    """
    def same_genome(q: str, s: str) -> Optional[float]:
        ev = pairs.get((q, s))
        qlen = query_lengths.get(q, 0)
        if ev is None or not ev.hsps or qlen <= 0:
            return None
        ident, covered = _best_hsp_identity(ev)
        return ident if ident >= min_identity and covered / qlen >= _SAME_GENOME_MIN_COV else None

    parent = {c: c for c in candidates}

    def find(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for q in candidates:
        for s in candidates:
            if q != s and same_genome(q, s) is not None:
                parent[find(q)] = find(s)
    roots = sorted({find(c) for c in candidates})
    group_of = {c: roots.index(find(c)) for c in candidates}

    group_bp: Dict[int, int] = defaultdict(int)
    for q in sorted(set(assigned) | set(candidates)):
        if q not in group_of:
            hits = [(ident, s) for s in candidates
                    if (ident := same_genome(q, s)) is not None]
            if not hits:
                continue
            group_of[q] = group_of[max(hits)[1]]
        group_bp[group_of[q]] += query_lengths.get(q, 0)
    return group_of, dict(group_bp)


def _best_hsp_identity(ev: _PairEvidence) -> Tuple[float, int]:
    """Identity and aligned bp over non-overlapping HSPs, longest first.

    Pooling every HSP would count matches between different repeat copies
    (frequent in mitogenomes), which lowers the identity of two copies of
    one genome.
    """
    taken: List[Tuple[int, int]] = []
    matches = 0.0
    aligned = 0
    for qs, qe, ident, length in sorted(ev.hsps, key=lambda h: (-h[3], -h[2])):
        if any(qs < e and s < qe for s, e in taken):
            continue
        taken.append((qs, qe))
        matches += ident * length
        aligned += length
    return (matches / aligned if aligned else 0.0), aligned


def repick_by_read_depth(
    calls: OrganelleCalls,
    depth_median: Dict[str, float],
    query_lengths: Dict[str, int],
) -> Dict[str, Tuple[Optional[str], str]]:
    """Re-choose complete organelle contigs by read depth when candidates
    come from more than one genome.

    The genome whose contigs hold the most read bases (median depth x length,
    summed over every contig assigned to it) is the sample's own; another
    species' organelle genome in the sample is at lower depth even where it
    has as much or more assembled sequence.  Updates ``calls`` in place and
    returns organelle -> (previous complete contig, new one) for changes.
    """
    swaps: Dict[str, Tuple[Optional[str], str]] = {}
    for organelle, candidates in calls.complete_candidates.items():
        group_of = calls.genome_groups.get(organelle, {})
        read_bases: Dict[int, float] = defaultdict(float)
        for contig, group in group_of.items():
            depth = depth_median.get(contig)
            if depth:
                read_bases[group] += depth * query_lengths.get(contig, 0)
        if not read_bases:
            continue
        main = max(read_bases, key=lambda g: read_bases[g])
        logger.info(
            f"{organelle}: read bases per genome: "
            + ", ".join(
                f"{'/'.join(c for c in candidates if group_of[c] == g)[:60]} "
                f"{read_bases.get(g, 0) / 1e6:,.1f} Mb"
                for g in sorted({group_of[c] for c in candidates})
            )
        )
        new = next((c for c in candidates if group_of[c] == main), None)
        if new is None:
            continue
        calls.other_genomes[organelle] = [c for c in candidates if group_of[c] != main]
        current = getattr(calls, f"{organelle}_contig")
        if new == current:
            continue
        if current is not None:
            old_hit = calls.hits[current]
            old_hit.status, old_hit.is_complete = "debris", False
            calls.debris.add(current)
        new_hit = calls.hits[new]
        new_hit.status, new_hit.is_complete = "complete", True
        calls.debris.discard(new)
        setattr(calls, f"{organelle}_contig", new)
        calls.set_ref_coverage[organelle] = new_hit.ref_coverage
        swaps[organelle] = (current, new)
        logger.warning(
            f"{organelle}: {new} replaces {current} as the complete contig; its genome "
            f"has the most read bases ({read_bases[main] / 1e6:,.1f} Mb)"
        )
    return swaps
