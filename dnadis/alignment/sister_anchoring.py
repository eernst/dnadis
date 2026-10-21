"""Sister-assembly anchoring: a second assignment pass in multi-assembly runs.

Contigs with few reference anchors (TE-dense sequence, regions the reference
lacks or has diverged from) are aligned to the confidently assigned contigs of
the other assemblies in the run ("sisters").  Each sister votes for the
reference chromosome its matching region belongs to; the votes decide the
contig's assignment, subject to these rules:

* A sister hit must be a primary alignment of at least ``min_aln`` bp at
  gap-compressed identity of at least ``min_identity``.  Gap-compressed
  identity (minimap2 ``de`` tag) is used because accessions differ by
  insertions (e.g. TE insertion polymorphisms) that dominate BLAST-style
  identity even when the shared sequence is nearly identical.
* A hit takes the label of the anchor's *major reference segment* at the hit
  position, not the anchor's whole-contig assignment.  Anchors that carry a
  translocation or misjoin would otherwise propagate their overall label to
  sequence from another chromosome.  Small off-target blocks are ignored, as
  they are as poorly anchored as the contigs being placed.
* A sister whose hits cover two chromosomes at ``min_cov`` or more each marks
  the contig as chimeric; no call is made.
* At least ``min_votes`` sisters must agree, they must make up at least
  ``min_agree`` of the voting sisters, and the most similar sister (highest
  identity on this contig) must be among them.  A single sister is not
  enough: one assembly can carry sequence private to its own lineage.
* When the contig has reference evidence (``min_ref_bp`` of chained alignment
  to a chromosome), the sister call must be one of those chromosomes, compared
  by chromosome number so the subgenome can change but the chromosome cannot.
  A call outside that set keeps the reference result and is reported as a
  conflict.  This guards against structural variants shared by the sisters
  but absent from the assembly being placed.

The pass reads each assembly's first-pass outputs (contig summary, evidence
summary, macro blocks) and returns per-assembly calls; assemblies with calls
that change an assignment are rerun with those calls applied (see
``cli.run_assembly(sister_calls=...)``).
"""
from __future__ import annotations

import csv
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Set, Tuple

from dnadis.alignment.external_tools import get_minimap2_exe
from dnadis.utils.io_utils import (
    merge_intervals,
    open_maybe_gzip,
    output_is_current,
    replace_if_changed,
    run_to_gzip,
)
from dnadis.utils.logging_config import get_logger
from dnadis.utils.reference_utils import is_nuclear_chromosome, split_chrom_subgenome

logger = get_logger(__name__)

# Separator between assembly name and contig name in the pooled candidate FASTA.
_SEP = "|"

# Contig classes eligible as candidates.  Organelles, rDNA, cobionts, circular
# elements and chromosome debris have their own evidence and are not revisited.
CANDIDATE_CLASSES = frozenset(
    {"chrom_assigned", "chrom_fragment", "chrom_unassigned", "debris", "unclassified"}
)
# Contig classes eligible as anchors (with enough reference anchoring).
ANCHOR_CLASSES = frozenset({"chrom_assigned", "chrom_fragment"})

# Decisions that change the contig's assignment.
APPLIED_DECISIONS = frozenset({"assigned", "reassigned"})

# A reference chromosome is one of an anchor's major segments when its blocks
# make up at least this fraction of the anchor's anchored bp (and _MAJOR_MIN_BP).
_MAJOR_FRAC = 0.10
_MAJOR_MIN_BP = 100_000


@dataclass
class SisterVote:
    """One sister assembly's vote for one candidate contig."""
    label: str  # reference ID, or "CHIMERIC"
    coverage: float  # fraction of the candidate covered by hits with this label
    identity: float  # aligned-bp-weighted gap-compressed identity of all hits
    reversed: Optional[bool]  # candidate orientation relative to the reference


@dataclass
class SisterCall:
    """Second-pass decision for one candidate contig."""
    assembly: str
    contig: str  # original contig name
    contig_len: int
    first_class: str
    first_ref_id: str  # first-pass assigned reference ("" if none)
    ref_candidates: Set[str]  # reference IDs with >= min_ref_bp chained evidence
    decision: str  # assigned|reassigned|confirmed|conflict_keep_ref|no_call_*|no_vote
    ref_id: str = ""  # sister majority label ("" if no votes)
    n_agree: int = 0
    n_votes: int = 0
    identity: Optional[float] = None  # median identity of agreeing sisters
    coverage: Optional[float] = None  # median coverage of agreeing sisters
    reversed: Optional[bool] = None
    closest_sister: str = ""
    votes: Dict[str, SisterVote] = field(default_factory=dict)

    @property
    def applies(self) -> bool:
        return self.decision in APPLIED_DECISIONS

    @property
    def votes_str(self) -> str:
        return f"{self.n_agree}/{self.n_votes}" if self.n_votes else ""


# ---------------------------------------------------------------------------
# First-pass outputs
# ---------------------------------------------------------------------------
def _ref_frac(row: dict) -> float:
    length = int(row["length"] or 0)
    return int(row["total_aligned_bp"] or 0) / length if length > 0 else 0.0


def read_contig_summary(path: Path) -> List[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def select_candidates(
    rows: Iterable[dict], min_len: int, max_ref_frac: float,
) -> List[dict]:
    """Contigs of at least *min_len* bp with less than *max_ref_frac* of their
    length aligned to the reference, in a revisitable class."""
    return [
        r for r in rows
        if r["classification"] in CANDIDATE_CLASSES
        and int(r["length"] or 0) >= min_len
        and _ref_frac(r) < max_ref_frac
    ]


def select_anchors(rows: Iterable[dict], max_ref_frac: float) -> List[str]:
    """Chromosome contigs with at least *max_ref_frac* aligned to the reference."""
    return [
        r["original_name"] for r in rows
        if r["classification"] in ANCHOR_CLASSES and _ref_frac(r) >= max_ref_frac
    ]


def load_ref_candidates(evidence_tsv: Path, min_bp: int) -> Dict[str, Set[str]]:
    """contig -> nuclear reference IDs with at least *min_bp* chained evidence."""
    out: Dict[str, Set[str]] = defaultdict(set)
    if not evidence_tsv.exists():
        return out
    with open(evidence_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if is_nuclear_chromosome(r["ref_id"]) and int(r["union_bp"] or 0) >= min_bp:
                out[r["contig"]].add(r["ref_id"])
    return out


def load_anchor_segments(
    macro_blocks_tsv: Path, contigs: Optional[Set[str]] = None,
) -> Dict[str, List[Tuple[int, int, str, str]]]:
    """contig -> sorted (qstart, qend, ref_id, strand) for its major references.

    A reference counts as major when its blocks sum to at least
    ``_MAJOR_FRAC`` of the contig's anchored bp and ``_MAJOR_MIN_BP``.
    """
    blocks: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
    if not macro_blocks_tsv.exists():
        return blocks
    with open(macro_blocks_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = r["contig"]
            if contigs is not None and c not in contigs:
                continue
            if not is_nuclear_chromosome(r["ref_id"]):
                continue
            blocks[c].append((int(r["qstart"]), int(r["qend"]), r["ref_id"], r["strand"]))
    for c, v in blocks.items():
        per_ref: Counter = Counter()
        for s, e, ref, _strand in v:
            per_ref[ref] += e - s
        total = sum(per_ref.values())
        keep = {
            ref for ref, bp in per_ref.items()
            if bp >= max(_MAJOR_MIN_BP, _MAJOR_FRAC * total)
        }
        blocks[c] = sorted(b for b in v if b[2] in keep)
    return blocks


def segment_label(
    segments: List[Tuple[int, int, str, str]], ts: int, te: int,
) -> Optional[Tuple[str, str]]:
    """(ref_id, strand) of the major segment overlapping [ts, te), else the nearest."""
    if not segments:
        return None
    overlap: Counter = Counter()
    strand_bp: Dict[str, Counter] = defaultdict(Counter)
    for s, e, ref, strand in segments:
        o = min(e, te) - max(s, ts)
        if o > 0:
            overlap[ref] += o
            strand_bp[ref][strand] += o
    if overlap:
        ref = overlap.most_common(1)[0][0]
        return ref, strand_bp[ref].most_common(1)[0][0]
    nearest = min(
        segments, key=lambda b: (b[0] - te) if b[0] >= te else (ts - b[1]),
    )
    return nearest[2], nearest[3]


# ---------------------------------------------------------------------------
# Inputs and alignment
# ---------------------------------------------------------------------------
def _write_fasta_subset(
    src: Path, out: Path, names: Set[str], prefix: str = "",
) -> int:
    """Append records of *src* named in *names* to the open-for-append *out*.

    Returns the number of records written.
    """
    n = 0
    with open_maybe_gzip(src, "rt") as fh, out.open("a") as dst:
        writing = False
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                writing = name in names
                if writing:
                    dst.write(f">{prefix}{name}\n")
                    n += 1
            elif writing:
                dst.write(line)
    return n


def write_inputs(
    assemblies: Dict[str, Path],
    candidates: Dict[str, Set[str]],
    anchors: Dict[str, Set[str]],
    work_dir: Path,
) -> Tuple[Optional[Path], Dict[str, Path]]:
    """Write the pooled candidate FASTA and one anchor FASTA per assembly.

    Returns (candidates_fasta or None if no candidates, {assembly: anchors_fasta}).
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    cand_fa = work_dir / "candidates.fa"
    tmp = cand_fa.with_name(cand_fa.name + ".tmp")
    tmp.unlink(missing_ok=True)
    tmp.touch()
    n_cand = 0
    for name in sorted(candidates):
        if candidates[name]:
            n_cand += _write_fasta_subset(
                assemblies[name], tmp, candidates[name], prefix=f"{name}{_SEP}",
            )
    replace_if_changed(tmp, cand_fa)

    anchor_fas: Dict[str, Path] = {}
    anchor_dir = work_dir / "anchors"
    anchor_dir.mkdir(exist_ok=True)
    for name in sorted(anchors):
        if not anchors[name]:
            continue
        fa = anchor_dir / f"{name}.fa"
        tmp = fa.with_name(fa.name + ".tmp")
        tmp.unlink(missing_ok=True)
        tmp.touch()
        if _write_fasta_subset(assemblies[name], tmp, anchors[name]):
            replace_if_changed(tmp, fa)
            anchor_fas[name] = fa
        else:
            tmp.unlink()
    return (cand_fa if n_cand else None), anchor_fas


def run_sister_alignment(
    anchors_fasta: Path, candidates_fasta: Path, paf_gz: Path, threads: int,
) -> Path:
    """Align pooled candidates to one sister's anchors (minimap2 asm10, CIGAR).

    Reuses *paf_gz* when it is newer than both inputs.  Plain-typed arguments
    so the call can be submitted to a SLURM executor.
    """
    if output_is_current(paf_gz, anchors_fasta, candidates_fasta):
        logger.info(f"Reusing cached sister alignment: {paf_gz}")
        return paf_gz
    mapper = get_minimap2_exe()
    if not mapper:
        raise RuntimeError("Neither mm2plus nor minimap2 found in PATH")
    cmd = [
        mapper, "-c", "-x", "asm10", "--secondary=no", "-t", str(threads),
        str(anchors_fasta), str(candidates_fasta),
    ]
    err = paf_gz.with_name(paf_gz.name.replace(".paf.gz", ".err"))
    run_to_gzip(cmd, paf_gz, err)
    return paf_gz


# ---------------------------------------------------------------------------
# Votes and decisions
# ---------------------------------------------------------------------------
def _gap_compressed_identity(tokens: List[str]) -> Optional[float]:
    for t in tokens[12:]:
        if t.startswith("de:f:"):
            return 1.0 - float(t[5:])
    return None


def parse_sister_votes(
    paf_gz: Path,
    sister: str,
    anchor_segments: Dict[str, List[Tuple[int, int, str, str]]],
    min_aln: int,
    min_identity: float,
    min_cov: float,
) -> Dict[str, SisterVote]:
    """candidate ("asm|contig") -> this sister's vote."""
    hits: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    strand_bp: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    qlen: Dict[str, int] = {}
    id_sum: Dict[str, List[float]] = defaultdict(lambda: [0.0, 0])
    with open_maybe_gzip(paf_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or "tp:A:P" not in f[12:]:
                continue
            q = f[0]
            if q.split(_SEP, 1)[0] == sister:
                continue
            blk = int(f[10])
            if blk < min_aln:
                continue
            ident = _gap_compressed_identity(f)
            if ident is None or ident < min_identity:
                continue
            lab = segment_label(anchor_segments.get(f[5], []), int(f[7]), int(f[8]))
            if lab is None:
                continue
            ref, anchor_strand = lab
            qs, qe = int(f[2]), int(f[3])
            hits[q][ref].append((qs, qe))
            # Candidate is forward relative to the reference when its strand
            # relative to the anchor matches the anchor's strand relative to it.
            fwd = (f[4] == "+") == (anchor_strand == "+")
            strand_bp[q][ref]["+" if fwd else "-"] += qe - qs
            qlen[q] = int(f[1])
            id_sum[q][0] += ident * blk
            id_sum[q][1] += blk

    votes: Dict[str, SisterVote] = {}
    for q, per_ref in hits.items():
        cov = sorted(
            ((merge_intervals(iv)[1] / qlen[q], ref) for ref, iv in per_ref.items()),
            reverse=True,
        )
        if cov[0][0] < min_cov:
            continue
        identity = id_sum[q][0] / id_sum[q][1]
        if len(cov) > 1 and cov[1][0] >= min_cov:
            votes[q] = SisterVote("CHIMERIC", cov[0][0], identity, None)
            continue
        ref = cov[0][1]
        sb = strand_bp[q][ref]
        rev = None if sb["+"] == sb["-"] else sb["-"] > sb["+"]
        votes[q] = SisterVote(ref, cov[0][0], identity, rev)
    return votes


def _chrom(ref_id: str) -> str:
    return split_chrom_subgenome(ref_id)[0]


def decide(
    votes: Dict[str, SisterVote],
    first_ref_id: str,
    ref_candidates: Set[str],
    min_agree: float,
    min_votes: int = 2,
) -> Tuple[str, str, int, Optional[str]]:
    """Return (decision, sister_ref_id, n_agree, closest_sister)."""
    if not votes:
        return "no_vote", "", 0, None
    counts = Counter(v.label for v in votes.values())
    top, n = counts.most_common(1)[0]
    closest = max(votes, key=lambda s: (votes[s].identity, votes[s].coverage))
    if top == "CHIMERIC":
        return "no_call_chimeric", "", n, closest
    if n / len(votes) < min_agree:
        return "no_call_disagree", top, n, closest
    if n < min_votes:
        return "no_call_few_votes", top, n, closest
    if votes[closest].label != top:
        return "no_call_closest", top, n, closest
    if ref_candidates and _chrom(top) not in {_chrom(r) for r in ref_candidates}:
        return "conflict_keep_ref", top, n, closest
    if top == first_ref_id:
        return "confirmed", top, n, closest
    return ("reassigned" if first_ref_id else "assigned"), top, n, closest


def build_call(
    assembly: str,
    row: dict,
    votes: Dict[str, SisterVote],
    ref_candidates: Set[str],
    min_agree: float,
    min_votes: int = 2,
) -> SisterCall:
    first_ref = row["assigned_ref_id"] if row["assigned_ref_id"] not in ("", "NA") else ""
    decision, ref_id, n_agree, closest = decide(
        votes, first_ref, ref_candidates, min_agree, min_votes,
    )
    call = SisterCall(
        assembly=assembly,
        contig=row["original_name"],
        contig_len=int(row["length"]),
        first_class=row["classification"],
        first_ref_id=first_ref,
        ref_candidates=set(ref_candidates),
        decision=decision,
        ref_id=ref_id,
        n_agree=n_agree,
        n_votes=len(votes),
        closest_sister=closest or "",
        votes=votes,
    )
    agreeing = [v for v in votes.values() if v.label == ref_id] if ref_id else []
    if agreeing:
        call.identity = median(v.identity for v in agreeing)
        call.coverage = median(v.coverage for v in agreeing)
        orient = Counter(v.reversed for v in agreeing if v.reversed is not None)
        if orient and (len(orient) == 1 or orient[True] != orient[False]):
            call.reversed = orient.most_common(1)[0][0]
    return call


def write_calls_tsv(calls: Iterable[SisterCall], path: Path) -> None:
    header = [
        "assembly", "contig", "length", "first_classification", "first_ref_id",
        "ref_candidates", "decision", "sister_ref_id", "sister_votes",
        "closest_sister", "sister_identity", "sister_coverage", "reversed",
        "per_sister",
    ]
    with path.open("w") as out:
        out.write("\t".join(header) + "\n")
        for c in calls:
            per = ",".join(
                f"{s}:{v.label}@{v.identity:.4f}/{v.coverage:.2f}"
                for s, v in sorted(c.votes.items())
            )
            out.write("\t".join([
                c.assembly, c.contig, str(c.contig_len), c.first_class,
                c.first_ref_id or "NA", ",".join(sorted(c.ref_candidates)),
                c.decision, c.ref_id or "NA", c.votes_str, c.closest_sister,
                f"{c.identity:.4f}" if c.identity is not None else "",
                f"{c.coverage:.3f}" if c.coverage is not None else "",
                "" if c.reversed is None else ("yes" if c.reversed else "no"),
                per,
            ]) + "\n")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def run_sister_anchoring(
    results,
    work_dir: Path,
    *,
    min_len: int,
    max_ref_frac: float,
    min_identity: float,
    min_aln: int,
    min_cov: float,
    min_agree: float,
    min_votes: int,
    min_ref_bp: int,
    threads: int,
    executor=None,
    cluster_config=None,
) -> Dict[str, Dict[str, SisterCall]]:
    """Run the second pass over first-pass *results* (``AssemblyResult`` list).

    Returns {assembly_name: {original_contig_name: SisterCall}} for every
    candidate, and writes ``sister_anchoring.tsv`` to *work_dir*.
    """
    assemblies = {r.assembly_name: Path(r.assembly_path) for r in results}
    rows_by_asm: Dict[str, List[dict]] = {}
    cand_rows: Dict[str, Dict[str, dict]] = {}
    anchors: Dict[str, Set[str]] = {}
    for r in results:
        summary = Path(str(r.outprefix) + ".contig_summary.tsv")
        if not summary.exists():
            logger.warning(f"Sister anchoring: no contig summary for {r.assembly_name}; skipped")
            continue
        rows = read_contig_summary(summary)
        rows_by_asm[r.assembly_name] = rows
        cand_rows[r.assembly_name] = {
            row["original_name"]: row
            for row in select_candidates(rows, min_len, max_ref_frac)
        }
        anchors[r.assembly_name] = set(select_anchors(rows, max_ref_frac))

    n_cand = sum(len(v) for v in cand_rows.values())
    n_anchor_asm = sum(1 for v in anchors.values() if v)
    logger.info(
        f"Sister anchoring: {n_cand} candidate contigs (>= {min_len:,} bp, "
        f"< {max_ref_frac:.0%} reference-aligned); {n_anchor_asm} assemblies with anchors"
    )
    calls: Dict[str, Dict[str, SisterCall]] = defaultdict(dict)
    if n_cand == 0 or n_anchor_asm < 2:
        return calls

    cand_fa, anchor_fas = write_inputs(
        assemblies, {a: set(v) for a, v in cand_rows.items()}, anchors, work_dir,
    )
    if cand_fa is None:
        return calls

    use_cluster = bool(cluster_config and cluster_config.enabled)
    futures = {}
    paf_dir = work_dir / "paf"
    paf_dir.mkdir(exist_ok=True)
    for sister, fa in anchor_fas.items():
        paf_gz = paf_dir / f"{sister}.paf.gz"
        if executor is None:
            futures[sister] = run_sister_alignment(fa, cand_fa, paf_gz, threads)
            continue
        spec = None
        n_threads = threads
        if use_cluster:
            from dnadis.utils.resource_estimation import estimate_sister_alignment_resources
            spec = estimate_sister_alignment_resources(fa, cluster_config)
            n_threads = spec.cores
        futures[sister] = executor.submit(
            run_sister_alignment, fa, cand_fa, paf_gz, n_threads, resource_spec=spec,
        )

    votes: Dict[str, Dict[str, SisterVote]] = defaultdict(dict)
    for sister, fut in futures.items():
        try:
            paf_gz = fut if isinstance(fut, Path) else fut.result()
        except Exception as e:
            logger.warning(f"Sister anchoring: alignment to {sister} failed: {e}")
            continue
        segs = load_anchor_segments(
            Path(str(next(r.outprefix for r in results if r.assembly_name == sister))
                 + ".macro_blocks.tsv"),
            anchors[sister],
        )
        for q, v in parse_sister_votes(
            paf_gz, sister, segs, min_aln, min_identity, min_cov,
        ).items():
            votes[q][sister] = v

    for r in results:
        a = r.assembly_name
        if not cand_rows.get(a):
            continue
        ref_cands = load_ref_candidates(
            Path(str(r.outprefix) + ".evidence_summary.tsv"), min_ref_bp,
        )
        for contig, row in cand_rows[a].items():
            calls[a][contig] = build_call(
                a, row, votes.get(f"{a}{_SEP}{contig}", {}),
                ref_cands.get(contig, set()), min_agree, min_votes,
            )

    all_calls = [c for a in sorted(calls) for c in calls[a].values()]
    tsv = work_dir / "sister_anchoring.tsv"
    write_calls_tsv(all_calls, tsv)
    tally = Counter(c.decision for c in all_calls)
    logger.info(
        "Sister anchoring decisions: "
        + ", ".join(f"{k}={v}" for k, v in sorted(tally.items()))
    )
    logger.done(f"Sister anchoring: {tsv}")
    return calls


def apply_sister_calls(
    best_ref: Dict[str, str],
    calls: Dict[str, SisterCall],
    ref_orig_to_norm: Optional[Dict[str, str]] = None,
) -> Set[str]:
    """Set ``best_ref`` from the calls that change an assignment.

    Returns the contigs whose assignment came from the sister pass.
    """
    applied: Set[str] = set()
    for contig, call in calls.items():
        if call.applies and call.ref_id:
            best_ref[contig] = (ref_orig_to_norm or {}).get(call.ref_id, call.ref_id)
            applied.add(contig)
    return applied
