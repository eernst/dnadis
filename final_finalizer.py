#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ----------------------------
# Small utilities
# ----------------------------
def open_maybe_gzip(path: Path, mode: str):
    """
    Open plain or gzipped text/binary files based on suffix.
    mode should be compatible with open()/gzip.open().
    """
    if path.suffix == ".gz" or path.suffix == ".bgz":
        return gzip.open(path, mode)
    return open(path, mode)


def have_exe(exe: str) -> bool:
    """Check if an executable exists in PATH using shutil.which()."""
    return shutil.which(exe) is not None


def normalize_ref_id(ref_id: str) -> str:
    """Normalize reference IDs to a lowercase 'chr' prefix when present.

    Handles various case patterns: Chr1 -> chr1, CHR1 -> chr1, etc.
    IDs already starting with lowercase 'chr' are unchanged.
    """
    if len(ref_id) >= 3 and ref_id[:3].upper() == "CHR" and not ref_id.startswith("chr"):
        return "chr" + ref_id[3:]
    return ref_id


def build_ref_id_maps(ref_ids) -> tuple[dict[str, str], dict[str, str]]:
    """Build bidirectional mappings between original and normalized reference IDs.

    Args:
        ref_ids: Iterable of original reference IDs from the reference FASTA.

    Returns:
        Tuple of (orig_to_norm, norm_to_orig) dictionaries.
        - orig_to_norm: Maps each original ID to its normalized form
        - norm_to_orig: Maps each normalized ID back to an original ID

    Note:
        If multiple original IDs normalize to the same value (e.g., "ChrM" and "Mt"
        both normalize to "chrM"), the first encountered ID is used for the
        reverse mapping (norm_to_orig). This preserves the reference's primary
        naming convention in output files.
    """
    orig_to_norm: dict[str, str] = {}
    norm_to_orig: dict[str, str] = {}
    for ref_id in ref_ids:
        norm = normalize_ref_id(ref_id)
        organelle = normalize_organelle_id(ref_id)
        if organelle:
            norm = organelle
        orig_to_norm[ref_id] = norm
        if norm not in norm_to_orig:
            norm_to_orig[norm] = ref_id
    return orig_to_norm, norm_to_orig


def normalize_ref_lengths(ref_lengths: Dict[str, int], orig_to_norm: Dict[str, str]) -> Dict[str, int]:
    normalized: Dict[str, int] = {}
    for ref_id, length in ref_lengths.items():
        norm_id = orig_to_norm.get(ref_id, normalize_ref_id(ref_id))
        if norm_id not in normalized:
            normalized[norm_id] = int(length)
        else:
            normalized[norm_id] = max(int(length), normalized[norm_id])
    return normalized


def read_fasta_lengths_with_map(
    fasta_path: Path,
) -> tuple[Dict[str, int], Dict[str, str], Dict[str, str]]:
    lengths = read_fasta_lengths(fasta_path)
    orig_to_norm, norm_to_orig = build_ref_id_maps(lengths.keys())
    normalized = normalize_ref_lengths(lengths, orig_to_norm)
    return normalized, orig_to_norm, norm_to_orig


def run_to_gzip(cmd: list[str], gz_out: Path, err_path: Path) -> None:
    """
    Run a command, streaming stdout into a gzipped file. Stderr -> err_path.
    """
    err_path.parent.mkdir(parents=True, exist_ok=True)

    with err_path.open("wb") as err_fh, gzip.open(gz_out, "wb") as gz_fh:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err_fh)
        try:
            assert proc.stdout is not None
            for chunk in iter(lambda: proc.stdout.read(1 << 20), b""):
                gz_fh.write(chunk)
        finally:
            try:
                if proc.stdout:
                    proc.stdout.close()
            finally:
                ret = proc.wait()

    if ret != 0:
        raise RuntimeError(f"Command failed with return code {ret} (see {err_path})")


# ----------------------------
# Minimap2/mm2plus detection and helpers
# ----------------------------
# mm2plus is a drop-in replacement for minimap2 with additional features.
# We prefer mm2plus if available, otherwise fall back to minimap2.
_MINIMAP2_EXE: Optional[str] = None


def get_minimap2_exe() -> Optional[str]:
    """Detect and cache the minimap2/mm2plus executable.

    Prefers mm2plus if available, falls back to minimap2.
    Returns None if neither is found.
    """
    global _MINIMAP2_EXE
    if _MINIMAP2_EXE is None:
        if have_exe("mm2plus"):
            _MINIMAP2_EXE = "mm2plus"
        elif have_exe("minimap2"):
            _MINIMAP2_EXE = "minimap2"
        else:
            _MINIMAP2_EXE = ""  # Mark as checked but not found
    return _MINIMAP2_EXE if _MINIMAP2_EXE else None


def run_minimap2(
    ref: Path,
    qry: Path,
    paf_out: Path,
    threads: int,
    preset: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    gzip_output: bool = False,
    err_path: Optional[Path] = None,
) -> bool:
    """Run minimap2/mm2plus alignment.

    Uses mm2plus if available, otherwise minimap2. Both are functionally
    equivalent for our use cases.

    Args:
        ref: Reference FASTA path
        qry: Query FASTA path
        paf_out: Output PAF path (will be gzipped if gzip_output=True)
        threads: Number of threads
        preset: Minimap2 preset (e.g., "asm5", "asm20")
        extra_args: Additional command-line arguments
        gzip_output: If True, compress output with gzip
        err_path: Path for stderr output (optional)

    Returns:
        True if alignment succeeded, False otherwise
    """
    # Check if output already exists
    if paf_out.exists():
        print(f"[info] PAF exists, reusing: {paf_out}", file=sys.stderr)
        return True

    mapper = get_minimap2_exe()
    if not mapper:
        print("[warn] Neither mm2plus nor minimap2 found in PATH", file=sys.stderr)
        return False

    # Build command
    cmd = [mapper, "-t", str(threads)]
    if preset:
        cmd += ["-x", preset]
    if extra_args:
        cmd.extend(extra_args)

    if gzip_output:
        # Output to stdout for gzip compression
        cmd += [str(ref), str(qry)]
    else:
        # Output directly to file
        cmd += ["-o", str(paf_out), str(ref), str(qry)]

    print(f"[info] Running {mapper} -> {paf_out}", file=sys.stderr)

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if gzip_output:
            run_to_gzip(cmd, paf_out, err_path or Path("/dev/null"))
        else:
            if err_path:
                with err_path.open("wb") as err_fh:
                    ret = subprocess.run(cmd, stderr=err_fh, check=False)
            else:
                ret = subprocess.run(cmd, stderr=subprocess.DEVNULL, check=False)
            if ret.returncode != 0:
                print(f"[warn] {mapper} failed with return code {ret.returncode}", file=sys.stderr)
                return False
    except Exception as e:
        print(f"[warn] {mapper} failed: {e}", file=sys.stderr)
        return False

    return True


# ----------------------------
# External tool helpers
# ----------------------------
def run_minimap2_synteny(
    ref: Path,
    qry: Path,
    paf_gz_out: Path,
    threads: int,
    preset: Optional[str],
    k: Optional[int],
    w: Optional[int],
    err_path: Path,
) -> None:
    """Run minimap2/mm2plus with parameterization tuned for cross-species synteny alignments.

    This uses specialized parameters (-c -f1 -B2 -O1,24 -s 10 --max-chain-skip=80 -z 2000,200)
    for genome-wide synteny analysis. Output is gzipped PAF.
    """
    if paf_gz_out.exists():
        print(f"[info] PAF.gz exists, reusing: {paf_gz_out}", file=sys.stderr)
        return

    mapper = get_minimap2_exe()
    if not mapper:
        raise RuntimeError("Neither mm2plus nor minimap2 found in PATH")

    cmd = [
        mapper,
        "-c",
        "-f1",
        "-B2",
        "-O1,24",
        "-s", "10",
        "--max-chain-skip=80",
        "-z", "2000,200",
        "-t", str(threads),
    ]
    if preset:
        cmd += ["-x", preset]
    if k is not None:
        cmd += ["-k", str(k)]
    if w is not None:
        cmd += ["-w", str(w)]
    cmd += [str(ref), str(qry)]

    print(
        f"[info] Running {mapper} -> {paf_gz_out} (stderr: {err_path}): " + " ".join(cmd),
        file=sys.stderr,
    )
    run_to_gzip(cmd, paf_gz_out, err_path)


def run_gffread_extract_proteins(
    gffread_exe: str,
    ref_fasta: Path,
    ref_gff3: Path,
    proteins_faa: Path,
    err_path: Path,
) -> None:
    """
    Use gffread to extract protein sequences from reference genome + GFF3.
      gffread -g ref.fa -y proteins.fa ref.gff3
    """
    if proteins_faa.exists():
        print(f"[info] Proteins FASTA exists, reusing: {proteins_faa}", file=sys.stderr)
        return

    if not have_exe(gffread_exe):
        raise RuntimeError(f"gffread executable not found/usable: {gffread_exe}")

    cmd = [gffread_exe, "-g", str(ref_fasta), "-y", str(proteins_faa), str(ref_gff3)]
    print(f"[info] Extracting proteins with gffread (stderr: {err_path}): " + " ".join(cmd), file=sys.stderr)

    err_path.parent.mkdir(parents=True, exist_ok=True)
    with err_path.open("wb") as err_fh:
        ret = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=err_fh)

    if ret != 0:
        raise RuntimeError(f"gffread failed with return code {ret} (see {err_path})")

    if not proteins_faa.exists() or proteins_faa.stat().st_size == 0:
        raise RuntimeError(f"gffread produced empty proteins FASTA: {proteins_faa}")


def run_miniprot(
    miniprot_exe: str,
    query_fasta: Path,
    proteins_faa: Path,
    paf_gz_out: Path,
    threads: int,
    extra_args: str,
    err_path: Path,
) -> None:
    """
    Align proteins to query genome with miniprot.
    Typical usage:
      miniprot -t N <target_genome.fa> <query_proteins.faa> > out.paf
    Here target = query genome, query = reference proteins.
    Output is gzipped PAF.
    """
    if paf_gz_out.exists():
        print(f"[info] miniprot PAF.gz exists, reusing: {paf_gz_out}", file=sys.stderr)
        return

    if not have_exe(miniprot_exe):
        raise RuntimeError(f"miniprot executable not found/usable: {miniprot_exe}")

    cmd = [miniprot_exe, "-t", str(threads)]
    if extra_args:
        cmd += extra_args.strip().split()
    cmd += [str(query_fasta), str(proteins_faa)]

    print(
        f"[info] Running miniprot -> {paf_gz_out} (stderr: {err_path}): " + " ".join(cmd),
        file=sys.stderr,
    )
    run_to_gzip(cmd, paf_gz_out, err_path)


# ----------------------------
# PAF helpers
# ----------------------------
def paf_tag_value(tokens: list[str], prefix: str) -> Optional[str]:
    for tok in tokens[12:]:
        if tok.startswith(prefix):
            # tok like "tp:A:P" -> last field "P"
            return tok.split(":")[-1]
    return None


def is_primary_only(tokens: list[str]) -> bool:
    """
    Keep only primary alignments (tp:A:P). Drop supplementary (tp:A:I) and secondary (tp:A:S).
    If tp tag is absent, conservatively keep.
    """
    tp = paf_tag_value(tokens, "tp:A:")
    if tp is None:
        return True
    return tp == "P"


# Nuclear chromosomes:
#   chr1            (single-set reference; no suffix)
#   chr1A, chr1P... (multi-set reference; suffix is any A-Z)
_NUCLEAR_SINGLE_RE = re.compile(r"^(?P<chrom>chr\d+)$")
_NUCLEAR_RE = re.compile(r"^(?P<chrom>chr\d+)(?P<sg>[A-Z])$")

# Organelles (plant convention)
_ORGANELLE = {"chrC", "chrM"}
_ORGANELLE_ALIASES = {
    "chrC": {"chrc", "cpdna", "chloroplast", "plastid", "pt"},
    "chrM": {"chrm", "mt", "mitochondrion", "mitochondria", "mitochondrial"},
}

# Reference naming presets (tried in order when pattern not provided)
_REF_ID_PATTERNS = [
    re.compile(r"^(?P<chrom>chr\d+)(?P<sg>[A-Z])$", re.IGNORECASE),
    re.compile(r"^(?P<chrom>chr\d+)(?P<sg>At|Dt)$", re.IGNORECASE),
    re.compile(r"^(?P<chrom>chr\d+)$", re.IGNORECASE),
]
_ACTIVE_REF_ID_PATTERNS: Optional[List[re.Pattern[str]]] = None


def normalize_organelle_id(ref_id: str) -> Optional[str]:
    if not ref_id:
        return None
    norm = normalize_ref_id(ref_id)
    if norm in _ORGANELLE:
        return norm
    norm_l = norm.lower()
    if norm_l in _ORGANELLE_ALIASES["chrC"]:
        return "chrC"
    if norm_l in _ORGANELLE_ALIASES["chrM"]:
        return "chrM"
    return None


def get_ref_id_patterns() -> List[re.Pattern[str]]:
    """Return active reference ID patterns (custom if set, else defaults)."""
    return _ACTIVE_REF_ID_PATTERNS if _ACTIVE_REF_ID_PATTERNS else _REF_ID_PATTERNS


def set_ref_id_patterns(patterns: Optional[List[re.Pattern[str]]]) -> None:
    """Set custom reference ID patterns. Pass None to reset to defaults.

    This is called once at startup from main() if --ref-id-pattern is provided.
    For testing, call with None to reset to default patterns.
    """
    global _ACTIVE_REF_ID_PATTERNS
    _ACTIVE_REF_ID_PATTERNS = patterns


def compile_ref_id_patterns(patterns: List[str]) -> List[re.Pattern[str]]:
    compiled: List[re.Pattern[str]] = []
    for pat in patterns:
        try:
            c = re.compile(pat, re.IGNORECASE)
        except re.error as exc:
            raise ValueError(f"Invalid --ref-id-pattern regex: {pat} ({exc})") from exc
        if "chrom" not in c.groupindex:
            raise ValueError(f"--ref-id-pattern must define a (?P<chrom>...) group: {pat}")
        compiled.append(c)
    return compiled


def get_min_nuclear_chrom_length(ref_lengths: Dict[str, int]) -> int:
    """
    Return the length of the smallest non-organelle reference chromosome.
    Returns 0 if no nuclear chromosomes found.
    """
    nuclear_lengths = [
        length for ref_id, length in ref_lengths.items()
        if normalize_organelle_id(ref_id) is None
    ]
    return min(nuclear_lengths) if nuclear_lengths else 0


def split_chrom_subgenome(
    ref_id: str,
    patterns: Optional[List[re.Pattern[str]]] = None,
) -> tuple[str, str]:
    """Split reference ID into chromosome and subgenome components.

    Args:
        ref_id: Reference sequence ID (e.g., "chr1A", "Chr5", "chrM")
        patterns: Optional list of compiled regex patterns to use. If None,
                  uses get_ref_id_patterns() (global/default patterns).

    Returns:
        Tuple of (chrom_id, subgenome). Subgenome is "NA" if not present.
        Preserves original case style (Chr vs chr) in output.
    """
    if not ref_id:
        return "", "NA"

    style = "Chr" if ref_id.startswith("Chr") else "chr" if ref_id.startswith("chr") else ""
    norm_id = normalize_ref_id(ref_id)

    organelle = normalize_organelle_id(norm_id)
    if organelle:
        chrom_id = organelle
        if style == "Chr" and chrom_id.startswith("chr"):
            chrom_id = "Chr" + chrom_id[3:]
        return chrom_id, "NA"

    active_patterns = patterns if patterns is not None else get_ref_id_patterns()
    for pat in active_patterns:
        m = pat.match(norm_id)
        if not m:
            continue
        chrom_id = m.groupdict().get("chrom", norm_id)
        sub = m.groupdict().get("sg", "NA")
        if style == "Chr" and chrom_id.startswith("chr"):
            chrom_id = "Chr" + chrom_id[3:]
        return chrom_id, sub

    return ref_id, "NA"


def read_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """
    Return dict: seqname -> length for a FASTA that may be gzipped.
    """
    lengths: dict[str, int] = {}
    name: Optional[str] = None
    seqlen = 0

    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seqlen
                name = line[1:].strip().split()[0]
                seqlen = 0
            else:
                seqlen += len(line.strip())

    if name is not None:
        lengths[name] = seqlen

    return lengths


def write_ref_lengths_tsv(ref_lengths_path: Path, ref_fasta: Path) -> None:
    lengths = read_fasta_lengths(ref_fasta)
    with ref_lengths_path.open("w") as out:
        out.write("ref_id\tchrom_id\tsubgenome\tref_len\n")
        for ref_id, L in sorted(lengths.items()):
            chrom_id, sub = split_chrom_subgenome(ref_id)
            out.write(f"{ref_id}\t{chrom_id}\t{sub}\t{L}\n")


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], int]:
    """
    Merge closed-open intervals [(s,e),...] where s<e.
    Returns merged list and total merged length.
    """
    if not intervals:
        return [], 0
    intervals = sorted(intervals)
    merged: list[tuple[int, int]] = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    total = sum(e - s for s, e in merged)
    return merged, total


# ----------------------------
# Legacy: primary-only stats (QA)
# ----------------------------
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


# ----------------------------
# GFF3 transcript coordinate mapping
# ----------------------------
def parse_gff3_transcript_coords(
    gff3_path: Path,
    ref_id_map: Optional[Dict[str, str]] = None,
):
    """
    Parse GFF3 to map transcript IDs (mRNA/transcript) -> (ref_chrom, start, end, strand)
    AND transcript IDs -> gene ID (from Parent=...; fallback to transcript ID).

    Returns:
      tx2loc:  dict[str] -> (chrom, start0, end0, strand)
      tx2gene: dict[str] -> gene_id
    """
    tx2loc: dict[str, tuple[str, int, int, str]] = {}
    tx2gene: dict[str, str] = {}

    def _attrs_to_dict(attr_str: str) -> dict[str, str]:
        d: dict[str, str] = {}
        for kv in attr_str.split(";"):
            kv = kv.strip()
            if not kv:
                continue
            if "=" in kv:
                k, v = kv.split("=", 1)
                d[k] = v
        return d

    with open_maybe_gzip(gff3_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = parts
            ftype_l = ftype.lower()
            if ftype_l not in ("mrna", "transcript"):
                continue

            try:
                s = int(start) - 1  # 0-based half-open
                e = int(end)
            except ValueError:
                continue

            ad = _attrs_to_dict(attrs)
            tid = ad.get("ID")
            if not tid:
                continue

            parent = ad.get("Parent", "")
            gene_id = parent.split(",")[0] if parent else tid

            ref_id = ref_id_map.get(seqid, normalize_ref_id(seqid)) if ref_id_map else seqid
            tx2loc[tid] = (ref_id, s, e, strand if strand in ("+", "-") else "+")
            tx2gene[tid] = gene_id

    return tx2loc, tx2gene


# ----------------------------
# Synteny-block chaining (shared with both modes)
# ----------------------------
@dataclass
class Block:
    qs: int
    qe: int
    rs: int
    re: int
    matches: int
    aln_len: int
    mapq: int
    strand: str
    gene_id: Optional[str] = None


@dataclass
class Chain:
    last_qe: int
    last_r: int
    diag0: int
    q_intervals: list[tuple[int, int]]
    matches_sum: int
    alnlen_sum: int
    gene_ids: Set[str]


@dataclass
class ChainEvidenceResult:
    """Result container for chain evidence and segment parsing functions."""
    qlens_from_paf: Dict[str, int]
    qr_union_bp: Dict[Tuple[str, str], int]
    qr_matches: Dict[Tuple[str, str], int]
    qr_alnlen: Dict[Tuple[str, str], int]
    qr_gene_count: Dict[Tuple[str, str], int]
    contig_total: Dict[str, int]
    contig_refs: Dict[str, Set[str]]
    qr_score_topk: Dict[Tuple[str, str], float]
    qr_weight_all: Dict[Tuple[str, str], float]
    qr_nchains_kept: Dict[Tuple[str, str], int]
    best_ref: Dict[str, str]
    best_score: Dict[str, float]
    best_bp: Dict[str, int]
    second_ref: Dict[str, str]
    second_score: Dict[str, float]
    second_bp: Dict[str, int]
    chain_segments_rows: List[Tuple]
    macro_block_rows: List[Tuple]
    chain_summary_rows: List[Tuple]


# ----------------------------
# Contig classification data structures
# ----------------------------
@dataclass
class ContigClassification:
    """Classification result for a single contig."""
    original_name: str
    new_name: str
    classification: str  # chrom, chrom_debris, organelle_complete, organelle_debris, rDNA, contaminant, unclassified
    reversed: bool
    contaminant_taxid: Optional[int]  # NCBI taxonomy ID for contaminants
    contaminant_sci: Optional[str]  # Scientific name for contaminants
    assigned_ref_id: Optional[str]
    ref_gene_proportion: Optional[float]
    contig_len: int


@dataclass
class BlastHitSummary:
    """Aggregated BLAST hits for a contig."""
    contig_name: str
    total_coverage: float
    best_hit_subject: str
    best_hit_taxid: Optional[int]
    best_hit_evalue: float
    total_aligned_bp: int


# ----------------------------
# FASTA sequence utilities
# ----------------------------
_COMPLEMENT = str.maketrans("ATCGatcgNnRrYySsWwKkMmBbDdHhVv", "TAGCtagcNnYyRrSsWwMmKkVvHhDdBb")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def read_fasta_sequences(fasta_path: Path) -> Dict[str, str]:
    """Read FASTA file and return dict of name -> sequence."""
    sequences: Dict[str, str] = {}
    name: Optional[str] = None
    seq_parts: List[str] = []

    with open_maybe_gzip(fasta_path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(seq_parts)
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if name is not None:
        sequences[name] = "".join(seq_parts)

    return sequences


def write_fasta(sequences: Dict[str, str], output_path: Path, wrap: int = 80) -> None:
    """Write sequences to FASTA file with line wrapping."""
    with output_path.open("w") as out:
        for name, seq in sequences.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i + wrap] + "\n")


def write_filtered_fasta(
    input_fasta: Path,
    output_fasta: Path,
    include_contigs: Set[str],
    wrap: int = 80,
) -> None:
    """Write a filtered FASTA containing only specified contigs.

    Streams through the input to avoid loading entire genome into memory.
    """
    with open_maybe_gzip(input_fasta, "rt") as fh, output_fasta.open("w") as out:
        writing = False
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                writing = name in include_contigs
                if writing:
                    out.write(line)
            elif writing:
                out.write(line)


def is_hifiasm_circular(contig_name: str) -> bool:
    """Check if contig name matches hifiasm circular pattern (ptg*c)."""
    return bool(re.match(r"ptg\d+c", contig_name))


# ----------------------------
# BLAST infrastructure
# ----------------------------
def run_makeblastdb(
    fasta_path: Path,
    db_path: Path,
    dbtype: str = "nucl",
    err_path: Optional[Path] = None,
) -> None:
    """Create BLAST database from FASTA file."""
    if not have_exe("makeblastdb"):
        raise RuntimeError("makeblastdb executable not found in PATH")

    # Check if database already exists
    db_files = list(db_path.parent.glob(f"{db_path.name}.n*"))
    if db_files:
        print(f"[info] BLAST database exists, reusing: {db_path}", file=sys.stderr)
        return

    db_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-out", str(db_path),
        "-dbtype", dbtype,
        "-hash_index",
    ]

    print(f"[info] Creating BLAST database: {db_path}", file=sys.stderr)

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)
        with err_path.open("wb") as err_fh:
            ret = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=err_fh)
    else:
        ret = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if ret != 0:
        raise RuntimeError(f"makeblastdb failed with return code {ret}")


def run_blastn_megablast(
    query_fasta: Path,
    db_paths: List[str],
    output_path: Path,
    threads: int,
    max_hsps: int = 500,
    max_target_seqs: int = 5,
    taxidlist: Optional[Path] = None,
    err_path: Optional[Path] = None,
) -> None:
    """Run blastn with megablast task.

    Output format: "6 std staxids" (standard 12 columns + taxonomy IDs)
    Supports gzipped query FASTA files by piping decompressed data via stdin.

    Args:
        query_fasta: Query FASTA file (may be gzipped)
        db_paths: List of BLAST database paths
        output_path: Output file path
        threads: Number of threads
        max_hsps: Maximum HSPs per subject
        max_target_seqs: Maximum target sequences per query
        taxidlist: Optional path to file containing taxids to restrict search
        err_path: Optional path for stderr output
    """
    if output_path.exists():
        print(f"[info] BLAST output exists, reusing: {output_path}", file=sys.stderr)
        return

    if not have_exe("blastn"):
        raise RuntimeError("blastn executable not found in PATH")

    # Join database paths with spaces
    db_str = " ".join(db_paths)

    # Check if query is gzipped
    is_gzipped = query_fasta.suffix in (".gz", ".bgz")

    # Build common command options
    base_cmd = [
        "blastn",
        "-num_threads", str(threads),
        "-task", "megablast",
        "-db", db_str,
        "-max_hsps", str(max_hsps),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", "6 std staxids",
        "-out", str(output_path),
    ]

    # Add taxid restriction if provided
    if taxidlist and taxidlist.exists():
        base_cmd += ["-taxidlist", str(taxidlist)]

    if is_gzipped:
        # Use stdin for gzipped files
        cmd = base_cmd + ["-query", "-"]
    else:
        cmd = base_cmd + ["-query", str(query_fasta)]

    print(f"[info] Running blastn megablast -> {output_path}", file=sys.stderr)

    if err_path:
        err_path.parent.mkdir(parents=True, exist_ok=True)

    if is_gzipped:
        # Use gunzip -c to decompress and pipe to BLAST
        gunzip_cmd = ["gunzip", "-c", str(query_fasta)]
        gunzip_proc = subprocess.Popen(gunzip_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        if err_path:
            with err_path.open("wb") as err_fh:
                ret = subprocess.call(cmd, stdin=gunzip_proc.stdout, stderr=err_fh)
        else:
            ret = subprocess.call(cmd, stdin=gunzip_proc.stdout, stderr=subprocess.DEVNULL)
        gunzip_proc.stdout.close()
        gunzip_proc.wait()
    else:
        if err_path:
            with err_path.open("wb") as err_fh:
                ret = subprocess.call(cmd, stderr=err_fh)
        else:
            ret = subprocess.call(cmd, stderr=subprocess.DEVNULL)

    if ret != 0:
        raise RuntimeError(f"blastn failed with return code {ret}")


def parse_blast_coverage(
    blast_path: Path,
    query_lengths: Dict[str, int],
) -> Dict[str, BlastHitSummary]:
    """Parse BLAST tabular output and compute per-query coverage.

    Expected format: "6 std staxids" (qseqid, sseqid, pident, length, mismatch,
    gapopen, qstart, qend, sstart, send, evalue, bitscore, staxids)

    Returns dict: contig_name -> BlastHitSummary
    """
    # Collect alignment intervals per query
    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    query_best_hit: Dict[str, Tuple[str, Optional[int], float]] = {}  # (subject, taxid, evalue)

    if not blast_path.exists() or blast_path.stat().st_size == 0:
        return {}

    with blast_path.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            qseqid = fields[0]
            sseqid = fields[1]
            try:
                qstart = int(fields[6])
                qend = int(fields[7])
                evalue = float(fields[10])
            except ValueError:
                continue

            # Parse taxid if present (column 13, may contain multiple semicolon-separated)
            taxid: Optional[int] = None
            if len(fields) >= 13 and fields[12]:
                taxid_str = fields[12].split(";")[0]
                try:
                    taxid = int(taxid_str)
                except ValueError:
                    taxid = None

            # Ensure qstart < qend
            if qstart > qend:
                qstart, qend = qend, qstart

            query_intervals[qseqid].append((qstart, qend))

            # Track best hit by evalue
            if qseqid not in query_best_hit or evalue < query_best_hit[qseqid][2]:
                query_best_hit[qseqid] = (sseqid, taxid, evalue)

    # Compute merged coverage for each query
    results: Dict[str, BlastHitSummary] = {}
    for qseqid, intervals in query_intervals.items():
        _, total_bp = merge_intervals(intervals)
        qlen = query_lengths.get(qseqid, 0)
        coverage = (total_bp / qlen) if qlen > 0 else 0.0

        best_subject, best_taxid, best_evalue = query_best_hit.get(qseqid, ("", None, float("inf")))

        results[qseqid] = BlastHitSummary(
            contig_name=qseqid,
            total_coverage=coverage,
            best_hit_subject=best_subject,
            best_hit_taxid=best_taxid,
            best_hit_evalue=best_evalue,
            total_aligned_bp=total_bp,
        )

    return results


# ----------------------------
# Organelle detection
# ----------------------------
def extract_organelle_from_ref(
    ref_fasta: Path,
    organelle_id: str,
    output_path: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> bool:
    """Extract a specific sequence (chrC or chrM) from reference FASTA.

    Returns True if found and extracted, False otherwise.
    """
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
            print(f"[info] Using user-provided chrC reference: {chrC_path}", file=sys.stderr)
        else:
            print(f"[warn] chrC reference not found: {chrC_ref_arg}", file=sys.stderr)
    else:
        # Try to extract from reference
        chrC_extracted = work_dir / "chrC_ref.fa"
        if extract_organelle_from_ref(ref_fasta, "chrC", chrC_extracted, ref_norm_to_orig):
            chrC_path = chrC_extracted
            print(f"[info] Extracted chrC from reference: {chrC_path}", file=sys.stderr)
        else:
            print("[info] No chrC found in reference FASTA", file=sys.stderr)

    # Handle mitochondrial reference
    if chrM_ref_arg:
        chrM_arg_path = Path(chrM_ref_arg)
        if chrM_arg_path.exists():
            chrM_path = chrM_arg_path
            print(f"[info] Using user-provided chrM reference: {chrM_path}", file=sys.stderr)
        else:
            print(f"[warn] chrM reference not found: {chrM_ref_arg}", file=sys.stderr)
    else:
        # Try to extract from reference
        chrM_extracted = work_dir / "chrM_ref.fa"
        if extract_organelle_from_ref(ref_fasta, "chrM", chrM_extracted, ref_norm_to_orig):
            chrM_path = chrM_extracted
            print(f"[info] Extracted chrM from reference: {chrM_path}", file=sys.stderr)
        else:
            print("[info] No chrM found in reference FASTA", file=sys.stderr)

    return chrC_path, chrM_path


def detect_organelles(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    chrC_ref: Optional[Path],
    chrM_ref: Optional[Path],
    work_dir: Path,
    threads: int,
    min_coverage: float,
    chrC_len_tol: float,
    chrM_len_tol: float,
) -> Tuple[Optional[str], Optional[str], Set[str]]:
    """Identify best organelle candidates.

    Criteria for full organelle:
    - BLAST coverage >= min_coverage
    - Length within tolerance of reference
    - Prefer hifiasm circular contigs (ptg*c pattern)

    Returns (chrC_contig, chrM_contig, debris_contigs)
    where debris_contigs are those with 50%+ coverage but not selected as main organelle.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    chrC_contig: Optional[str] = None
    chrM_contig: Optional[str] = None
    debris_contigs: Set[str] = set()

    # Get reference lengths
    chrC_len = 0
    chrM_len = 0
    if chrC_ref:
        chrC_lengths = read_fasta_lengths(chrC_ref)
        chrC_len = list(chrC_lengths.values())[0] if chrC_lengths else 0
    if chrM_ref:
        chrM_lengths = read_fasta_lengths(chrM_ref)
        chrM_len = list(chrM_lengths.values())[0] if chrM_lengths else 0

    # Combine organelle references for BLAST if available
    if not chrC_ref and not chrM_ref:
        print("[info] No organelle references available, skipping organelle detection", file=sys.stderr)
        return None, None, set()

    # Create combined organelle reference for BLAST
    # Important: rename sequences to "chrC" and "chrM" so BLAST subject IDs match expected names
    combined_ref = work_dir / "organelle_refs.fa"
    combined_seqs: Dict[str, str] = {}
    if chrC_ref:
        chrC_seqs = read_fasta_sequences(chrC_ref)
        if chrC_seqs:
            # Take the first (or only) sequence and rename to "chrC"
            first_seq = next(iter(chrC_seqs.values()))
            combined_seqs["chrC"] = first_seq
    if chrM_ref:
        chrM_seqs = read_fasta_sequences(chrM_ref)
        if chrM_seqs:
            # Take the first (or only) sequence and rename to "chrM"
            first_seq = next(iter(chrM_seqs.values()))
            combined_seqs["chrM"] = first_seq
    write_fasta(combined_seqs, combined_ref)

    # Create BLAST database
    organelle_db = work_dir / "organelle_refs"
    run_makeblastdb(combined_ref, organelle_db, err_path=work_dir / "makeblastdb.err")

    # Run BLAST
    # Note: organelle genomes (especially chloroplasts with inverted repeats) can produce
    # many HSPs per query-subject pair. Use higher value to capture all alignments.
    blast_out = work_dir / "organelle_blast.txt"
    run_blastn_megablast(
        query_fasta=query_fasta,
        db_paths=[str(organelle_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=1000,
        err_path=work_dir / "blastn.err",
    )

    # Parse BLAST results and compute coverage per (query, subject) pair
    query_subject_intervals: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)

    if blast_out.exists() and blast_out.stat().st_size > 0:
        with blast_out.open("r") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue

                qseqid = fields[0]
                sseqid = fields[1]
                try:
                    qstart = int(fields[6])
                    qend = int(fields[7])
                except ValueError:
                    continue

                if qstart > qend:
                    qstart, qend = qend, qstart

                query_subject_intervals[(qseqid, sseqid)].append((qstart, qend))

    # Compute coverage and identify candidates
    chrC_candidates: List[Tuple[str, float, int, bool]] = []  # (contig, coverage, length, is_circular)
    chrM_candidates: List[Tuple[str, float, int, bool]] = []

    for (qseqid, sseqid), intervals in query_subject_intervals.items():
        _, total_bp = merge_intervals(intervals)
        qlen = query_lengths.get(qseqid, 0)
        coverage = (total_bp / qlen) if qlen > 0 else 0.0

        is_circular = is_hifiasm_circular(qseqid)

        if sseqid == "chrC" and chrC_len > 0:
            len_diff = abs(qlen - chrC_len) / chrC_len
            if coverage >= min_coverage and len_diff <= chrC_len_tol:
                chrC_candidates.append((qseqid, coverage, qlen, is_circular))
            elif coverage >= 0.50:
                debris_contigs.add(qseqid)

        elif sseqid == "chrM" and chrM_len > 0:
            len_diff = abs(qlen - chrM_len) / chrM_len
            if coverage >= min_coverage and len_diff <= chrM_len_tol:
                chrM_candidates.append((qseqid, coverage, qlen, is_circular))
            elif coverage >= 0.50:
                debris_contigs.add(qseqid)

    # Select best candidates (prefer circular, then by length similarity)
    def select_best(candidates: List[Tuple[str, float, int, bool]], ref_len: int) -> Optional[str]:
        if not candidates:
            return None
        # Sort by: circular (desc), then length difference (asc)
        candidates.sort(key=lambda x: (-int(x[3]), abs(x[2] - ref_len)))
        return candidates[0][0]

    chrC_contig = select_best(chrC_candidates, chrC_len)
    chrM_contig = select_best(chrM_candidates, chrM_len)

    # Remove selected organelles from debris
    if chrC_contig:
        debris_contigs.discard(chrC_contig)
        print(f"[info] Selected chrC candidate: {chrC_contig}", file=sys.stderr)
    if chrM_contig:
        debris_contigs.discard(chrM_contig)
        print(f"[info] Selected chrM candidate: {chrM_contig}", file=sys.stderr)

    return chrC_contig, chrM_contig, debris_contigs


# ----------------------------
# rDNA detection
# ----------------------------
def prepare_rdna_reference(
    rdna_ref_arg: Optional[str],
    script_dir: Path,
) -> Optional[Path]:
    """Prepare rDNA reference FASTA.

    If rdna_ref_arg is None or 'default', use data/athal-45s-ref.fa.
    Otherwise, use the provided path.
    """
    if rdna_ref_arg is None or rdna_ref_arg.lower() == "default":
        default_path = script_dir / "data" / "athal-45s-ref.fa"
        if default_path.exists():
            print(f"[info] Using default rDNA reference: {default_path}", file=sys.stderr)
            return default_path
        else:
            print(f"[warn] Default rDNA reference not found: {default_path}", file=sys.stderr)
            return None

    rdna_path = Path(rdna_ref_arg)
    if rdna_path.exists():
        print(f"[info] Using user-provided rDNA reference: {rdna_path}", file=sys.stderr)
        return rdna_path

    print(f"[warn] rDNA reference not found: {rdna_ref_arg}", file=sys.stderr)
    return None


def detect_rdna_contigs(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    rdna_ref: Path,
    work_dir: Path,
    threads: int,
    min_coverage: float,
    exclude_contigs: Set[str],
) -> Set[str]:
    """Identify contigs with significant rDNA content.

    Returns set of contig names with coverage >= min_coverage.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Create BLAST database for rDNA
    rdna_db = work_dir / "rdna_ref"
    run_makeblastdb(rdna_ref, rdna_db, err_path=work_dir / "makeblastdb_rdna.err")

    # Run BLAST
    # rDNA arrays are highly repetitive; use higher max_hsps to capture full coverage
    blast_out = work_dir / "rdna_blast.txt"
    run_blastn_megablast(
        query_fasta=query_fasta,
        db_paths=[str(rdna_db)],
        output_path=blast_out,
        threads=threads,
        max_hsps=100,
        err_path=work_dir / "blastn_rdna.err",
    )

    # Parse results
    blast_results = parse_blast_coverage(blast_out, query_lengths)

    rdna_contigs: Set[str] = set()
    for contig, summary in blast_results.items():
        if contig in exclude_contigs:
            continue
        if summary.total_coverage >= min_coverage:
            rdna_contigs.add(contig)
            print(f"[info] rDNA contig: {contig} (coverage={summary.total_coverage:.2f})", file=sys.stderr)

    return rdna_contigs


# ----------------------------
# Chromosome debris detection
# ----------------------------
def detect_chromosome_debris(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    chromosome_contigs: Set[str],
    work_dir: Path,
    threads: int,
    min_coverage: float = 0.80,
    min_identity: float = 0.90,
    exclude_contigs: Set[str] = None,
) -> Set[str]:
    """Identify contigs that are assembly artifacts: near-identical copies of chromosome contigs.

    MOTIVATION: Genome assemblers (especially hifiasm) can produce multiple representations
    of the same genomic region - haplotigs, bubble branches, or other duplicates. These
    contigs are nearly identical to portions of the primary chromosome contigs and should
    be classified as debris rather than left as "unclassified" or mistakenly flagged as
    contaminants. This detection method aligns candidate contigs against the *assembled*
    chromosome contigs (not the reference) to catch these assembly-specific artifacts.

    This complements reference-based debris detection (classify_debris_and_unclassified),
    which catches sequences with homology to the reference genome but may miss:
    - Non-coding duplicates (no protein hits)
    - Assembly artifacts specific to this assembly (not in reference)

    Uses minimap2/mm2plus with asm5 preset (assembly-to-assembly, high identity) to find
    contigs with high coverage (>=80%) and high identity (>=90%) to chromosome contigs.

    Args:
        query_fasta: Path to query FASTA containing all contigs
        query_lengths: Dict of contig name -> length
        chromosome_contigs: Set of contigs already classified as chromosomes
        work_dir: Working directory for intermediate files
        threads: Number of threads for alignment
        min_coverage: Min fraction of query aligned to classify as debris (default 0.80)
        min_identity: Min alignment identity to classify as debris (default 0.90)
        exclude_contigs: Set of contigs to exclude from debris detection

    Returns:
        Set of contig names classified as chromosome_debris
    """
    if exclude_contigs is None:
        exclude_contigs = set()

    work_dir.mkdir(parents=True, exist_ok=True)

    if not chromosome_contigs:
        print("[info] No chromosome contigs for debris detection", file=sys.stderr)
        return set()

    # Check for minimap2/mm2plus
    if not get_minimap2_exe():
        print("[warn] Neither minimap2 nor mm2plus found, skipping chromosome debris detection", file=sys.stderr)
        return set()

    # Create FASTA of chromosome contigs as reference (streaming to avoid loading entire genome)
    chrs_fasta = work_dir / "chromosome_contigs.fa"
    if not chrs_fasta.exists():
        write_filtered_fasta(query_fasta, chrs_fasta, chromosome_contigs)
    print(f"[info] Created chromosome reference for debris detection: {len(chromosome_contigs)} contigs", file=sys.stderr)

    # Create FASTA of candidate contigs to test (streaming)
    excluded = chromosome_contigs | exclude_contigs
    candidate_contigs = set(query_lengths.keys()) - excluded

    if not candidate_contigs:
        print("[info] No candidate contigs for debris detection", file=sys.stderr)
        return set()

    candidates_fasta = work_dir / "debris_candidates.fa"
    if not candidates_fasta.exists():
        write_filtered_fasta(query_fasta, candidates_fasta, candidate_contigs)

    # Run minimap2/mm2plus: candidates (query) vs chromosomes (reference)
    paf_out = work_dir / "debris_vs_chrs.paf"
    success = run_minimap2(
        ref=chrs_fasta,
        qry=candidates_fasta,
        paf_out=paf_out,
        threads=threads,
        preset="asm5",  # assembly-to-assembly preset, high identity
        extra_args=["-c"],  # output CIGAR
    )
    if not success:
        return set()

    # Parse PAF and compute coverage + identity per query
    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    query_matches: Dict[str, int] = defaultdict(int)
    query_alnlen: Dict[str, int] = defaultdict(int)

    with paf_out.open("r") as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            qname = fields[0]
            try:
                qs = int(fields[2])
                qe = int(fields[3])
                matches = int(fields[9])
                aln_len = int(fields[10])
            except ValueError:
                continue

            if qs > qe:
                qs, qe = qe, qs

            query_intervals[qname].append((qs, qe))
            query_matches[qname] += matches
            query_alnlen[qname] += aln_len

    # Identify debris: high coverage and identity
    debris_contigs: Set[str] = set()
    for qname, intervals in query_intervals.items():
        qlen = query_lengths.get(qname, 0)
        if qlen == 0:
            continue

        _, total_bp = merge_intervals(intervals)
        coverage = total_bp / qlen

        total_matches = query_matches[qname]
        total_alnlen = query_alnlen[qname]
        identity = (total_matches / total_alnlen) if total_alnlen > 0 else 0.0

        if coverage >= min_coverage and identity >= min_identity:
            debris_contigs.add(qname)
            print(
                f"[info] Chromosome debris: {qname} (cov={coverage:.2f}, ident={identity:.2f})",
                file=sys.stderr,
            )

    print(f"[info] Chromosome debris contigs: {len(debris_contigs)}", file=sys.stderr)
    return debris_contigs


# ----------------------------
# Contaminant detection (centrifuger)
# ----------------------------
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
) -> Dict[str, Tuple[int, str]]:
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
        Dict mapping contig_name -> (taxid, scientific_name)
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
    contaminants: Dict[str, Tuple[int, str]] = {}

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

            contaminants[contig_name] = (taxid, sci_name)
            print(
                f"[info] Contaminant: {contig_name} (taxid={taxid}, {sci_name}, score={score}, cov={coverage:.2f})",
                file=sys.stderr,
            )

    print(f"[info] Contaminant contigs: {len(contaminants)}", file=sys.stderr)
    return contaminants


# ----------------------------
# Orientation determination
# ----------------------------
def compute_orientation_votes(
    macro_block_rows: List[Tuple],
    contig: str,
    assigned_ref_id: str,
) -> Tuple[int, int]:
    """Count forward/reverse orientation votes from macro blocks.

    Each macro block provides one vote based on its strand.

    Returns (fwd_count, rev_count)
    """
    fwd_count = 0
    rev_count = 0

    for row in macro_block_rows:
        # Row format: (contig, contig_len, ref_id, chrom_id, subgenome, strand, ...)
        if len(row) < 6:
            continue
        row_contig = row[0]
        row_ref_id = row[2]
        row_strand = row[5]

        if row_contig != contig:
            continue
        if row_ref_id != assigned_ref_id:
            continue

        if row_strand == "+":
            fwd_count += 1
        elif row_strand == "-":
            rev_count += 1

    return fwd_count, rev_count


def determine_contig_orientations(
    macro_block_rows: List[Tuple],
    best_ref: Dict[str, str],
    chromosome_contigs: Set[str],
) -> Dict[str, bool]:
    """Determine which chromosome contigs need to be reverse-complemented.

    Only chromosome-assigned contigs are subject to reorientation based on synteny
    block strand votes. Non-chromosome contigs (debris, contaminants, unclassified)
    are left in their original orientation since they lack reliable synteny evidence.

    Returns dict: contig_name -> should_reverse_complement (only for chromosome contigs)
    """
    orientations: Dict[str, bool] = {}

    for contig in chromosome_contigs:
        assigned_ref = best_ref.get(contig, "")
        if not assigned_ref:
            orientations[contig] = False
            continue

        fwd, rev = compute_orientation_votes(macro_block_rows, contig, assigned_ref)

        # Reverse if more reverse votes than forward
        should_reverse = rev > fwd
        orientations[contig] = should_reverse

        if should_reverse:
            print(f"[info] Contig {contig} will be reverse-complemented (fwd={fwd}, rev={rev})", file=sys.stderr)

    return orientations


# ----------------------------
# Gene count statistics
# ----------------------------
def count_genes_per_ref_chrom(
    gff3_path: Path,
    ref_id_map: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """Count genes per reference chromosome from GFF3.

    Counts unique gene IDs per chromosome.

    Returns dict: ref_id -> gene_count
    """
    chrom_genes: Dict[str, Set[str]] = defaultdict(set)

    with open_maybe_gzip(gff3_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _src, ftype, _start, _end, _score, _strand, _phase, attrs = parts

            if ftype.lower() != "gene":
                continue

            # Parse gene ID from attributes
            gene_id = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("ID="):
                    gene_id = kv[3:]
                    break

            if gene_id:
                ref_id = ref_id_map.get(seqid, normalize_ref_id(seqid)) if ref_id_map else seqid
                chrom_genes[ref_id].add(gene_id)

    return {chrom: len(genes) for chrom, genes in chrom_genes.items()}


def compute_mean_genes_per_Mbp(
    qr_gene_count: Dict[Tuple[str, str], int],
    query_lengths: Dict[str, int],
    chromosome_contigs: Set[str],
    best_ref: Dict[str, str],
) -> float:
    """Compute mean genes per Mbp across chromosome-assigned contigs."""
    total_genes = 0
    total_bp = 0

    for contig in chromosome_contigs:
        ref_id = best_ref.get(contig, "")
        if not ref_id:
            continue
        gene_count = qr_gene_count.get((contig, ref_id), 0)
        contig_len = query_lengths.get(contig, 0)

        total_genes += gene_count
        total_bp += contig_len

    if total_bp == 0:
        return 0.0

    mean_gpmbp = total_genes / (total_bp / 1_000_000.0)
    print(f"[info] Mean genes per Mbp (query): {mean_gpmbp:.2f} ({total_genes} genes / {total_bp} bp)", file=sys.stderr)
    return mean_gpmbp


# ----------------------------
# Reference-based debris classification
# ----------------------------
def classify_debris_and_unclassified(
    remaining_contigs: Set[str],
    query_fasta: Path,
    query_lengths: Dict[str, int],
    ref_fasta: Path,
    chrs_fasta: Optional[Path],
    qr_gene_count: Dict[Tuple[str, str], int],
    work_dir: Path,
    threads: int,
    min_coverage: float,
    min_protein_hits: int,
) -> Tuple[Set[str], Set[str]]:
    """Classify remaining contigs as debris or unclassified based on reference homology.

    MOTIVATION: After chromosome debris detection (which catches assembly artifacts that
    are near-identical to assembled chromosomes), some contigs may still have homology to
    the reference genome without being classified. These include:
    - Divergent haplotigs that don't align well at nucleotide level to assembled chromosomes
      but retain protein-coding content (detected via miniprot hits)
    - Sequences from genomic regions not represented in the primary chromosome contigs
      but present in the reference (detected via reference nucleotide alignment)

    This reference-based approach complements chromosome debris detection by catching
    sequences with ancestral/evolutionary homology rather than assembly-specific duplicates.

    Classification criteria:
    - Debris: contigs with 50%+ nucleotide alignment coverage to reference OR 2+ miniprot hits
    - Unclassified: everything else (potential novel sequences, contaminants missed by
      earlier screening, or highly divergent sequences)

    Args:
        remaining_contigs: Set of contigs not yet classified
        query_fasta: Path to query FASTA
        query_lengths: Dict of contig name -> length
        ref_fasta: Path to reference FASTA
        chrs_fasta: Optional path to assembled chromosome FASTA (currently unused)
        qr_gene_count: Dict of (contig, ref_id) -> gene count from miniprot
        work_dir: Working directory
        threads: Number of threads
        min_coverage: Min alignment coverage fraction for debris classification
        min_protein_hits: Min miniprot protein hits for debris classification

    Returns:
        Tuple of (debris_contigs, unclassified_contigs)
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    debris_contigs: Set[str] = set()
    unclassified_contigs: Set[str] = set()

    # Check protein hits for each remaining contig
    contigs_with_protein_hits: Set[str] = set()
    for (contig, _), gene_count in qr_gene_count.items():
        if contig in remaining_contigs and gene_count >= min_protein_hits:
            contigs_with_protein_hits.add(contig)

    # For contigs without enough protein hits, check nucleotide alignment
    contigs_needing_alignment = remaining_contigs - contigs_with_protein_hits

    if contigs_needing_alignment:
        # Write subset of contigs to align
        subset_fasta = work_dir / "debris_candidates.fa"
        query_seqs = read_fasta_sequences(query_fasta)
        subset_seqs = {k: v for k, v in query_seqs.items() if k in contigs_needing_alignment}
        if subset_seqs:
            write_fasta(subset_seqs, subset_fasta)

            # Align against reference using minimap2/mm2plus if available
            if get_minimap2_exe() and subset_fasta.exists():
                paf_out = work_dir / "debris_vs_ref.paf"

                run_minimap2(
                    ref=ref_fasta,
                    qry=subset_fasta,
                    paf_out=paf_out,
                    threads=threads,
                    preset="asm20",
                )

                # Parse PAF to get coverage
                if paf_out.exists():
                    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
                    with paf_out.open("r") as fh:
                        for line in fh:
                            fields = line.rstrip("\n").split("\t")
                            if len(fields) < 12:
                                continue
                            qname = fields[0]
                            try:
                                qs = int(fields[2])
                                qe = int(fields[3])
                            except ValueError:
                                continue
                            if qs > qe:
                                qs, qe = qe, qs
                            query_intervals[qname].append((qs, qe))

                    for qname, intervals in query_intervals.items():
                        _, total_bp = merge_intervals(intervals)
                        qlen = query_lengths.get(qname, 0)
                        coverage = (total_bp / qlen) if qlen > 0 else 0.0
                        if coverage >= min_coverage:
                            debris_contigs.add(qname)

    # Add contigs with protein hits to debris
    debris_contigs.update(contigs_with_protein_hits)

    # Everything else is unclassified
    unclassified_contigs = remaining_contigs - debris_contigs

    print(f"[info] Debris contigs: {len(debris_contigs)}", file=sys.stderr)
    print(f"[info] Unclassified contigs: {len(unclassified_contigs)}", file=sys.stderr)

    return debris_contigs, unclassified_contigs


# ----------------------------
# Contig naming
# ----------------------------
def generate_contig_names(
    classifications: List[ContigClassification],
    query_lengths: Dict[str, int],
    add_subgenome_suffix: Optional[str],
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """Generate new contig names based on classification.

    - Chromosome contigs: inherit reference name (e.g., chr5T)
    - Multiple contigs to same ref: add _1, _2 suffix by descending length
    - Non-chromosome: contig_N (N increases with decreasing length)
    """
    name_mapping: Dict[str, str] = {}

    # Group chromosome contigs by assigned reference
    ref_to_contigs: Dict[str, List[str]] = defaultdict(list)
    non_chrom_contigs: List[str] = []

    for clf in classifications:
        if clf.classification == "chrom" and clf.assigned_ref_id:
            ref_to_contigs[clf.assigned_ref_id].append(clf.original_name)
        else:
            non_chrom_contigs.append(clf.original_name)

    # Name chromosome contigs
    for ref_id, contigs in ref_to_contigs.items():
        # Sort by length descending
        contigs.sort(key=lambda c: query_lengths.get(c, 0), reverse=True)

        base_name = ref_norm_to_orig.get(ref_id, ref_id) if ref_norm_to_orig else ref_id
        # Add subgenome suffix if reference doesn't have one and user requested it
        if add_subgenome_suffix:
            chrom_id, sub = split_chrom_subgenome(base_name)
            if sub == "NA":
                base_name = f"{chrom_id}{add_subgenome_suffix}"

        if len(contigs) == 1:
            name_mapping[contigs[0]] = base_name
        else:
            for i, contig in enumerate(contigs, start=1):
                name_mapping[contig] = f"{base_name}_{i}"

    # Name non-chromosome contigs
    # Sort by length descending
    non_chrom_contigs.sort(key=lambda c: query_lengths.get(c, 0), reverse=True)
    for i, contig in enumerate(non_chrom_contigs, start=1):
        name_mapping[contig] = f"contig_{i}"

    return name_mapping


# ----------------------------
# FASTA output
# ----------------------------
def write_classified_fastas(
    query_fasta: Path,
    classifications: List[ContigClassification],
    contig_orientations: Dict[str, bool],
    output_prefix: Path,
) -> Dict[str, Path]:
    """Write 6 classified FASTA files.

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
        if clf.classification == "chrom":
            category = "chrs"
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

        # Sort by chromosome number or contig number
        def sort_key(name: str):
            # Extract numeric parts for sorting
            m = re.match(r"[cC]hr(\d+)", name)
            if m:
                return (0, int(m.group(1)), name)
            m = re.match(r"contig_(\d+)", name)
            if m:
                return (1, int(m.group(1)), name)
            return (2, 0, name)

        sorted_seqs = dict(sorted(out_seqs.items(), key=lambda x: sort_key(x[0])))
        write_fasta(sorted_seqs, output_path)

        print(f"[done] {category}: {output_path} ({len(sorted_seqs)} contigs)", file=sys.stderr)

    return output_paths


# ----------------------------
# Classification pipeline
# ----------------------------
def classify_all_contigs(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    best_ref: Dict[str, str],
    chr_like_minlen: int,
    ev: ChainEvidenceResult,
    ref_gene_counts: Dict[str, int],
    chrC_contig: Optional[str],
    chrM_contig: Optional[str],
    organelle_debris: Set[str],
    rdna_contigs: Set[str],
    contaminants: Dict[str, Tuple[int, str]],
    chromosome_debris: Set[str],
    other_debris: Set[str],
    add_subgenome_suffix: Optional[str],
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> List[ContigClassification]:
    """Classify all contigs and generate classifications list.

    Args:
        query_fasta: Path to query FASTA
        query_lengths: Dict of contig name -> length
        best_ref: Dict of contig -> best reference assignment
        chr_like_minlen: Minimum length for chromosome classification
        ev: Chain evidence results
        ref_gene_counts: Dict of ref_id -> gene count
        chrC_contig: Chloroplast contig name
        chrM_contig: Mitochondrial contig name
        organelle_debris: Set of organelle debris contigs
        rdna_contigs: Set of rDNA contigs
        contaminants: Dict of contig -> (taxid, scientific_name)
        chromosome_debris: Set of chromosome debris contigs (from chr-vs-chr alignment)
        other_debris: Set of other debris contigs (from ref/protein alignment)
        add_subgenome_suffix: Optional subgenome suffix to add
    """

    classifications: List[ContigClassification] = []
    classified_contigs: Set[str] = set()

    # 1. Chromosome contigs (status OK and length >= chr_like_minlen)
    for contig, ref_id in best_ref.items():
        if not ref_id:
            continue
        contig_len = query_lengths.get(contig, 0)
        if contig_len < chr_like_minlen:
            continue

        # Check if it passed synteny gates
        gene_count = ev.qr_gene_count.get((contig, ref_id), 0)
        ref_total_genes = ref_gene_counts.get(ref_id, 0)
        gene_proportion = (gene_count / ref_total_genes) if ref_total_genes > 0 else None

        classifications.append(ContigClassification(
            original_name=contig,
            new_name="",  # Will be filled later
            classification="chrom",
            reversed=False,  # Will be filled later
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id=ref_id,
            ref_gene_proportion=gene_proportion,
            contig_len=contig_len,
        ))
        classified_contigs.add(contig)

    # 2. Organelle contigs
    if chrC_contig and chrC_contig not in classified_contigs:
        classifications.append(ContigClassification(
            original_name=chrC_contig,
            new_name="chrC",
            classification="organelle_complete",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chrC",
            ref_gene_proportion=None,
            contig_len=query_lengths.get(chrC_contig, 0),
        ))
        classified_contigs.add(chrC_contig)

    if chrM_contig and chrM_contig not in classified_contigs:
        classifications.append(ContigClassification(
            original_name=chrM_contig,
            new_name="chrM",
            classification="organelle_complete",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chrM",
            ref_gene_proportion=None,
            contig_len=query_lengths.get(chrM_contig, 0),
        ))
        classified_contigs.add(chrM_contig)

    # Organelle debris
    for contig in organelle_debris:
        if contig not in classified_contigs:
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="organelle_debris",
                reversed=False,
                contaminant_taxid=None,
            contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
            ))
            classified_contigs.add(contig)

    # 3. rDNA contigs
    for contig in rdna_contigs:
        if contig not in classified_contigs:
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="rDNA",
                reversed=False,
                contaminant_taxid=None,
            contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
            ))
            classified_contigs.add(contig)

    # 4. Contaminants
    for contig, (taxid, sci_name) in contaminants.items():
        if contig not in classified_contigs:
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="contaminant",
                reversed=False,
                contaminant_taxid=taxid,
                contaminant_sci=sci_name,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
            ))
            classified_contigs.add(contig)

    # 5. Chromosome debris (from chr-vs-chr alignment detection)
    for contig in chromosome_debris:
        if contig not in classified_contigs:
            ref_id = best_ref.get(contig, "")
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="chrom_debris",
                reversed=False,
                contaminant_taxid=None,
            contaminant_sci=None,
                assigned_ref_id=ref_id if ref_id else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
            ))
            classified_contigs.add(contig)

    # 6. Other debris (from reference/protein alignment detection)
    for contig in other_debris:
        if contig not in classified_contigs:
            # Check if this contig has synteny support
            ref_id = best_ref.get(contig, "")
            classification = "chrom_debris" if ref_id else "debris"

            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification=classification,
                reversed=False,
                contaminant_taxid=None,
            contaminant_sci=None,
                assigned_ref_id=ref_id if ref_id else None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
            ))
            classified_contigs.add(contig)

    # 7. Unclassified (everything else)
    for contig in query_lengths.keys():
        if contig not in classified_contigs:
            classifications.append(ContigClassification(
                original_name=contig,
                new_name="",
                classification="unclassified",
                reversed=False,
                contaminant_taxid=None,
            contaminant_sci=None,
                assigned_ref_id=None,
                ref_gene_proportion=None,
                contig_len=query_lengths.get(contig, 0),
            ))

    # Generate names
    name_mapping = generate_contig_names(
        classifications,
        query_lengths,
        add_subgenome_suffix,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    # Update classifications with names
    for clf in classifications:
        if not clf.new_name:
            clf.new_name = name_mapping.get(clf.original_name, clf.original_name)

    return classifications


# ----------------------------
# Enhanced summary TSV
# ----------------------------
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

    New columns added:
    - original_name
    - classification
    - reversed
    - contaminant_taxid (NCBI taxonomy ID for contaminants)
    - contaminant_sci (scientific name for contaminants)
    - length (after status)
    - ref_gene_proportion
    - genes_per_Mbp (for chromosome-assigned contigs)
    """
    # Build lookup by original name
    clf_lookup = {clf.original_name: clf for clf in classifications}

    header = [
        "contig",
        "original_name",
        "classification",
        "reversed",
        "contaminant_taxid",
        "contaminant_sci",
        "assigned_subgenome",
        "assigned_ref_id",
        "assigned_chrom_id",
        "status",
        "length",
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

            assigned_ref_id = str(best_ref.get(q, "") or "")
            bs = float(best_score.get(q, 0.0) or 0.0)
            sr = float(second_score.get(q, 0.0) or 0.0)

            best_union_bp = int(best_bp.get(q, 0) or 0)
            second_ref_id = str(second_ref.get(q, "") or "")
            second_union_bp = int(second_bp.get(q, 0) or 0)

            _mtmp, m_b, al_b, ident_b, dist_b = fetch_metrics(q, assigned_ref_id)
            _stmp, m_s, al_s, ident_s, dist_s = fetch_metrics(q, second_ref_id)

            if not assigned_ref_id or bs <= 0.0:
                assigned_chrom_id, assigned_subgenome = ("", "NA")
                status = "NO_HITS"
                best_ref_union_frac = 0.0
                score_ratio = 0.0
                assigned_ref_id_out = "NA"
            else:
                assigned_ref_id_out = _denormalize_ref_id(assigned_ref_id, ref_norm_to_orig)
                assigned_chrom_id, assigned_subgenome = split_chrom_subgenome(assigned_ref_id_out)
                best_ref_union_frac = (best_union_bp / contig_len) if contig_len > 0 else 0.0
                score_ratio = (bs / sr) if sr > 0 else (float("inf") if bs > 0 else 0.0)

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
            reversed_val = "yes" if contig_orientations.get(q, False) else "no"
            contaminant_taxid = str(clf.contaminant_taxid) if clf and clf.contaminant_taxid is not None else ""
            contaminant_sci = clf.contaminant_sci if clf and clf.contaminant_sci else ""
            ref_gene_proportion = clf.ref_gene_proportion if clf and clf.ref_gene_proportion is not None else ""

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
                        reversed_val,
                        contaminant_taxid,
                        contaminant_sci,
                        str(assigned_subgenome),
                        str(assigned_ref_id_out),
                        str(assigned_chrom_id),
                        str(status),
                        str(int(contig_len)),  # length
                        f"{bs:.3f}",
                        f"{sr:.3f}",
                        (f"{score_ratio:.3f}" if not math.isinf(score_ratio) else "Inf"),
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


def _diag_for_block(b: Block) -> int:
    return (b.rs - b.qs) if b.strand == "+" else (b.rs + b.qs)


def _rpos_for_chain_compat(b: Block) -> int:
    return b.rs if b.strand == "+" else -b.re


def _finalize_chain(chain: Chain):
    merged, qbp = merge_intervals(chain.q_intervals)
    return merged, qbp, chain.matches_sum, chain.alnlen_sum


def _chain_weight(qbp: int, matches: int, alnlen: int, mode: str) -> float:
    if alnlen <= 0:
        return 0.0
    ident = matches / alnlen
    if mode == "matches":
        return float(matches)
    if mode == "qbp_ident":
        return float(qbp) * ident
    if mode == "matches_ident":
        return float(matches) * ident
    raise ValueError(f"Unknown chain score mode: {mode}")


def _tp_keep(fields: list[str], assign_tp: str) -> bool:
    if assign_tp == "ALL":
        return True
    tp = paf_tag_value(fields, "tp:A:")
    if tp is None:
        return True
    if assign_tp == "P":
        return tp == "P"
    return tp in ("P", "I")


def parse_paf_chain_evidence_and_segments(
    paf_gz_path: Path,
    contig_lengths: dict,
    assign_minlen: int,
    assign_minmapq: int,
    assign_tp: str,
    chain_q_gap: int,
    chain_r_gap: int,
    chain_diag_slop: int,
    assign_min_ident: float,
    assign_chain_topk: int,
    assign_chain_score: str,
    assign_chain_min_bp: int,
    assign_ref_score: str,
    ref_id_map: Optional[Dict[str, str]] = None,
) -> ChainEvidenceResult:
    """
    Original nucleotide-chain mode (mm2plus PAF):
      contig=query, ref_id=reference sequence (typically chrNA/chrNP/etc.).
    """
    if assign_chain_topk < 1:
        raise ValueError("--assign-chain-topk must be >= 1")

    blocks = defaultdict(list)  # (q, ref_id, strand) -> [Block,...]
    qlens_from_paf: dict[str, int] = {}

    with open_maybe_gzip(paf_gz_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            if not _tp_keep(fields, assign_tp):
                continue

            qname = fields[0]
            try:
                qlen = int(fields[1])
                qs = int(fields[2])
                qe = int(fields[3])
                strand = fields[4]
                ref_id_raw = fields[5]
                rs = int(fields[7])
                re_ = int(fields[8])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            ref_id = ref_id_map.get(ref_id_raw, normalize_ref_id(ref_id_raw)) if ref_id_map else normalize_ref_id(ref_id_raw)

            if aln_len < assign_minlen:
                continue
            if mapq < assign_minmapq:
                continue

            if qe < qs:
                qs, qe = qe, qs
            if qe <= qs:
                continue

            if re_ < rs:
                rs, re_ = re_, rs

            qlens_from_paf.setdefault(qname, qlen)
            blocks[(qname, ref_id, strand)].append(
                Block(
                    qs=qs,
                    qe=qe,
                    rs=rs,
                    re=re_,
                    matches=matches,
                    aln_len=aln_len,
                    mapq=mapq,
                    strand=strand,
                    gene_id=None,
                )
            )

    return _chains_to_evidence_and_segments(
        blocks=blocks,
        contig_lengths=contig_lengths,
        qlens_from_paf=qlens_from_paf,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
        segments_strand_from_blocks=True,
    )


def parse_miniprot_synteny_evidence_and_segments(
    miniprot_paf_gz: Path,
    contig_lengths: dict,
    tx2loc: dict,
    tx2gene: dict,
    assign_minlen: int,
    assign_minmapq: int,
    chain_q_gap: int,
    chain_r_gap: int,
    chain_diag_slop: int,
    assign_min_ident: float,
    assign_chain_topk: int,
    assign_chain_score: str,
    assign_chain_min_bp: int,
    assign_ref_score: str,
) -> ChainEvidenceResult:
    """
    Protein-anchor mode:
      - miniprot PAF has: qname=protein/transcript, tname=query contig
      - Map qname -> ref_id(ref chrom) + (start,end,strand) from GFF3 transcript features
      - Build blocks in (query-contig coordinate) vs (reference transcript genomic coordinate) space
    """
    if assign_chain_topk < 1:
        raise ValueError("--assign-chain-topk must be >= 1")

    blocks = defaultdict(list)  # (contig, ref_id, strand) -> [Block,...]
    qlens_from_paf: dict[str, int] = {}  # contig -> tlen from PAF

    with open_maybe_gzip(miniprot_paf_gz, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            prot = fields[0]
            prot = prot.replace("mRNA:", "").replace("transcript:", "")
            tname = fields[5]  # target sequence = query contig

            try:
                tlen = int(fields[6])
                ts = int(fields[7])
                te = int(fields[8])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue

            # Map protein/transcript -> reference location
            prot_keys = [prot, prot.split()[0], prot.split("|")[0]]
            loc = None
            gene_id = None
            for pk in prot_keys:
                if pk in tx2loc:
                    loc = tx2loc.get(pk)
                    gene_id = tx2gene.get(pk) or pk
                    break
            if loc is None:
                continue

            ref_id, rs0, re0, _rstrand = loc

            if te < ts:
                ts, te = te, ts
            if te <= ts:
                continue

            # In protein mode, interpret assign_minlen as minimum TARGET SPAN on the contig.
            tspan = te - ts
            if tspan < assign_minlen:
                continue
            if mapq < assign_minmapq:
                continue

            ident = (matches / aln_len) if aln_len > 0 else 0.0
            if ident < assign_min_ident:
                continue

            # Chaining in increasing "ref coordinate" space; force '+'.
            chain_strand = "+"

            qlens_from_paf.setdefault(tname, tlen)
            blocks[(tname, ref_id, chain_strand)].append(
                Block(
                    qs=ts,
                    qe=te,
                    rs=rs0,
                    re=re0,
                    matches=matches,
                    aln_len=aln_len,
                    mapq=mapq,
                    strand=chain_strand,
                    gene_id=gene_id,
                )
            )

    return _chains_to_evidence_and_segments(
        blocks=blocks,
        contig_lengths=contig_lengths,
        qlens_from_paf=qlens_from_paf,
        chain_q_gap=chain_q_gap,
        chain_r_gap=chain_r_gap,
        chain_diag_slop=chain_diag_slop,
        assign_min_ident=assign_min_ident,
        assign_chain_topk=assign_chain_topk,
        assign_chain_score=assign_chain_score,
        assign_chain_min_bp=assign_chain_min_bp,
        assign_ref_score=assign_ref_score,
        segments_strand_from_blocks=False,  # forced '+'
    )


def _chains_to_evidence_and_segments(
    blocks,
    contig_lengths,
    qlens_from_paf,
    chain_q_gap,
    chain_r_gap,
    chain_diag_slop,
    assign_min_ident,
    assign_chain_topk,
    assign_chain_score,
    assign_chain_min_bp,
    assign_ref_score: str,
    segments_strand_from_blocks: bool,
) -> ChainEvidenceResult:
    # Evidence across kept chains
    qr_intervals = defaultdict(list)  # (q, ref_id) -> list[(qs,qe)]
    qr_matches = defaultdict(int)
    qr_alnlen = defaultdict(int)
    qr_gene_ids = defaultdict(set)  # (q, ref_id) -> set(gene_id)

    # Assignment weights
    qr_chain_weights = defaultdict(list)
    qr_best_chain_qbp = defaultdict(int)
    qr_best_chain_ident = defaultdict(float)
    qr_best_chain_weight = defaultdict(float)
    qr_nchains_kept = defaultdict(int)
    qr_weight_all = defaultdict(float)

    chain_segments_rows = []
    macro_block_rows = []

    for (q, ref_id, strand), blks in blocks.items():
        blks.sort(key=lambda b: (b.qs, b.qe))
        chains: list[Chain] = []

        for b in blks:
            diag = _diag_for_block(b)
            rpos = _rpos_for_chain_compat(b)

            best_i = None
            best_penalty = None

            for i, ch in enumerate(chains):
                qgap = b.qs - ch.last_qe
                if qgap > chain_q_gap:
                    continue

                rgap = rpos - ch.last_r
                if rgap < -chain_r_gap or rgap > chain_r_gap:
                    continue

                if abs(diag - ch.diag0) > chain_diag_slop:
                    continue

                penalty = abs(qgap) + abs(rgap) + (abs(diag - ch.diag0) // 10)
                if best_penalty is None or penalty < best_penalty:
                    best_penalty = penalty
                    best_i = i

            if best_i is None:
                gset: set[str] = set()
                if b.gene_id:
                    gset.add(str(b.gene_id))
                chains.append(
                    Chain(
                        last_qe=b.qe,
                        last_r=rpos,
                        diag0=diag,
                        q_intervals=[(b.qs, b.qe)],
                        matches_sum=b.matches,
                        alnlen_sum=b.aln_len,
                        gene_ids=gset,
                    )
                )
            else:
                ch = chains[best_i]
                if b.gene_id:
                    ch.gene_ids.add(str(b.gene_id))
                ch.last_qe = max(ch.last_qe, b.qe)
                ch.last_r = rpos
                ch.q_intervals.append((b.qs, b.qe))
                ch.matches_sum += b.matches
                ch.alnlen_sum += b.aln_len

        qlen = int(contig_lengths.get(q, 0) or qlens_from_paf.get(q, 0) or 0)
        chrom_id, sub = split_chrom_subgenome(ref_id)

        for chain_id, ch in enumerate(chains, start=1):
            merged, qbp, msum, alnsum = _finalize_chain(ch)

            if qbp <= 0 or alnsum <= 0 or not merged:
                continue
            if qbp < assign_chain_min_bp:
                continue

            ident = msum / alnsum
            if ident < assign_min_ident:
                continue

            w = _chain_weight(qbp=qbp, matches=msum, alnlen=alnsum, mode=assign_chain_score)
            if w <= 0:
                continue

            seg_strand = strand if segments_strand_from_blocks else "+"

            # Macro block spanning the full merged chain extent on the contig
            qstart = min(s for s, _e in merged)
            qend = max(e for _s, e in merged)
            qspan = max(0, qend - qstart)
            nseg = len(merged)
            gene_count_chain = len(ch.gene_ids) if ch.gene_ids else 0

            macro_block_rows.append(
                (
                    q,
                    qlen,
                    ref_id,
                    chrom_id,
                    sub,
                    seg_strand,
                    chain_id,
                    qstart,
                    qend,
                    qspan,
                    qbp,
                    msum,
                    alnsum,
                    f"{ident:.6f}",
                    f"{w:.3f}",
                    nseg,
                    gene_count_chain,
                )
            )

            key_qr = (q, ref_id)
            if ch.gene_ids:
                qr_gene_ids[key_qr].update(ch.gene_ids)
            qr_intervals[key_qr].extend(merged)
            qr_matches[key_qr] += msum
            qr_alnlen[key_qr] += alnsum

            qr_chain_weights[key_qr].append(w)
            qr_weight_all[key_qr] += w
            qr_nchains_kept[key_qr] += 1

            if w > qr_best_chain_weight[key_qr]:
                qr_best_chain_weight[key_qr] = w
                qr_best_chain_qbp[key_qr] = qbp
                qr_best_chain_ident[key_qr] = ident

            for (s, e) in merged:
                chain_segments_rows.append((q, qlen, chrom_id, sub, seg_strand, chain_id, s, e))

    # finalize per-(q,ref_id) union bp
    qr_union_bp = defaultdict(int)
    for key, ivs in qr_intervals.items():
        _m, total = merge_intervals(ivs)
        qr_union_bp[key] = total

    contig_total = defaultdict(int)
    contig_refs = defaultdict(set)
    for (q, ref_id), ubp in qr_union_bp.items():
        contig_total[q] += ubp
        contig_refs[q].add(ref_id)

    qr_score_topk = defaultdict(float)
    for key, ws in qr_chain_weights.items():
        ws_sorted = sorted(ws, reverse=True)
        qr_score_topk[key] = float(sum(ws_sorted[:assign_chain_topk]))

    if assign_ref_score == "all":
        qr_ref_score = qr_weight_all
    elif assign_ref_score == "topk":
        qr_ref_score = qr_score_topk
    else:
        raise ValueError(f"Unknown assign_ref_score: {assign_ref_score}")

    best_ref = defaultdict(str)
    best_score = defaultdict(float)
    best_bp = defaultdict(int)

    second_ref = defaultdict(str)
    second_score = defaultdict(float)
    second_bp = defaultdict(int)

    for (q, ref_id), score in qr_ref_score.items():
        ubp = int(qr_union_bp.get((q, ref_id), 0) or 0)
        if score > best_score[q]:
            second_score[q] = best_score[q]
            second_ref[q] = best_ref[q]
            second_bp[q] = best_bp[q]

            best_score[q] = score
            best_ref[q] = ref_id
            best_bp[q] = ubp
        elif score > second_score[q]:
            second_score[q] = score
            second_ref[q] = ref_id
            second_bp[q] = ubp

    chain_summary_rows = []
    for (q, ref_id), _ubp in qr_union_bp.items():
        ubp = int(qr_union_bp.get((q, ref_id), 0) or 0)
        msum = int(qr_matches.get((q, ref_id), 0) or 0)
        alnsum = int(qr_alnlen.get((q, ref_id), 0) or 0)
        ident = (msum / alnsum) if alnsum > 0 else 0.0
        gene_count_ref = len(qr_gene_ids.get((q, ref_id), set()))

        w_all = float(qr_weight_all.get((q, ref_id), 0.0) or 0.0)
        w_topk = float(qr_score_topk.get((q, ref_id), 0.0) or 0.0)
        nch = int(qr_nchains_kept.get((q, ref_id), 0) or 0)
        bqbp = int(qr_best_chain_qbp.get((q, ref_id), 0) or 0)
        bident = float(qr_best_chain_ident.get((q, ref_id), 0.0) or 0.0)
        bw = float(qr_best_chain_weight.get((q, ref_id), 0.0) or 0.0)

        chrom_id, sub = split_chrom_subgenome(ref_id)
        chain_summary_rows.append(
            (
                q,
                ref_id,
                chrom_id,
                sub,
                ubp,
                msum,
                alnsum,
                ident,
                gene_count_ref,
                nch,
                w_topk,
                w_all,
                bw,
                bqbp,
                bident,
            )
        )

    qr_gene_count = {k: len(v) for k, v in qr_gene_ids.items()}

    return ChainEvidenceResult(
        qlens_from_paf=dict(qlens_from_paf),
        qr_union_bp=dict(qr_union_bp),
        qr_matches=dict(qr_matches),
        qr_alnlen=dict(qr_alnlen),
        qr_gene_count=qr_gene_count,
        contig_total=dict(contig_total),
        contig_refs=dict(contig_refs),
        qr_score_topk=dict(qr_score_topk),
        qr_weight_all=dict(qr_weight_all),
        qr_nchains_kept=dict(qr_nchains_kept),
        best_ref=dict(best_ref),
        best_score=dict(best_score),
        best_bp=dict(best_bp),
        second_ref=dict(second_ref),
        second_score=dict(second_score),
        second_bp=dict(second_bp),
        chain_segments_rows=chain_segments_rows,
        macro_block_rows=macro_block_rows,
        chain_summary_rows=chain_summary_rows,
    )


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
                score_ratio = (bs / sr) if sr > 0 else (float("inf") if bs > 0 else 0.0)

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
                        (f"{score_ratio:.3f}" if not math.isinf(score_ratio) else "Inf"),
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


def have_rscript() -> bool:
    try:
        subprocess.run(["Rscript", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception:
        return False


def _read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8") as fh:
        return fh.read()


def run_plot(
    summary_tsv: Path,
    ref_lengths_tsv: Path,
    segments_tsv: Path,
    chain_summary_tsv: Path,
    macro_blocks_tsv: Path,
    outprefix: Path,
    chr_like_minlen: int,
    plot_title_suffix: str,
):
    if not have_rscript():
        print("[warn] --plot specified, but Rscript not found in PATH; skipping plot.", file=sys.stderr)
        return

    plot_pdf = Path(str(outprefix) + ".subgenome_assignment_overview.pdf")
    r_script_path = Path(str(outprefix) + ".subgenome_assignment_overview.R")

    # Template expected alongside this Python file
    tmpl_path = Path(__file__).resolve().with_name("subgenome_assignment_overview.tmpl.R")
    if not tmpl_path.exists():
        raise FileNotFoundError(
            f"Missing R template: {tmpl_path}\n"
            "Create it next to the Python script (see instructions)."
        )

    tmpl = _read_text(tmpl_path)

    # Very simple placeholder replacement (safe + robust)
    def esc(s: Path | str) -> str:
        # ensure backslashes don't accidentally escape quotes on Windows
        return str(s).replace("\\", "/")

    filled = (
        tmpl.replace("__SUMMARY__", esc(summary_tsv))
        .replace("__REF__", esc(ref_lengths_tsv))
        .replace("__SEGMENTS__", esc(segments_tsv))
        .replace("__EVIDENCE__", esc(chain_summary_tsv))
        .replace("__MACRO__", esc(macro_blocks_tsv))
        .replace("__OUTPDF__", esc(plot_pdf))
        .replace("__CHRLIKE__", str(int(chr_like_minlen)))
        .replace("__SUFFIX__", str(plot_title_suffix).replace('"', '\\"'))
    )

    with r_script_path.open("w", encoding="utf-8") as rf:
        rf.write(filled)

    print(f"[info] Running Rscript to generate plot: {plot_pdf}", file=sys.stderr)
    try:
        subprocess.run(["Rscript", str(r_script_path)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[warn] Rscript failed with code {e.returncode}; plot not generated.", file=sys.stderr)
    else:
        print(f"[done] Plot written to: {plot_pdf}", file=sys.stderr)


# ----------------------------
# Main
# ----------------------------
def main():
    p = argparse.ArgumentParser(
        description="Genome-wide subgenome distance analysis using minimap2 chains or miniprot protein-anchored synteny blocks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # =========================================================================
    # Required arguments
    # =========================================================================
    req = p.add_argument_group("Required arguments")
    req.add_argument("-r", "--ref", required=True, help="Reference assembly FASTA (plain or gzipped; chrNA/chrNP/etc.)")
    req.add_argument("-q", "--query", required=True, help="Query assembly FASTA (plain or gzipped; contigs/scaffolds)")
    req.add_argument("-o", "--outprefix", required=True, help="Output prefix (no extension)")
    req.add_argument(
        "--ref-gff3",
        required=True,
        help="Reference GFF3 with protein-coding gene annotations. "
        "Used to extract proteins (gffread) and run protein-anchor synteny blocks (miniprot).",
    )

    # =========================================================================
    # Common options
    # =========================================================================
    common = p.add_argument_group("Common options")
    common.add_argument("-t", "--threads", type=int, default=8, help="Threads for minimap2/miniprot [8]")
    common.add_argument("--plot", action="store_true", help="Generate overview plots with R/ggplot2")
    common.add_argument(
        "-C", "--chr-like-minlen", type=int, default=None,
        help="Minimum contig length (bp) to be considered chromosome-like. Default: 80%% of smallest nuclear ref chromosome.",
    )
    common.add_argument(
        "--add-subgenome-suffix", type=str, default=None,
        help="Add this suffix to contig names when reference has no subgenome identifiers (e.g., 'A' -> chr5A)",
    )
    common.add_argument(
        "--ref-id-pattern",
        action="append",
        default=None,
        help=(
            "Regex for reference IDs with named group 'chrom' and optional 'sg'. "
            "May be provided multiple times; tried in order. "
            "Example: --ref-id-pattern '^(?P<chrom>chr\\d+)(?P<sg>[AB])$' for chr1A/chr1B style. "
            "If unset, presets handle chrN[A-Z], chrN(At|Dt), and chrN patterns."
        ),
    )

    # =========================================================================
    # Pipeline phase toggles
    # =========================================================================
    toggles = p.add_argument_group("Pipeline phase toggles")
    toggles.add_argument(
        "--nt-synteny", action="store_true",
        help="Run nucleotide-based synteny alignments (minimap2/mm2plus) for QA/diagnostics.",
    )
    toggles.add_argument("--skip-organelles", action="store_true", help="Skip organelle detection.")
    toggles.add_argument("--skip-rdna", action="store_true", help="Skip rDNA detection.")
    toggles.add_argument("--skip-contaminants", action="store_true", help="Skip contaminant detection.")

    # =========================================================================
    # External tool paths
    # =========================================================================
    tools = p.add_argument_group("External tool paths")
    tools.add_argument("--gffread", default="gffread", help="gffread executable [gffread]")
    tools.add_argument("--miniprot", default="miniprot", help="miniprot executable [miniprot]")
    tools.add_argument("--miniprot-args", default="", help='Extra args passed to miniprot (e.g. "-G 200k") [none]')

    # =========================================================================
    # Reference inputs for classification
    # =========================================================================
    refs = p.add_argument_group("Reference inputs for classification")
    refs.add_argument(
        "--chrC-ref", type=str, default=None,
        help="Chloroplast reference FASTA. Default: extract chrC from --ref if present.",
    )
    refs.add_argument(
        "--chrM-ref", type=str, default=None,
        help="Mitochondrion reference FASTA. Default: extract chrM from --ref if present.",
    )
    refs.add_argument(
        "--rdna-ref", type=str, default=None,
        help="rDNA reference FASTA, or 'default' for bundled Arabidopsis 45S reference.",
    )
    refs.add_argument(
        "--centrifuger-idx", type=str, default=None,
        help="Path prefix to centrifuger index for contaminant screening.",
    )

    # =========================================================================
    # Thresholds: Chromosome assignment
    # =========================================================================
    chr_thresh = p.add_argument_group("Thresholds: Chromosome assignment")
    chr_thresh.add_argument(
        "--assign-min-frac", type=float, default=0.10,
        help="Min union-bp/contig_len to accept assignment [0.10]",
    )
    chr_thresh.add_argument(
        "--assign-min-ratio", type=float, default=1.25,
        help="Min best/second score ratio for unambiguous assignment [1.25]",
    )
    chr_thresh.add_argument(
        "-P", "--chimera-primary-frac", type=float, default=0.8,
        help="Primary fraction threshold for chimera detection [0.8]",
    )
    chr_thresh.add_argument(
        "-S", "--chimera-secondary-frac", type=float, default=0.2,
        help="Secondary fraction threshold for chimera detection [0.2]",
    )
    chr_thresh.add_argument(
        "--low-ref-span-threshold", type=float, default=0.75,
        help="Fraction of ref chrom that must be spanned by syntenic blocks [0.75]",
    )

    # =========================================================================
    # Thresholds: Synteny block building
    # =========================================================================
    block_thresh = p.add_argument_group("Thresholds: Synteny block building")
    block_thresh.add_argument(
        "--assign-minlen", type=int, default=10_000,
        help="Min alignment length for block building [10000]",
    )
    block_thresh.add_argument(
        "--assign-minmapq", type=int, default=0,
        help="Min MAPQ for block building [0]",
    )
    block_thresh.add_argument(
        "--assign-min-ident", type=float, default=0.0,
        help="Min identity proxy (matches/aln_len) [0.0]",
    )
    block_thresh.add_argument(
        "--assign-tp", choices=["P", "PI", "ALL"], default="PI",
        help="(nt-synteny) Which tp tags: P=primary, PI=primary+supp, ALL=all [PI]",
    )

    # =========================================================================
    # Thresholds: Chain scoring
    # =========================================================================
    chain_thresh = p.add_argument_group("Thresholds: Chain scoring")
    chain_thresh.add_argument(
        "--chain-q-gap", type=int, default=200_000,
        help="Max query gap between blocks in chain [200000]",
    )
    chain_thresh.add_argument(
        "--chain-r-gap", type=int, default=400_000,
        help="Max ref gap between blocks in chain [400000]",
    )
    chain_thresh.add_argument(
        "--chain-diag-slop", type=int, default=150_000,
        help="Max diagonal drift in chain [150000]",
    )
    chain_thresh.add_argument(
        "--assign-chain-min-bp", type=int, default=0,
        help="Min chain union-bp to keep [0]",
    )
    chain_thresh.add_argument(
        "--assign-chain-score", choices=["matches", "qbp_ident", "matches_ident"], default="matches",
        help="How to score each chain [matches]",
    )
    chain_thresh.add_argument(
        "--assign-chain-topk", type=int, default=3,
        help="Sum top-K chain scores per contig/ref [3]",
    )
    chain_thresh.add_argument(
        "--assign-ref-score", choices=["topk", "all"], default="all",
        help="Score per ref: topk or all kept chains [all]",
    )

    # =========================================================================
    # Thresholds: Protein-anchor (miniprot) gates
    # =========================================================================
    prot_thresh = p.add_argument_group("Thresholds: Protein-anchor (miniprot) gates")
    prot_thresh.add_argument(
        "--miniprot-min-genes", type=int, default=3,
        help="Min unique genes for best ref assignment [3]",
    )
    prot_thresh.add_argument(
        "--miniprot-min-segments", type=int, default=5,
        help="Min number of synteny segments [5]",
    )
    prot_thresh.add_argument(
        "--miniprot-min-span-frac", type=float, default=0.20,
        help="Min span fraction of contig [0.20]",
    )
    prot_thresh.add_argument(
        "--miniprot-min-span-bp", type=int, default=50_000,
        help="Min absolute span in bp [50000]",
    )

    # =========================================================================
    # Thresholds: Organelle detection
    # =========================================================================
    org_thresh = p.add_argument_group("Thresholds: Organelle detection")
    org_thresh.add_argument(
        "--organelle-min-cov", type=float, default=0.80,
        help="Min query coverage for organelle candidate [0.80]",
    )
    org_thresh.add_argument(
        "--chrC-len-tolerance", type=float, default=0.05,
        help="Length tolerance for chrC candidates (fraction) [0.05]",
    )
    org_thresh.add_argument(
        "--chrM-len-tolerance", type=float, default=0.20,
        help="Length tolerance for chrM candidates (fraction) [0.20]",
    )

    # =========================================================================
    # Thresholds: rDNA detection
    # =========================================================================
    rdna_thresh = p.add_argument_group("Thresholds: rDNA detection")
    rdna_thresh.add_argument(
        "--rdna-min-cov", type=float, default=0.50,
        help="Min query coverage for rDNA classification [0.50]",
    )

    # =========================================================================
    # Thresholds: Chromosome debris detection
    # =========================================================================
    chr_debris_thresh = p.add_argument_group("Thresholds: Chromosome debris detection")
    chr_debris_thresh.add_argument(
        "--chr-debris-min-cov", type=float, default=0.80,
        help="Min query coverage vs assembled chromosomes [0.80]",
    )
    chr_debris_thresh.add_argument(
        "--chr-debris-min-identity", type=float, default=0.90,
        help="Min alignment identity vs assembled chromosomes [0.90]",
    )

    # =========================================================================
    # Thresholds: Contaminant detection
    # =========================================================================
    contam_thresh = p.add_argument_group("Thresholds: Contaminant detection")
    contam_thresh.add_argument(
        "--contaminant-min-score", type=int, default=150,
        help="Min centrifuger score for contaminant classification [150]",
    )

    # =========================================================================
    # Thresholds: Reference-based debris detection
    # =========================================================================
    debris_thresh = p.add_argument_group("Thresholds: Reference-based debris detection")
    debris_thresh.add_argument(
        "--debris-min-cov", type=float, default=0.50,
        help="Min alignment coverage vs reference for debris [0.50]",
    )
    debris_thresh.add_argument(
        "--debris-min-protein-hits", type=int, default=2,
        help="Min miniprot protein hits for debris [2]",
    )

    # =========================================================================
    # Minimap2 tuning (for --nt-synteny mode)
    # =========================================================================
    mm2_tune = p.add_argument_group("Minimap2 tuning (for --nt-synteny mode)")
    mm2_tune.add_argument("-x", "--preset", default="asm20", help="minimap2 preset [asm20]")
    mm2_tune.add_argument("-k", "--kmer", type=int, default=None, help="minimap2 k-mer size")
    mm2_tune.add_argument("-w", "--window", type=int, default=None, help="minimap2 minimizer window size")
    mm2_tune.add_argument(
        "-L", "--aln-minlen", type=int, default=10_000,
        help="Min alignment length in PRIMARY-only QA stats [10000]",
    )

    args = p.parse_args()

    if args.ref_id_pattern:
        set_ref_id_patterns(compile_ref_id_patterns(args.ref_id_pattern))

    # --- Phase 1: Reading reference and query FASTA ---
    print("[info] Phase 1: Reading reference and query FASTA", file=sys.stderr)
    ref = Path(args.ref)
    qry = Path(args.query)
    outprefix = Path(args.outprefix)
    ref_lengths_norm, ref_orig_to_norm, ref_norm_to_orig = read_fasta_lengths_with_map(ref)
    qry_lengths = read_fasta_lengths(qry)
    print(f"[info] Query contigs: {len(qry_lengths)}", file=sys.stderr)

    # Outputs (shared)
    ref_lengths_tsv = Path(str(outprefix) + ".ref_lengths.tsv")
    segments_tsv = Path(str(outprefix) + ".segments.tsv")
    chain_summary_tsv = Path(str(outprefix) + ".evidence_summary.tsv")
    macro_blocks_tsv = Path(str(outprefix) + ".macro_blocks.tsv")

    # minimap2 QA outputs
    paf_gz = Path(str(outprefix) + ".all_vs_ref.paf.gz")
    err_log = Path(str(outprefix) + ".alignment.err")
    stats_tsv = Path(str(outprefix) + ".contig_ref.stats.tsv")

    # miniprot outputs
    proteins_faa = Path(str(outprefix) + ".ref_proteins.fa")
    gffread_err = Path(str(outprefix) + ".gffread.err")
    miniprot_paf_gz = Path(str(outprefix) + ".miniprot.paf.gz")
    miniprot_err = Path(str(outprefix) + ".miniprot.err")

    if not ref_lengths_tsv.exists():
        print(f"[info] Writing reference lengths: {ref_lengths_tsv}", file=sys.stderr)
        write_ref_lengths_tsv(ref_lengths_tsv, ref)
    else:
        print(f"[info] Reference lengths TSV exists, reusing: {ref_lengths_tsv}", file=sys.stderr)

    # --- Phase 2: Nucleotide synteny analysis (optional QA) ---
    print("[info] Phase 2: Nucleotide synteny analysis (optional QA)", file=sys.stderr)
    if args.nt_synteny:
        run_minimap2_synteny(ref, qry, paf_gz, args.threads, args.preset, args.kmer, args.window, err_log)
        (
            aln_sum_p,
            match_sum_p,
            _alnlen_sum_p,
            _ctotal_p,
            _crefs_p,
            _cqlen_p,
            _best_aln_p,
            _best_ref_p,
            _second_aln_p,
            _second_ref_p,
        ) = parse_paf_primary(paf_gz, args.aln_minlen)
        write_stats_tsv(stats_tsv, aln_sum_p, match_sum_p)
        print(f"[done] Stats (primary QA): {stats_tsv}", file=sys.stderr)
    else:
        print("[info] --nt-synteny not set: skipping nucleotide alignments + primary-only QA stats.", file=sys.stderr)

    # --- Phase 3: Protein-anchor synteny analysis ---
    print("[info] Phase 3: Protein-anchor synteny analysis", file=sys.stderr)
    ref_gff3 = Path(args.ref_gff3)
    print(f"[info] Protein-anchor mode enabled (ref GFF3: {ref_gff3})", file=sys.stderr)

    tx2loc, tx2gene = parse_gff3_transcript_coords(ref_gff3, ref_id_map=ref_orig_to_norm)
    if not (tx2loc and tx2gene):
        raise RuntimeError(f"No transcript features found in GFF3 (mRNA/transcript with ID=...). File: {ref_gff3}")

    run_gffread_extract_proteins(args.gffread, ref, ref_gff3, proteins_faa, gffread_err)
    run_miniprot(args.miniprot, qry, proteins_faa, miniprot_paf_gz, args.threads, args.miniprot_args, miniprot_err)

    ev = parse_miniprot_synteny_evidence_and_segments(
        miniprot_paf_gz=miniprot_paf_gz,
        contig_lengths=qry_lengths,
        tx2loc=tx2loc,
        tx2gene=tx2gene,
        assign_minlen=args.assign_minlen,
        assign_minmapq=args.assign_minmapq,
        chain_q_gap=args.chain_q_gap,
        chain_r_gap=args.chain_r_gap,
        chain_diag_slop=args.chain_diag_slop,
        assign_min_ident=args.assign_min_ident,
        assign_chain_topk=args.assign_chain_topk,
        assign_chain_score=args.assign_chain_score,
        assign_chain_min_bp=args.assign_chain_min_bp,
        assign_ref_score=args.assign_ref_score,
    )

    # --- Protein-anchor anti-contamination gate ---
    seg_count, span_bp = build_segment_support_from_rows(ev.chain_segments_rows)

    def _passes_gate(
        q: str,
        ref_id: str,
        clen: int,
        min_segments: int,
        min_span_bp: int,
        min_span_frac: float,
        min_genes: int,
    ) -> tuple[bool, int, int, float, int]:
        key = (q, ref_id)
        nseg = int(seg_count.get(key, 0) or 0)
        spbp = int(span_bp.get(key, 0) or 0)
        spfrac = (spbp / clen) if clen > 0 else 0.0
        ngen = int(ev.qr_gene_count.get(key, 0) or 0)

        ok = (nseg >= min_segments) and (spbp >= min_span_bp) and (spfrac >= min_span_frac) and (ngen >= min_genes)
        return ok, nseg, spbp, spfrac, ngen

    # Gate-aware reranking over candidate refs
    # Make mutable copies for reranking
    best_ref = dict(ev.best_ref)
    best_score = dict(ev.best_score)
    best_bp = dict(ev.best_bp)
    second_ref = dict(ev.second_ref)
    second_score = dict(ev.second_score)
    second_bp = dict(ev.second_bp)

    n_demoted = 0
    n_switched = 0
    qr_ref_score = ev.qr_weight_all if args.assign_ref_score == "all" else ev.qr_score_topk

    candidates_by_contig = defaultdict(list)
    for (q, ref_id), sc in qr_ref_score.items():
        if sc > 0:
            candidates_by_contig[q].append((float(sc), ref_id))

    print(f"[info] Assignment score per ref: {args.assign_ref_score}", file=sys.stderr)

    for q in qry_lengths.keys():
        clen = int(qry_lengths.get(q, 0) or ev.qlens_from_paf.get(q, 0) or 0)
        if clen <= 0:
            continue

        cands = candidates_by_contig.get(q, [])
        if not cands:
            best_ref[q] = ""
            best_score[q] = 0.0
            best_bp[q] = 0
            second_ref[q] = ""
            second_score[q] = 0.0
            second_bp[q] = 0
            continue

        cands.sort(key=lambda x: (x[0], int(ev.qr_union_bp.get((q, x[1]), 0) or 0)), reverse=True)
        orig_best = cands[0][1]

        passing = []
        for sc, ref_id in cands:
            ok, _nseg, _spbp, _spfrac, _ngen = _passes_gate(
                q,
                ref_id,
                clen,
                args.miniprot_min_segments,
                args.miniprot_min_span_bp,
                args.miniprot_min_span_frac,
                args.miniprot_min_genes,
            )
            if ok:
                passing.append((sc, ref_id))
            if len(passing) >= 2:
                break

        if not passing:
            best_ref[q] = ""
            best_score[q] = 0.0
            best_bp[q] = 0
            second_ref[q] = ""
            second_score[q] = 0.0
            second_bp[q] = 0
            n_demoted += 1
            continue

        chosen_score, chosen_ref_id = passing[0]
        best_ref[q] = chosen_ref_id
        best_score[q] = float(chosen_score)
        best_bp[q] = int(ev.qr_union_bp.get((q, chosen_ref_id), 0) or 0)

        if len(passing) > 1:
            second_score_val, second_ref_id = passing[1]
            second_ref[q] = second_ref_id
            second_score[q] = float(second_score_val)
            second_bp[q] = int(ev.qr_union_bp.get((q, second_ref_id), 0) or 0)
        else:
            second_ref[q] = ""
            second_score[q] = 0.0
            second_bp[q] = 0

        if chosen_ref_id != orig_best:
            n_switched += 1

    print(f"[info] Protein-anchor gate-aware rerank: demoted={n_demoted} switched={n_switched}", file=sys.stderr)

    plot_suffix = "synteny (miniprot)"
    print(f"[done] Proteins:         {proteins_faa}", file=sys.stderr)
    print(f"[done] miniprot PAF.gz:  {miniprot_paf_gz}", file=sys.stderr)
    print(f"[done] miniprot err:     {miniprot_err}", file=sys.stderr)

    # Write segments + evidence summary TSVs
    write_chain_segments_tsv(segments_tsv, ev.chain_segments_rows, ref_norm_to_orig=ref_norm_to_orig)
    write_chain_summary_tsv(chain_summary_tsv, ev.chain_summary_rows, ref_norm_to_orig=ref_norm_to_orig)
    print(f"[done] Segments TSV:      {segments_tsv}", file=sys.stderr)
    print(f"[done] Evidence TSV:      {chain_summary_tsv}", file=sys.stderr)

    # Write macro blocks TSV
    write_macro_blocks_tsv(macro_blocks_tsv, ev.macro_block_rows, ref_norm_to_orig=ref_norm_to_orig)
    print(f"[done] Macro blocks TSV:  {macro_blocks_tsv}", file=sys.stderr)

    # --- Compute chr_like_minlen if not provided ---
    ref_lengths = ref_lengths_norm
    if args.chr_like_minlen is None:
        min_nuclear = get_min_nuclear_chrom_length(ref_lengths)
        chr_like_minlen = int(min_nuclear * 0.8) if min_nuclear > 0 else 1_000_000
        print(f"[info] chr_like_minlen not specified, using 80% of smallest nuclear chrom: {chr_like_minlen}", file=sys.stderr)
    else:
        chr_like_minlen = args.chr_like_minlen

    # --- Identify chromosome contigs early (before BLAST phases) ---
    # This allows creating a filtered FASTA for faster organelle/rDNA/contaminant detection
    chromosome_contigs: Set[str] = set()
    for contig, ref_id in best_ref.items():
        if ref_id and qry_lengths.get(contig, 0) >= chr_like_minlen:
            chromosome_contigs.add(contig)

    # Create work directory for classification outputs
    work_dir = outprefix.parent / f"{outprefix.name}_classification"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Create filtered FASTA with only non-chromosome contigs for BLAST-based detection
    # This significantly speeds up organelle/rDNA/contaminant BLAST searches
    non_chrom_contigs = set(qry_lengths.keys()) - chromosome_contigs
    non_chrom_fasta = work_dir / "non_chromosome_contigs.fa"
    if non_chrom_contigs:
        print(f"[info] Creating filtered FASTA for BLAST ({len(non_chrom_contigs)} non-chromosome contigs)", file=sys.stderr)
        write_filtered_fasta(qry, non_chrom_fasta, non_chrom_contigs)
        blast_query_fasta = non_chrom_fasta
        blast_query_lengths = {k: v for k, v in qry_lengths.items() if k in non_chrom_contigs}
    else:
        # Edge case: all contigs are chromosomes
        blast_query_fasta = qry
        blast_query_lengths = qry_lengths

    # --- Phase 4: Organelle detection ---
    chrC_contig: Optional[str] = None
    chrM_contig: Optional[str] = None
    organelle_debris: Set[str] = set()

    if not args.skip_organelles:
        print("[info] Phase 4: Organelle detection", file=sys.stderr)
        chrC_ref, chrM_ref = prepare_organelle_references(
            ref_fasta=ref,
            chrC_ref_arg=getattr(args, 'chrC_ref', None),
            chrM_ref_arg=getattr(args, 'chrM_ref', None),
            work_dir=work_dir / "organelles",
            ref_norm_to_orig=ref_norm_to_orig,
        )

        if chrC_ref or chrM_ref:
            chrC_contig, chrM_contig, organelle_debris = detect_organelles(
                query_fasta=blast_query_fasta,
                query_lengths=blast_query_lengths,
                chrC_ref=chrC_ref,
                chrM_ref=chrM_ref,
                work_dir=work_dir / "organelles",
                threads=args.threads,
                min_coverage=args.organelle_min_cov,
                chrC_len_tol=getattr(args, 'chrC_len_tolerance', 0.05),
                chrM_len_tol=getattr(args, 'chrM_len_tolerance', 0.20),
            )
    else:
        print("[info] Phase 4: Skipping organelle detection (--skip-organelles)", file=sys.stderr)

    # --- Phase 5: rDNA detection ---
    rdna_contigs: Set[str] = set()
    already_classified = {chrC_contig, chrM_contig} | organelle_debris
    already_classified.discard(None)

    if not args.skip_rdna:
        print("[info] Phase 5: rDNA detection", file=sys.stderr)
        script_dir = Path(__file__).resolve().parent
        rdna_ref = prepare_rdna_reference(args.rdna_ref, script_dir)

        if rdna_ref:
            rdna_contigs = detect_rdna_contigs(
                query_fasta=blast_query_fasta,
                query_lengths=blast_query_lengths,
                rdna_ref=rdna_ref,
                work_dir=work_dir / "rdna",
                threads=args.threads,
                min_coverage=args.rdna_min_cov,
                exclude_contigs=already_classified,
            )
    else:
        print("[info] Phase 5: Skipping rDNA detection (--skip-rdna)", file=sys.stderr)

    # --- Phase 6: Chromosome debris detection ---
    # TWO-PRONGED DEBRIS DETECTION STRATEGY:
    # 1. Chromosome debris (this phase): Aligns candidates against assembled chromosome
    #    contigs to find near-identical duplicates (haplotigs, bubble branches, assembly
    #    artifacts). Uses high thresholds (80% coverage, 90% identity) to catch sequences
    #    that are essentially copies of already-assembled regions.
    #
    # 2. Reference-based debris (Phase 8): Aligns remaining candidates against the
    #    reference genome and checks for protein hits. Catches divergent sequences with
    #    ancestral homology that wouldn't match the assembled chromosomes well but still
    #    have detectable reference similarity.
    #
    # Both methods are needed because:
    # - Chromosome debris catches assembly-specific artifacts (not in reference)
    # - Reference-based catches divergent haplotigs and non-coding homologs
    # - Neither alone covers all cases
    print("[info] Phase 6: Chromosome debris detection", file=sys.stderr)
    already_classified = already_classified | rdna_contigs

    chromosome_debris: Set[str] = set()
    if chromosome_contigs:
        chromosome_debris = detect_chromosome_debris(
            query_fasta=qry,
            query_lengths=qry_lengths,
            chromosome_contigs=chromosome_contigs,
            work_dir=work_dir / "chr_debris",
            threads=args.threads,
            min_coverage=args.chr_debris_min_cov,
            min_identity=args.chr_debris_min_identity,
            exclude_contigs=already_classified,
        )
        already_classified = already_classified | chromosome_debris

    # --- Phase 7: Contaminant detection ---
    contaminants: Dict[str, Tuple[int, str]] = {}

    if not args.skip_contaminants and args.centrifuger_idx:
        print("[info] Phase 7: Contaminant detection", file=sys.stderr)
        # Create residual FASTA excluding all debris for faster contaminant screening
        residual_contigs = set(qry_lengths.keys()) - already_classified - chromosome_contigs
        if residual_contigs:
            residual_fasta = work_dir / "residual_for_contaminant_screen.fa"
            write_filtered_fasta(qry, residual_fasta, residual_contigs)
            residual_lengths = {k: v for k, v in qry_lengths.items() if k in residual_contigs}

            contaminants = detect_contaminants(
                query_fasta=residual_fasta,
                query_lengths=residual_lengths,
                centrifuger_idx=args.centrifuger_idx,
                work_dir=work_dir / "contaminants",
                threads=args.threads,
                min_score=args.contaminant_min_score,
                exclude_contigs=set(),  # already filtered
            )
    else:
        print("[info] Phase 7: Skipping contaminant detection", file=sys.stderr)

    # --- Phase 8: Debris/unclassified classification ---
    print("[info] Phase 8: Debris/unclassified classification", file=sys.stderr)

    # Remaining contigs to classify (after all other classifications)
    already_classified = already_classified | set(contaminants.keys()) | chromosome_contigs
    remaining_contigs = set(qry_lengths.keys()) - already_classified

    # Also check remaining contigs for reference/protein-based debris
    additional_debris, _unclassified = classify_debris_and_unclassified(
        remaining_contigs=remaining_contigs,
        query_fasta=qry,
        query_lengths=qry_lengths,
        ref_fasta=ref,
        chrs_fasta=None,  # Will be created later
        qr_gene_count=ev.qr_gene_count,
        work_dir=work_dir / "debris",
        threads=args.threads,
        min_coverage=args.debris_min_cov,
        min_protein_hits=args.debris_min_protein_hits,
    )
    other_debris = additional_debris  # debris from ref/protein alignment (not chromosome-vs-chromosome)

    # --- Phase 9: Gene count statistics ---
    print("[info] Phase 9: Gene count statistics", file=sys.stderr)
    ref_gene_counts = count_genes_per_ref_chrom(ref_gff3, ref_id_map=ref_orig_to_norm)
    compute_mean_genes_per_Mbp(
        qr_gene_count=ev.qr_gene_count,
        query_lengths=qry_lengths,
        chromosome_contigs=chromosome_contigs,
        best_ref=best_ref,
    )

    # --- Phase 10: Orientation determination ---
    print("[info] Phase 10: Orientation determination", file=sys.stderr)
    contig_orientations = determine_contig_orientations(
        macro_block_rows=ev.macro_block_rows,
        best_ref=best_ref,
        chromosome_contigs=chromosome_contigs,
    )

    # --- Phase 11: Classify all contigs ---
    print("[info] Phase 11: Classifying all contigs", file=sys.stderr)
    classifications = classify_all_contigs(
        query_fasta=qry,
        query_lengths=qry_lengths,
        best_ref=best_ref,
        chr_like_minlen=chr_like_minlen,
        ev=ev,
        ref_gene_counts=ref_gene_counts,
        chrC_contig=chrC_contig,
        chrM_contig=chrM_contig,
        organelle_debris=organelle_debris,
        rdna_contigs=rdna_contigs,
        contaminants=contaminants,
        chromosome_debris=chromosome_debris,
        other_debris=other_debris,
        add_subgenome_suffix=args.add_subgenome_suffix,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    # Update orientation info in classifications
    for clf in classifications:
        clf.reversed = contig_orientations.get(clf.original_name, False)

    # --- Phase 12: Write FASTA outputs ---
    print("[info] Phase 12: Writing FASTA outputs", file=sys.stderr)
    write_classified_fastas(
        query_fasta=qry,
        classifications=classifications,
        contig_orientations=contig_orientations,
        output_prefix=outprefix,
    )

    # --- Phase 13: Write enhanced summary TSV ---
    print("[info] Phase 13: Writing enhanced summary TSV", file=sys.stderr)
    summary_tsv = Path(str(outprefix) + ".contig_summary.tsv")

    write_contig_summary_tsv(
        output_path=summary_tsv,
        classifications=classifications,
        contig_orientations=contig_orientations,
        all_contig_lengths=qry_lengths,
        contig_total=ev.contig_total,
        contig_refs=ev.contig_refs,
        qlens_from_paf=ev.qlens_from_paf,
        qr_union_bp=ev.qr_union_bp,
        qr_matches=ev.qr_matches,
        qr_alnlen=ev.qr_alnlen,
        qr_gene_count=ev.qr_gene_count,
        best_ref=best_ref,
        best_score=best_score,
        best_bp=best_bp,
        second_ref=second_ref,
        second_score=second_score,
        second_bp=second_bp,
        chr_like_minlen=chr_like_minlen,
        chimera_primary_frac=args.chimera_primary_frac,
        chimera_secondary_frac=args.chimera_secondary_frac,
        assign_min_frac=args.assign_min_frac,
        assign_min_ratio=args.assign_min_ratio,
        ref_norm_to_orig=ref_norm_to_orig,
    )

    print(f"[done] Summary:           {summary_tsv}", file=sys.stderr)
    print(f"[done] Ref lengths:       {ref_lengths_tsv}", file=sys.stderr)

    # Print classification summary
    clf_counts: Dict[str, int] = defaultdict(int)
    for clf in classifications:
        clf_counts[clf.classification] += 1
    print("[done] Classification summary:", file=sys.stderr)
    for cat, count in sorted(clf_counts.items()):
        print(f"       {cat}: {count}", file=sys.stderr)

    if args.plot:
        run_plot(
            summary_tsv,
            ref_lengths_tsv,
            segments_tsv,
            chain_summary_tsv,
            macro_blocks_tsv,
            outprefix,
            chr_like_minlen,
            plot_suffix,
        )


if __name__ == "__main__":
    main()
