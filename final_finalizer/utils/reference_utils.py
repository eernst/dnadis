#!/usr/bin/env python3
"""
Reference utilities for final_finalizer.

Contains functions for reference ID normalization, GFF3 parsing, and PAF helpers.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from final_finalizer.utils.io_utils import open_maybe_gzip
from final_finalizer.utils.sequence_utils import read_fasta_lengths


# ----------------------------
# Reference ID normalization
# ----------------------------
def normalize_ref_id(ref_id: str) -> str:
    """Normalize reference IDs to a lowercase 'chr' prefix when present.

    Handles various case patterns: Chr1 -> chr1, CHR1 -> chr1, etc.
    IDs already starting with lowercase 'chr' are unchanged.
    """
    if len(ref_id) >= 3 and ref_id[:3].upper() == "CHR" and not ref_id.startswith("chr"):
        return "chr" + ref_id[3:]
    return ref_id


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


def _is_nuclear_chromosome(ref_id: str) -> bool:
    """Check if a reference ID matches a nuclear chromosome pattern.

    Returns True only for IDs that positively match chromosome patterns
    (e.g., chr1, chr5A, Chr12Dt). Returns False for organelles, scaffolds,
    unplaced contigs, or any other non-chromosome sequences.
    """
    if normalize_organelle_id(ref_id) is not None:
        return False

    norm_id = normalize_ref_id(ref_id)
    for pat in get_ref_id_patterns():
        if pat.match(norm_id):
            return True
    return False


def get_min_nuclear_chrom_length(ref_lengths: Dict[str, int]) -> int:
    """Return the length of the smallest nuclear chromosome.

    Only considers reference sequences that positively match chromosome
    naming patterns (e.g., chr1, chr5A). Excludes organelles, scaffolds,
    unplaced contigs, and other non-chromosome sequences.

    Returns 0 if no nuclear chromosomes found.
    """
    nuclear_lengths = [
        length for ref_id, length in ref_lengths.items()
        if _is_nuclear_chromosome(ref_id)
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


# ----------------------------
# GFF3 utilities
# ----------------------------
def filter_gff3_by_ref(
    gff3_path: Path,
    ref_ids: Set[str],
    output_path: Path,
    ref_norm_to_orig: Optional[Dict[str, str]] = None,
) -> Path:
    """Write a filtered GFF3 containing only records whose seqid is in ref_ids.

    Comment/pragma lines are preserved. Returns the output path.
    """
    if output_path.exists():
        return output_path

    with open_maybe_gzip(gff3_path, "rt") as fh_in, output_path.open("w") as fh_out:
        for line in fh_in:
            if not line or line.startswith("#"):
                fh_out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid = parts[0]
            if seqid in ref_ids:
                fh_out.write(line)
                continue

            norm_id = normalize_ref_id(seqid)
            organelle = normalize_organelle_id(seqid)
            if organelle:
                norm_id = organelle

            if ref_norm_to_orig and norm_id in ref_norm_to_orig:
                parts[0] = ref_norm_to_orig[norm_id]
                fh_out.write("\t".join(parts) + "\n")

    return output_path


def write_ref_lengths_tsv(ref_lengths_path: Path, ref_fasta: Path) -> None:
    lengths = read_fasta_lengths(ref_fasta)
    with ref_lengths_path.open("w") as out:
        out.write("ref_id\tchrom_id\tsubgenome\tref_len\n")
        for ref_id, L in sorted(lengths.items()):
            chrom_id, sub = split_chrom_subgenome(ref_id)
            out.write(f"{ref_id}\t{chrom_id}\t{sub}\t{L}\n")


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
