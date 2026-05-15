"""Synthetic genome rearrangements for testing chromosome assignment.

Applies a random sequence of structural operations to a reference FASTA
and records ground-truth metadata for each operation.  The permuted
genome can then be fed to the dnadis pipeline; the detected rearrangement
calls in the output TSV are compared against the recorded operations.

Operations supported:
* inversion
* reciprocal_translocation
* whole_arm_translocation
* fusion
* fission
"""
from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class Rearrangement:
    """Ground-truth metadata for one synthetic rearrangement operation.

    Coordinates use the REFERENCE-side numbering (i.e. positions on the
    original chromosomes before the operation was applied), which is what
    dnadis reports in its detected calls.
    """
    rearrangement_type: str  # 'inversion', 'reciprocal_translocation',
                             # 'whole_arm_translocation', 'fusion', 'fission'
    primary_ref: str
    partner_ref: Optional[str]
    ref_breakpoint: int
    partner_breakpoint: Optional[int]
    span_bp: int
    # Names of the query contigs that carry the rearranged material in the
    # permuted assembly.  For the simple operations below this is just one
    # or two contig names; tracking it makes downstream validation easier.
    query_contigs: List[str] = field(default_factory=list)


def load_genome(path: Path) -> Dict[str, str]:
    """Load a FASTA into an ordered {seqname: sequence} dict."""
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    genome: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    genome[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        genome[name] = "".join(chunks)
    return genome


def write_genome(genome: Dict[str, str], path: Path, line_width: int = 80) -> None:
    """Write a {seqname: sequence} dict to a FASTA file."""
    with path.open("w") as fh:
        for name, seq in genome.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width])
                fh.write("\n")


# ---------------------------------------------------------------------------
# Operations.  Each mutates `genome` in place and returns a Rearrangement
# (or list of them, for reciprocal translocations).
# ---------------------------------------------------------------------------
def apply_inversion(
    genome: Dict[str, str], chrom: str, start: int, end: int,
) -> Rearrangement:
    seq = genome[chrom]
    if start < 0 or end > len(seq) or start >= end:
        raise ValueError(f"inversion bounds out of range: {chrom}[{start}:{end}]")
    inverted = revcomp(seq[start:end])
    genome[chrom] = seq[:start] + inverted + seq[end:]
    return Rearrangement(
        rearrangement_type="inversion",
        primary_ref=chrom,
        partner_ref=None,
        ref_breakpoint=start,
        partner_breakpoint=None,
        span_bp=end - start,
        query_contigs=[chrom],
    )


def apply_reciprocal_translocation(
    genome: Dict[str, str], chrom1: str, brk1: int, chrom2: str, brk2: int,
) -> List[Rearrangement]:
    """Swap arms after the breakpoints between chrom1 and chrom2.

    Post-condition (in REF coords):
      chrom1 carries chrom1_orig[:brk1] + chrom2_orig[brk2:]
      chrom2 carries chrom2_orig[:brk2] + chrom1_orig[brk1:]
    """
    seq1, seq2 = genome[chrom1], genome[chrom2]
    if brk1 <= 0 or brk1 >= len(seq1) or brk2 <= 0 or brk2 >= len(seq2):
        raise ValueError("reciprocal-translocation breakpoints must be interior")
    new1 = seq1[:brk1] + seq2[brk2:]
    new2 = seq2[:brk2] + seq1[brk1:]
    genome[chrom1] = new1
    genome[chrom2] = new2
    return [
        Rearrangement(
            rearrangement_type="reciprocal_translocation",
            primary_ref=chrom1, partner_ref=chrom2,
            ref_breakpoint=brk1, partner_breakpoint=brk2,
            span_bp=(len(seq1) - brk1) + (len(seq2) - brk2),
            query_contigs=[chrom1, chrom2],
        ),
    ]


def apply_whole_arm_translocation(
    genome: Dict[str, str], donor: str, donor_brk: int, acceptor: str,
) -> Rearrangement:
    """Move donor_orig[donor_brk:] onto the end of acceptor.

    Post-condition (in REF coords):
      donor carries donor_orig[:donor_brk]   (loses its right arm)
      acceptor carries acceptor_orig + donor_orig[donor_brk:]
    """
    donor_seq = genome[donor]
    acceptor_seq = genome[acceptor]
    if donor_brk <= 0 or donor_brk >= len(donor_seq):
        raise ValueError("whole-arm translocation breakpoint must be interior")
    moved = donor_seq[donor_brk:]
    genome[donor] = donor_seq[:donor_brk]
    genome[acceptor] = acceptor_seq + moved
    return Rearrangement(
        rearrangement_type="whole_arm_translocation",
        primary_ref=donor, partner_ref=acceptor,
        ref_breakpoint=donor_brk, partner_breakpoint=len(acceptor_seq),
        span_bp=len(moved),
        query_contigs=[donor, acceptor],
    )


def apply_fusion(
    genome: Dict[str, str], chrom1: str, chrom2: str, new_name: Optional[str] = None,
) -> Rearrangement:
    """Concatenate chrom2 onto chrom1; both originals removed from the dict."""
    seq1, seq2 = genome[chrom1], genome[chrom2]
    fused_name = new_name or f"fus_{chrom1}_{chrom2}"
    fused_seq = seq1 + seq2
    # Build a new ordered dict that preserves chromosome order, with the
    # fused chromosome placed where chrom1 was and chrom2 removed.
    new_genome: Dict[str, str] = {}
    for name, seq in genome.items():
        if name == chrom1:
            new_genome[fused_name] = fused_seq
        elif name == chrom2:
            continue
        else:
            new_genome[name] = seq
    genome.clear()
    genome.update(new_genome)
    return Rearrangement(
        rearrangement_type="fusion",
        primary_ref=chrom1, partner_ref=chrom2,
        ref_breakpoint=len(seq1),  # junction position (within fused contig)
        partner_breakpoint=0,
        span_bp=len(fused_seq),
        query_contigs=[fused_name],
    )


def apply_fission(
    genome: Dict[str, str], chrom: str, brk: int,
    name_a: Optional[str] = None, name_b: Optional[str] = None,
) -> Rearrangement:
    """Split chrom at brk into two new chromosomes; original removed."""
    seq = genome[chrom]
    if brk <= 0 or brk >= len(seq):
        raise ValueError("fission breakpoint must be interior")
    a = name_a or f"fis_{chrom}_a"
    b = name_b or f"fis_{chrom}_b"
    seq_a, seq_b = seq[:brk], seq[brk:]
    new_genome: Dict[str, str] = {}
    for name, s in genome.items():
        if name == chrom:
            new_genome[a] = seq_a
            new_genome[b] = seq_b
        else:
            new_genome[name] = s
    genome.clear()
    genome.update(new_genome)
    return Rearrangement(
        rearrangement_type="fission",
        primary_ref=chrom, partner_ref=None,
        ref_breakpoint=brk, partner_breakpoint=None,
        span_bp=len(seq),
        query_contigs=[a, b],
    )


# ---------------------------------------------------------------------------
# Random generator
# ---------------------------------------------------------------------------
_ORGANELLE_PATTERNS = ("chrm", "chrc", "chrpt", "chrmt")
OP_TYPES = (
    "inversion",
    "reciprocal_translocation",
    "whole_arm_translocation",
    "fusion",
    "fission",
)


def _is_nuclear(name: str) -> bool:
    n = name.lower()
    return not any(p in n for p in _ORGANELLE_PATTERNS)


def generate_random_assembly(
    reference_path: Path,
    output_path: Path,
    seed: int,
    n_rearrangements: int = 2,
    min_chrom_length: int = 10_000_000,
    min_segment_length: int = 2_000_000,
    op_types: tuple = OP_TYPES,
) -> List[Rearrangement]:
    """Apply a random sequence of rearrangements to a reference FASTA.

    Returns the list of ground-truth Rearrangement records (in application
    order).  The permuted FASTA is written to ``output_path``.

    Chromosomes shorter than ``min_chrom_length`` and any matching common
    organelle patterns (chrM, chrC, etc.) are excluded from the candidate
    pool to keep operations on chromosome-scale sequences.
    """
    rng = random.Random(seed)
    genome = load_genome(reference_path)
    eligible = [c for c, s in genome.items()
                if len(s) >= min_chrom_length and _is_nuclear(c)]
    if not eligible:
        raise ValueError(f"no chromosomes ≥ {min_chrom_length} bp in {reference_path}")

    rearrangements: List[Rearrangement] = []
    for _ in range(n_rearrangements):
        # Refresh live set each iteration since fusion/fission rename
        # chromosomes and reduce/expand the candidate pool.
        live = [c for c in genome
                if _is_nuclear(c) and len(genome[c]) >= min_chrom_length]
        if not live:
            break
        # Pick an op compatible with current live set.
        compatible = list(op_types)
        if len(live) < 2:
            compatible = [op for op in compatible
                          if op not in ("reciprocal_translocation",
                                        "whole_arm_translocation", "fusion")]
        if not compatible:
            break
        op = rng.choice(compatible)

        if op == "inversion":
            chrom = rng.choice(live)
            seq_len = len(genome[chrom])
            lo = min_segment_length
            hi = seq_len - 2 * min_segment_length
            if hi <= lo:
                continue
            start = rng.randint(lo, hi)
            end = rng.randint(start + min_segment_length,
                              seq_len - min_segment_length)
            rearrangements.append(apply_inversion(genome, chrom, start, end))

        elif op == "reciprocal_translocation":
            c1, c2 = rng.sample(live, 2)
            brk1 = rng.randint(min_segment_length, len(genome[c1]) - min_segment_length)
            brk2 = rng.randint(min_segment_length, len(genome[c2]) - min_segment_length)
            rearrangements.extend(
                apply_reciprocal_translocation(genome, c1, brk1, c2, brk2))

        elif op == "whole_arm_translocation":
            donor, acceptor = rng.sample(live, 2)
            brk = rng.randint(min_segment_length,
                              len(genome[donor]) - min_segment_length)
            rearrangements.append(
                apply_whole_arm_translocation(genome, donor, brk, acceptor))

        elif op == "fusion":
            c1, c2 = rng.sample(live, 2)
            rearrangements.append(apply_fusion(genome, c1, c2))

        elif op == "fission":
            chrom = rng.choice(live)
            seq_len = len(genome[chrom])
            brk = rng.randint(min_segment_length, seq_len - min_segment_length)
            rearrangements.append(apply_fission(genome, chrom, brk))

    write_genome(genome, output_path)
    return rearrangements
