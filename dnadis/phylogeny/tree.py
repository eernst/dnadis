"""Outgroup resolution and IQ-TREE invocation for the species tree.

The 'auto' outgroup heuristic picks the taxon with the lowest mean
identity to the rest of the supermatrix.  That is not a biologically
conclusive way to root a species tree, so we log a clear warning when
the user requests it.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from dnadis.alignment.external_tools import run_iqtree
from dnadis.utils.logging_config import get_logger

logger = get_logger("phylogeny")


def _read_aligned_fasta_ordered(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur_name: Optional[str] = None
    cur: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_name is not None:
                    seqs[cur_name] = "".join(cur)
                cur_name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.strip())
        if cur_name is not None:
            seqs[cur_name] = "".join(cur)
    return seqs


def _pair_identity(a: str, b: str) -> float:
    """Fraction of aligned columns where a and b agree, ignoring gap-gap columns."""
    if len(a) != len(b):
        return 0.0
    matches = 0
    informative = 0
    for ca, cb in zip(a, b):
        if ca == "-" and cb == "-":
            continue
        informative += 1
        if ca == cb and ca != "-":
            matches += 1
    if informative == 0:
        return 0.0
    return matches / informative


def auto_outgroup(supermatrix_fasta: Path) -> Optional[str]:
    """Pick the taxon with the lowest mean pairwise identity vs the others.

    Returns ``None`` if the supermatrix has fewer than three taxa
    (a two-taxon tree has no informative outgroup to choose).
    """
    seqs = _read_aligned_fasta_ordered(supermatrix_fasta)
    names = list(seqs.keys())
    if len(names) < 3:
        return None
    mean_ids: Dict[str, float] = {}
    for i, ni in enumerate(names):
        ids = []
        for j, nj in enumerate(names):
            if i == j:
                continue
            ids.append(_pair_identity(seqs[ni], seqs[nj]))
        mean_ids[ni] = sum(ids) / len(ids) if ids else 0.0
    return min(mean_ids, key=lambda n: mean_ids[n])


def resolve_outgroup(
    requested: str,
    leaf_labels: List[str],
    reference_labels,
    supermatrix_fasta: Path,
) -> Optional[str]:
    """Map a user-facing ``--phylo-outgroup`` value to a leaf label.

    ``reference_labels`` may be a single string (legacy), a list of
    strings (polyploid reference, one leaf per subgenome), or ``None``.
    When ``--phylo-outgroup=reference`` matches more than one reference
    leaf, all present reference leaves are joined with commas — IQ-TREE
    accepts that form for an outgroup clade.

    Returns ``None`` when the tree should be unrooted (requested == 'none'),
    or when the requested taxon is not present.
    """
    if reference_labels is None:
        ref_list: List[str] = []
    elif isinstance(reference_labels, str):
        ref_list = [reference_labels]
    else:
        ref_list = list(reference_labels)

    req = (requested or "none").strip().lower()
    if req == "none":
        return None
    if req == "reference":
        if not ref_list:
            logger.warning(
                "--phylo-outgroup=reference requested but reference is not in the tree; "
                "leaving tree unrooted"
            )
            return None
        present = [lbl for lbl in ref_list if lbl in leaf_labels]
        if not present:
            logger.warning(
                f"Reference leaves {ref_list} not in tree; leaving unrooted"
            )
            return None
        if len(present) == 1:
            return present[0]
        logger.info(
            f"--phylo-outgroup=reference: polyploid reference; using all "
            f"{len(present)} reference subgenomes as the outgroup clade: {present}"
        )
        return ",".join(present)
    if req == "auto":
        logger.warning(
            "--phylo-outgroup=auto picks the most-divergent taxon by alignment identity; "
            "this is NOT a biologically conclusive way to root a species tree."
        )
        return auto_outgroup(supermatrix_fasta)
    # Fall-through: requested matches a specific leaf label (case-sensitive)
    if requested in leaf_labels:
        return requested
    # Try matching as an assembly name to any of its subgenome leaves
    matches = [lbl for lbl in leaf_labels if lbl == requested or lbl.startswith(f"{requested}_")]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        logger.warning(
            f"--phylo-outgroup={requested!r} is ambiguous (matches {matches}); "
            f"specify the exact leaf label (e.g. {matches[0]}); leaving unrooted"
        )
        return None
    logger.warning(
        f"--phylo-outgroup={requested!r} does not match any leaf in the tree "
        f"(leaves: {leaf_labels}); leaving unrooted"
    )
    return None


def build_tree(
    supermatrix_fasta: Path,
    out_prefix: Path,
    threads: int,
    max_mem_gb: int,
    bootstrap: int,
    alrt: int,
    models: str,
    outgroup_label: Optional[str],
    err_path: Optional[Path] = None,
) -> Optional[Path]:
    """Run IQ-TREE and return the path to the ``.treefile`` on success."""
    ok = run_iqtree(
        supermatrix=supermatrix_fasta,
        prefix=out_prefix,
        threads=threads,
        max_mem_gb=max_mem_gb,
        bootstrap=bootstrap,
        alrt=alrt,
        models=models,
        outgroup=outgroup_label,
        err_path=err_path,
    )
    if not ok:
        return None
    treefile = Path(f"{out_prefix}.treefile")
    return treefile if treefile.exists() else None
