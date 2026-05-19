"""Parse compleasm BUSCO outputs and group per-subgenome single-copy hits.

A "leaf" in the species tree is an ``(assembly, subgenome)`` pair.  For each
leaf we determine which BUSCOs are single-copy *within that leaf*, then take
the intersection across all retained leaves to build a list of genes usable
for tree construction.

The reference genome contributes a single ``(reference_name, "ref")`` leaf
(unless ``--phylo-skip-reference`` is set), so allopolyploid assemblies and
the reference are placed on the same footing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from dnadis.models import CompleasmResult, ContigClassification
from dnadis.utils.io_utils import file_exists_and_valid
from dnadis.utils.logging_config import get_logger
from dnadis.utils.reference_utils import split_chrom_subgenome

logger = get_logger("phylogeny")


# Compleasm assigns a *global* status to each BUSCO based on how many hits
# were found across the whole assembly: Complete = 1 hit, Duplicated = ≥2
# hits, Fragmented = partial, Missing = none.  For per-leaf single-copy
# counting we ignore the global label entirely (a BUSCO duplicated globally
# can still be single-copy *within* a subgenome) and instead count how many
# rows land on contigs belonging to each leaf.  Fragmented hits are excluded
# because compleasm's translated_protein.fasta does not contain reliable
# protein sequence for them.
_INCLUDE_STATUSES = {"Complete", "Duplicated"}


@dataclass(frozen=True)
class LeafId:
    """Identifier for one taxon in the species tree.

    A leaf is an ``(assembly, ref_subgenome, query_subgenome)`` triple:

    * ``ref_subgenome`` is the suffix encoded in the assigned reference
      chromosome name (chr3A → ``A``, chr3P → ``P``).  Polyploid references
      whose subgenomes diverged before the speciation events under study
      contribute multiple leaves per query assembly.
    * ``query_subgenome`` is the within-query duplication group inferred
      by :func:`classification.classifier.infer_query_subgenomes` (often
      ``None`` for the primary haplotype and ``"B"`` for a secondary one).

    Both subgenome fields are ``None`` for the reference leaf and for
    haploid queries on monoploid references.
    """
    assembly: str
    ref_subgenome: Optional[str] = None
    query_subgenome: Optional[str] = None

    @property
    def label(self) -> str:
        parts = [self.assembly]
        if self.ref_subgenome:
            parts.append(self.ref_subgenome)
        if self.query_subgenome:
            parts.append(self.query_subgenome)
        return "_".join(parts)


@dataclass
class BuscoHit:
    """One row in a compleasm full_table_busco_format.tsv file."""
    busco_id: str
    status: str
    contig: Optional[str]
    start: Optional[int]
    end: Optional[int]
    strand: Optional[str]
    score: Optional[float]
    length: Optional[int]


@dataclass
class LeafBuscos:
    """Per-leaf BUSCO presence summary."""
    leaf: LeafId
    hits_by_busco: Dict[str, List[BuscoHit]] = field(default_factory=dict)
    translated_protein_path: Optional[Path] = None
    n_lineage_total: int = 0

    def single_copy_genes(self) -> set:
        """BUSCO IDs with exactly one hit on a contig in this leaf.

        Compleasm's global ``Duplicated`` status is irrelevant here: a BUSCO
        with one hit per reference subgenome is locally single-copy in each
        of those subgenome leaves even though compleasm calls it duplicated
        across the whole assembly.  ``Fragmented``/``Missing`` rows were
        already excluded when this leaf was built.
        """
        return {bid for bid, hits in self.hits_by_busco.items() if len(hits) == 1}

    def completeness_pct(self) -> float:
        """Per-leaf single-copy BUSCO percentage (0-100)."""
        if self.n_lineage_total <= 0:
            return 0.0
        return 100.0 * len(self.single_copy_genes()) / float(self.n_lineage_total)


def parse_full_table(path: Path) -> List[BuscoHit]:
    """Parse a compleasm full_table_busco_format.tsv file.

    Skips comment lines starting with ``#``.  Missing-status rows have only
    two columns (busco_id, status) and produce a hit with ``contig=None``.
    """
    if not file_exists_and_valid(path):
        return []

    hits: List[BuscoHit] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            bid, status = cols[0], cols[1]
            if status == "Missing" or len(cols) < 8:
                hits.append(BuscoHit(bid, status, None, None, None, None, None, None))
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[5]
                score = float(cols[6])
                length = int(cols[7])
            except (ValueError, IndexError):
                hits.append(BuscoHit(bid, status, cols[2], None, None, None, None, None))
                continue
            hits.append(BuscoHit(bid, status, cols[2], start, end, strand, score, length))
    return hits


def _leaf_for_contig(
    contig_to_leaf: Dict[str, LeafId],
    contig: Optional[str],
) -> Optional[LeafId]:
    if contig is None:
        return None
    return contig_to_leaf.get(contig)


def build_assembly_leaves(
    assembly_name: str,
    classifications: Iterable[ContigClassification],
    compleasm_result: CompleasmResult,
) -> List[LeafBuscos]:
    """Split a query assembly's BUSCO hits into per-subgenome leaves.

    Hits on contigs whose ``classification != 'chrom_assigned'`` are dropped
    (we only consider chromosome-anchored BUSCOs).  Each contig contributes
    to a leaf keyed by the *reference* subgenome encoded in
    ``assigned_ref_id`` (chr3A → ``A``) and the *query* subgenome assigned
    by polyploid duplication inference.  For a polyploid query against a
    polyploid reference this produces up to ``n_ref_sg × n_query_sg`` leaves
    per assembly.
    """
    if compleasm_result.full_table_path is None:
        logger.warning(
            f"No full_table for assembly {assembly_name!r}; skipping in phylogeny"
        )
        return []

    contig_to_leaf: Dict[str, LeafId] = {}
    for clf in classifications:
        if clf.classification != "chrom_assigned":
            continue
        ref_sg: Optional[str] = None
        if clf.assigned_ref_id:
            _, ref_sg = split_chrom_subgenome(clf.assigned_ref_id)
        leaf = LeafId(
            assembly=assembly_name,
            ref_subgenome=ref_sg,
            query_subgenome=clf.query_subgenome,
        )
        # Index by both the post-rename and the original name so we can join
        # against compleasm hits regardless of which contig name was passed
        # to compleasm.
        contig_to_leaf[clf.new_name] = leaf
        contig_to_leaf[clf.original_name] = leaf

    hits = parse_full_table(compleasm_result.full_table_path)
    leaves_by_id: Dict[LeafId, LeafBuscos] = {}
    for h in hits:
        if h.status not in _INCLUDE_STATUSES:
            continue
        leaf = _leaf_for_contig(contig_to_leaf, h.contig)
        if leaf is None:
            continue
        lb = leaves_by_id.setdefault(
            leaf,
            LeafBuscos(
                leaf=leaf,
                translated_protein_path=compleasm_result.translated_protein_path,
                n_lineage_total=compleasm_result.n_total,
            ),
        )
        lb.hits_by_busco.setdefault(h.busco_id, []).append(h)

    return list(leaves_by_id.values())


def build_reference_leaf(
    reference_name: str,
    compleasm_result: CompleasmResult,
) -> Optional[LeafBuscos]:
    """Build a single leaf representing the reference genome."""
    if compleasm_result.full_table_path is None:
        logger.warning("Reference compleasm has no full_table; cannot include reference in tree")
        return None
    leaf = LeafId(assembly=reference_name)
    lb = LeafBuscos(
        leaf=leaf,
        translated_protein_path=compleasm_result.translated_protein_path,
        n_lineage_total=compleasm_result.n_total,
    )
    for h in parse_full_table(compleasm_result.full_table_path):
        if h.status not in _INCLUDE_STATUSES:
            continue
        lb.hits_by_busco.setdefault(h.busco_id, []).append(h)
    return lb


def filter_leaves_by_completeness(
    leaves: List[LeafBuscos],
    min_pct: float,
) -> Tuple[List[LeafBuscos], List[Tuple[LeafBuscos, float]]]:
    """Drop leaves below ``min_pct`` single-copy BUSCO completeness.

    Returns ``(kept, dropped_with_pct)``.
    """
    kept: List[LeafBuscos] = []
    dropped: List[Tuple[LeafBuscos, float]] = []
    for lb in leaves:
        pct = lb.completeness_pct()
        if pct >= min_pct:
            kept.append(lb)
        else:
            dropped.append((lb, pct))
    return kept, dropped


def shared_single_copy_genes(leaves: List[LeafBuscos]) -> List[str]:
    """BUSCOs that are single-copy Complete in every leaf (sorted)."""
    if not leaves:
        return []
    common = set(leaves[0].single_copy_genes())
    for lb in leaves[1:]:
        common &= lb.single_copy_genes()
        if not common:
            return []
    return sorted(common)
