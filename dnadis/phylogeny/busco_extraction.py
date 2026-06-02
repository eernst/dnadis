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

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from dnadis.models import CompleasmResult, ContigClassification
from dnadis.utils.io_utils import file_exists_and_valid
from dnadis.utils.logging_config import get_logger
from dnadis.utils.reference_utils import is_nuclear_chromosome, split_chrom_subgenome

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
        # ``split_chrom_subgenome`` returns the literal string ``"NA"`` for
        # monoploid references; treat it the same as None so the leaf label
        # for a monoploid reference stays bare ``reference_name``.
        if self.ref_subgenome and self.ref_subgenome != "NA":
            parts.append(self.ref_subgenome)
        if self.query_subgenome and self.query_subgenome != "NA":
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
            if ref_sg == "NA":
                ref_sg = None
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


def build_reference_leaves(
    reference_name: str,
    compleasm_result: CompleasmResult,
) -> List[LeafBuscos]:
    """Build one leaf per reference subgenome.

    Reference contig names directly encode subgenome identity (chr3A → "A",
    chr3P → "P"); we split the reference's full_table by
    :func:`split_chrom_subgenome` of each row's contig.

    * Monoploid reference: every row has ``ref_sg = None`` → one leaf
      labeled with just ``reference_name``.
    * Polyploid composed reference (e.g. three haploidized assemblies
      stitched together with chr*A / chr*P / chr*T suffixes): one leaf per
      subgenome (``reference_name_A``, ``reference_name_P``, etc.).  This
      mirrors how query assemblies are split — without it a polyploid
      reference looks like every BUSCO is c-fold duplicated and the
      ``--phylo-min-busco-completeness`` gate drops the whole reference.
    """
    if compleasm_result.full_table_path is None:
        logger.warning("Reference compleasm has no full_table; cannot include reference in tree")
        return []

    leaves_by_id: Dict[LeafId, LeafBuscos] = {}
    for h in parse_full_table(compleasm_result.full_table_path):
        if h.status not in _INCLUDE_STATUSES:
            continue
        if h.contig is None:
            continue
        _, ref_sg = split_chrom_subgenome(h.contig)
        # split_chrom_subgenome returns "NA" (string) for monoploid refs;
        # normalise that to None so the leaf label stays bare "reference_name"
        # rather than "reference_name_NA".
        if ref_sg == "NA":
            ref_sg = None
        leaf_id = LeafId(assembly=reference_name, ref_subgenome=ref_sg)
        lb = leaves_by_id.setdefault(
            leaf_id,
            LeafBuscos(
                leaf=leaf_id,
                translated_protein_path=compleasm_result.translated_protein_path,
                n_lineage_total=compleasm_result.n_total,
            ),
        )
        lb.hits_by_busco.setdefault(h.busco_id, []).append(h)
    return list(leaves_by_id.values())


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


# ---------------------------------------------------------------------------
# Outgroup helpers
#
# When the user designates a specific query assembly as the outgroup, we
# relax the chrom_assigned filter so distant outgroups whose contigs largely
# don't align to the reference can still anchor the tree.  This carries two
# risks the helpers below try to detect:
#
#   * polyploid outgroup without subgenome separation: c >= 2 hits per
#     BUSCO, but very few outgroup contigs map to reference chromosomes —
#     subgenome identity is unrecoverable, so the outgroup can't contribute
#     meaningful per-subgenome leaves.
#   * partially-mappable polyploid outgroup: c >= 2 and one subgenome has
#     dense reference coverage while another doesn't.  We keep the well-
#     mapped subgenome(s) and drop the rest (per user judgement that the
#     non-mapping subgenome is too hard to distinguish from debris).
# ---------------------------------------------------------------------------


def infer_assembly_ploidy(hits: List[BuscoHit]) -> int:
    """Estimate ploidy ``c`` from per-BUSCO hit-count distribution.

    Returns the median number of present-status hits per BUSCO across all
    BUSCOs that compleasm detected.  ``c == 1`` indicates a haploid /
    collapsed assembly; ``c >= 2`` indicates polyploidy.  Returns ``1`` when
    there are no present-status hits (the caller has bigger problems than
    ploidy estimation).
    """
    per_busco: Dict[str, int] = defaultdict(int)
    for h in hits:
        if h.status in _INCLUDE_STATUSES:
            per_busco[h.busco_id] += 1
    if not per_busco:
        return 1
    counts = sorted(per_busco.values())
    return counts[len(counts) // 2]


def ref_assignment_quality_by_subgenome(
    classifications: Iterable[ContigClassification],
    ref_lengths_norm: Dict[str, int],
) -> Dict[Optional[str], float]:
    """Per-reference-subgenome assignment coverage for one query assembly.

    For each ref subgenome (chr1A → ``A``, chr1P → ``P``, etc.) returns the
    fraction of that subgenome's reference chromosomes that have at least
    one query contig assigned to them.  A distant outgroup that picks up
    7 / 63 chromosome assignments across three subgenomes scores low on
    every subgenome (~0.03–0.10); a close outgroup that maps cleanly to
    one subgenome scores ~1.0 on that subgenome and ~0 on the others.

    Organelles and non-nuclear references are excluded from both numerator
    and denominator.
    """
    sg_to_assigned: Dict[Optional[str], set] = defaultdict(set)
    for clf in classifications:
        if clf.classification != "chrom_assigned" or not clf.assigned_ref_id:
            continue
        _, sg = split_chrom_subgenome(clf.assigned_ref_id)
        sg_to_assigned[sg].add(clf.assigned_ref_id)

    sg_to_total: Dict[Optional[str], int] = defaultdict(int)
    for ref_id in ref_lengths_norm:
        if not is_nuclear_chromosome(ref_id):
            continue
        _, sg = split_chrom_subgenome(ref_id)
        sg_to_total[sg] += 1

    out: Dict[Optional[str], float] = {}
    for sg, assigned in sg_to_assigned.items():
        total = sg_to_total.get(sg, 0)
        out[sg] = (len(assigned) / total) if total > 0 else 0.0
    return out


def build_outgroup_only_leaves(
    assembly_name: str,
    compleasm_result: CompleasmResult,
) -> List[LeafBuscos]:
    """Build leaves for a phylogeny-only outgroup (``--phylo-outgroup-only``).

    Such an assembly runs a minimal pipeline (compleasm only) and has no
    classifications, so subgenome identity cannot come from reference
    assignment.  Instead we read it straight from the outgroup's *own* contig
    names via :func:`split_chrom_subgenome`: contigs sharing an alpha suffix
    (chr1A / chr2A → ``A``, chr1B → ``B``) group into one leaf per suffix,
    placed in :attr:`LeafId.query_subgenome`.  Contigs without a recognizable
    subgenome suffix pool into a single leaf labeled bare ``assembly_name`` —
    the common case for a distant outgroup that doesn't align to the reference.

    Unlike :func:`build_outgroup_leaves`, this never inspects reference
    assignment or infers ploidy, so it works for an outgroup whose contigs
    largely don't map to the reference at all.
    """
    if compleasm_result.full_table_path is None:
        logger.warning(
            f"Outgroup {assembly_name!r} has no compleasm full_table; "
            f"excluding from phylogeny"
        )
        return []

    leaves_by_id: Dict[LeafId, LeafBuscos] = {}
    for h in parse_full_table(compleasm_result.full_table_path):
        if h.status not in _INCLUDE_STATUSES:
            continue
        if h.contig is None:
            continue
        _, query_sg = split_chrom_subgenome(h.contig)
        # "NA" (string) means no subgenome suffix; normalise to None so the
        # leaf label stays bare "assembly_name" rather than "assembly_name_NA".
        if query_sg == "NA":
            query_sg = None
        leaf_id = LeafId(assembly=assembly_name, query_subgenome=query_sg)
        lb = leaves_by_id.setdefault(
            leaf_id,
            LeafBuscos(
                leaf=leaf_id,
                translated_protein_path=compleasm_result.translated_protein_path,
                n_lineage_total=compleasm_result.n_total,
            ),
        )
        lb.hits_by_busco.setdefault(h.busco_id, []).append(h)
    return list(leaves_by_id.values())


def build_outgroup_leaves(
    assembly_name: str,
    classifications: Iterable[ContigClassification],
    compleasm_result: CompleasmResult,
    ref_lengths_norm: Dict[str, int],
    min_ref_assignment: float = 0.5,
) -> Tuple[List[LeafBuscos], str]:
    """Build leaves for a designated outgroup with relaxed chrom-assigned rules.

    Returns ``(leaves, status)`` where ``status`` is one of:

    * ``"haploid_pooled"`` — c = 1; all single-copy BUSCOs in the outgroup
      contribute to a single leaf regardless of whether their contig was
      chromosome-assigned.
    * ``"polyploid_filtered"`` — c >= 2 and at least one reference subgenome
      has assignment quality >= ``min_ref_assignment``; only those well-
      assigned subgenomes contribute leaves.
    * ``"polyploid_unusable"`` — c >= 2 but no subgenome is well-assigned;
      caller should warn and skip the outgroup.
    """
    if compleasm_result.full_table_path is None:
        return [], "polyploid_unusable"

    hits = parse_full_table(compleasm_result.full_table_path)
    ploidy = infer_assembly_ploidy(hits)

    if ploidy <= 1:
        leaf = LeafId(assembly=assembly_name)
        lb = LeafBuscos(
            leaf=leaf,
            translated_protein_path=compleasm_result.translated_protein_path,
            n_lineage_total=compleasm_result.n_total,
        )
        for h in hits:
            if h.status not in _INCLUDE_STATUSES:
                continue
            lb.hits_by_busco.setdefault(h.busco_id, []).append(h)
        return [lb], "haploid_pooled"

    quality = ref_assignment_quality_by_subgenome(classifications, ref_lengths_norm)
    usable_sgs = {
        sg for sg, q in quality.items()
        if sg is not None and q >= min_ref_assignment
    }
    if not usable_sgs:
        return [], "polyploid_unusable"

    contig_to_leaf: Dict[str, LeafId] = {}
    for clf in classifications:
        if clf.classification != "chrom_assigned" or not clf.assigned_ref_id:
            continue
        _, ref_sg = split_chrom_subgenome(clf.assigned_ref_id)
        if ref_sg == "NA":
            ref_sg = None
        if ref_sg not in usable_sgs:
            continue
        leaf = LeafId(
            assembly=assembly_name,
            ref_subgenome=ref_sg,
            query_subgenome=clf.query_subgenome,
        )
        contig_to_leaf[clf.new_name] = leaf
        contig_to_leaf[clf.original_name] = leaf

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

    return list(leaves_by_id.values()), "polyploid_filtered"
