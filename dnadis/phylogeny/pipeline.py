"""Cross-assembly species-tree pipeline (top-level orchestrator).

Runs after all per-assembly pipelines have completed.  Submits any
heavy compute (reference compleasm, IQ-TREE) through the executor passed
in so SLURM submission still works when ``--cluster`` is set.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from dnadis.detection.compleasm import run_compleasm
from dnadis.models import (
    AssemblyResult,
    CompleasmResult,
    PhylogenyResult,
    ReferenceContext,
)
from dnadis.phylogeny.busco_extraction import (
    LeafBuscos,
    build_assembly_leaves,
    build_reference_leaf,
    filter_leaves_by_completeness,
    shared_single_copy_genes,
)
from dnadis.phylogeny.supermatrix import (
    align_and_trim_genes,
    concatenate_supermatrix,
    write_per_gene_fastas,
)
from dnadis.phylogeny.tree import build_tree, resolve_outgroup
from dnadis.utils.logging_config import get_logger

logger = get_logger("phylogeny")

# A tree needs at least four taxa to be biologically interesting and at least
# a handful of genes to give bootstrap values any meaning.
_MIN_LEAVES = 3
_MIN_GENES = 4


def run_phylogeny(
    args,
    ref_ctx: ReferenceContext,
    results: List[AssemblyResult],
    output_dir: Path,
    executor,
    cluster_config,
    reference_name: str,
) -> Optional[PhylogenyResult]:
    """Build a species tree from per-leaf single-copy BUSCOs.

    Returns ``None`` when the tree cannot be built (insufficient compleasm
    data, missing tools, too few taxa or shared genes).
    """
    logger.phase("=== Cross-assembly phylogenetics ===")

    phylo_dir = output_dir / "phylogeny"
    phylo_dir.mkdir(parents=True, exist_ok=True)

    # Submit reference compleasm early so it overlaps with per-gene MSA work.
    ref_compleasm_future = None
    if not args.phylo_skip_reference and args.compleasm_lineage:
        ref_compleasm_dir = output_dir / "reference" / "compleasm"
        use_cluster = cluster_config.enabled if cluster_config else False
        ref_compleasm_spec = None
        ref_threads = args.threads
        if use_cluster:
            from dnadis.utils.resource_estimation import estimate_compleasm_resources
            ref_compleasm_spec = estimate_compleasm_resources(
                ref_ctx.ref, cluster_config,
            )
            ref_threads = ref_compleasm_spec.cores
        ref_compleasm_future = executor.submit(
            run_compleasm,
            fasta=ref_ctx.ref,
            output_dir=ref_compleasm_dir,
            lineage=args.compleasm_lineage,
            threads=ref_threads,
            library_path=args.compleasm_library,
            compleasm_exe=args.compleasm_path,
            resource_spec=ref_compleasm_spec if use_cluster else None,
        )
    elif not args.compleasm_lineage:
        logger.info("Phylogeny skipped: no --compleasm-lineage configured")
        return None

    leaves: List[LeafBuscos] = []
    for r in results:
        if r.compleasm_chrs is None or r.compleasm_chrs.full_table_path is None:
            logger.info(
                f"Assembly {r.assembly_name!r} has no compleasm full_table — "
                f"excluding from phylogeny"
            )
            continue
        leaves.extend(build_assembly_leaves(r.assembly_name, r.classifications, r.compleasm_chrs))

    ref_leaf: Optional[LeafBuscos] = None
    ref_label: Optional[str] = None
    if ref_compleasm_future is not None:
        try:
            ref_compleasm: Optional[CompleasmResult] = ref_compleasm_future.result()
        except Exception as e:
            logger.warning(f"Reference compleasm failed: {e}")
            ref_compleasm = None
        if ref_compleasm is not None:
            ref_leaf = build_reference_leaf(reference_name, ref_compleasm)
            if ref_leaf is not None:
                leaves.append(ref_leaf)
                ref_label = ref_leaf.leaf.label

    if not leaves:
        logger.warning("Phylogeny: no usable leaves; skipping")
        return None

    kept, dropped = filter_leaves_by_completeness(leaves, args.phylo_min_busco_completeness)
    for lb, pct in dropped:
        logger.warning(
            f"Phylogeny: dropping leaf {lb.leaf.label!r} (single-copy BUSCO {pct:.1f}% "
            f"below threshold {args.phylo_min_busco_completeness}%)"
        )

    if len(kept) < _MIN_LEAVES:
        logger.warning(
            f"Phylogeny: only {len(kept)} leaves pass the completeness filter "
            f"(need ≥ {_MIN_LEAVES}); skipping tree"
        )
        return None

    shared = shared_single_copy_genes(kept)
    if len(shared) < _MIN_GENES:
        logger.warning(
            f"Phylogeny: only {len(shared)} BUSCOs are single-copy across all "
            f"{len(kept)} leaves (need ≥ {_MIN_GENES}); skipping tree"
        )
        return None
    logger.info(
        f"Phylogeny: {len(kept)} leaves, {len(shared)} shared single-copy BUSCOs"
    )

    per_gene_dir = phylo_dir / "per_gene"
    gene_fastas, leaf_labels = write_per_gene_fastas(kept, shared, per_gene_dir)
    if not gene_fastas:
        logger.warning("Phylogeny: no per-gene FASTAs written; skipping")
        return None

    msa_dir = phylo_dir / "msa"
    use_cluster = cluster_config.enabled if cluster_config else False
    trimmed = align_and_trim_genes(
        gene_fastas, msa_dir,
        threads=args.threads,
        executor=executor,
        use_cluster=use_cluster,
        cluster_config=cluster_config,
        chunk_size=args.phylo_msa_chunk_size,
        inner_threads=args.phylo_msa_inner_threads,
    )
    if not trimmed:
        logger.warning("Phylogeny: no per-gene MSAs produced; skipping")
        return None

    supermatrix_fasta = phylo_dir / "supermatrix.faa"
    sm = concatenate_supermatrix(trimmed, leaf_labels, supermatrix_fasta)
    if sm.n_genes_used == 0 or sm.n_columns == 0:
        logger.warning("Phylogeny: empty supermatrix; skipping")
        return None
    logger.info(
        f"Phylogeny supermatrix: {sm.n_genes_used} genes, {sm.n_columns} columns, "
        f"{len(sm.leaf_order)} taxa"
    )

    outgroup_label = resolve_outgroup(
        args.phylo_outgroup,
        sm.leaf_order,
        ref_label,
        supermatrix_fasta,
    )

    iqtree_prefix = phylo_dir / "species_tree"
    iqtree_err = phylo_dir / "iqtree.err"

    iqtree_threads = args.phylo_iqtree_threads
    if use_cluster and cluster_config is not None:
        from dnadis.utils.distributed import ResourceSpec, clamp_resources
        iqtree_spec = clamp_resources(
            ResourceSpec(
                cores=iqtree_threads,
                memory_gb=float(args.phylo_max_mem_gb),
                time_minutes=args.phylo_iqtree_time_minutes,
                job_name="dnadis_iqtree",
            ),
            cluster_config,
        )
        iqtree_threads = iqtree_spec.cores  # may have been clamped down
        logger.info(
            f"Submitting IQ-TREE to SLURM: cores={iqtree_spec.cores}, "
            f"mem={iqtree_spec.memory_gb:g}G, time={iqtree_spec.time_minutes}min"
        )
        treefile = executor.submit(
            build_tree,
            supermatrix_fasta=supermatrix_fasta,
            out_prefix=iqtree_prefix,
            threads=iqtree_threads,
            max_mem_gb=args.phylo_max_mem_gb,
            bootstrap=args.phylo_bootstrap_reps,
            alrt=args.phylo_alrt_reps,
            models=args.phylo_models,
            outgroup_label=outgroup_label,
            err_path=iqtree_err,
            resource_spec=iqtree_spec,
        ).result()
    else:
        treefile = build_tree(
            supermatrix_fasta=supermatrix_fasta,
            out_prefix=iqtree_prefix,
            threads=iqtree_threads,
            max_mem_gb=args.phylo_max_mem_gb,
            bootstrap=args.phylo_bootstrap_reps,
            alrt=args.phylo_alrt_reps,
            models=args.phylo_models,
            outgroup_label=outgroup_label,
            err_path=iqtree_err,
        )
    if treefile is None:
        logger.warning("Phylogeny: IQ-TREE failed; no tree produced")
        return None

    logger.done(f"Species tree: {treefile}")
    return PhylogenyResult(
        supermatrix_fasta=supermatrix_fasta,
        treefile=treefile,
        leaf_order=sm.leaf_order,
        n_genes_used=sm.n_genes_used,
        n_columns=sm.n_columns,
        outgroup=outgroup_label,
        iqtree_prefix=iqtree_prefix,
        dropped_leaves=[(lb.leaf.label, pct) for lb, pct in dropped],
    )
