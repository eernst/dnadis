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
    build_outgroup_leaves,
    build_reference_leaves,
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


def submit_reference_compleasm(
    args,
    ref_ctx: ReferenceContext,
    output_dir: Path,
    executor,
    cluster_config,
):
    """Submit the reference compleasm job at the start of the executor scope.

    Called from ``main()`` right after the executor opens, in parallel with
    the per-assembly pipelines, so the reference run doesn't become a
    straggler that only starts once every assembly has finished.  Returns
    a future whose ``.result()`` yields a :class:`CompleasmResult` (or
    ``None``), or ``None`` if the reference run is not needed.
    """
    if args.skip_phylogeny:
        return None
    if not args.compleasm_lineage or args.skip_compleasm:
        return None
    if args.phylo_skip_reference:
        return None

    ref_compleasm_dir = output_dir / "reference" / "compleasm"
    use_cluster = cluster_config.enabled if cluster_config else False
    ref_compleasm_spec = None
    ref_threads = args.threads
    if use_cluster:
        from dnadis.utils.resource_estimation import estimate_compleasm_resources
        ref_compleasm_spec = estimate_compleasm_resources(ref_ctx.ref, cluster_config)
        ref_threads = ref_compleasm_spec.cores
    return executor.submit(
        run_compleasm,
        fasta=ref_ctx.ref,
        output_dir=ref_compleasm_dir,
        lineage=args.compleasm_lineage,
        threads=ref_threads,
        library_path=args.compleasm_library,
        compleasm_exe=args.compleasm_path,
        resource_spec=ref_compleasm_spec if use_cluster else None,
    )


def run_phylogeny(
    args,
    ref_ctx: ReferenceContext,
    results: List[AssemblyResult],
    output_dir: Path,
    executor,
    cluster_config,
    reference_name: str,
    ref_compleasm_future=None,
) -> Optional[PhylogenyResult]:
    """Build a species tree from per-leaf single-copy BUSCOs.

    ``ref_compleasm_future`` is the future returned by
    :func:`submit_reference_compleasm` at the top of the executor scope;
    pass ``None`` here only when running without a reference leaf.

    Returns ``None`` when the tree cannot be built (insufficient compleasm
    data, missing tools, too few taxa or shared genes).
    """
    logger.phase("=== Cross-assembly phylogenetics ===")

    phylo_dir = output_dir / "phylogeny"
    phylo_dir.mkdir(parents=True, exist_ok=True)

    if not args.compleasm_lineage:
        logger.info("Phylogeny skipped: no --compleasm-lineage configured")
        return None

    outgroup_arg = (args.phylo_outgroup or "none").strip()
    # Special outgroup handling only applies when the user names a specific
    # query assembly.  "none"/"reference"/"auto" follow the standard path.
    outgroup_asm_name: Optional[str] = None
    if outgroup_arg.lower() not in {"none", "reference", "auto", ""}:
        for r in results:
            if r.assembly_name == outgroup_arg:
                outgroup_asm_name = r.assembly_name
                break
        if outgroup_asm_name is None:
            logger.warning(
                f"--phylo-outgroup={outgroup_arg!r} does not name any query assembly; "
                f"outgroup will be resolved against built leaves after the tree is built."
            )

    leaves: List[LeafBuscos] = []
    for r in results:
        compleasm_for_phylo = r.compleasm_full or r.compleasm_chrs
        if compleasm_for_phylo is None or compleasm_for_phylo.full_table_path is None:
            logger.info(
                f"Assembly {r.assembly_name!r} has no compleasm full_table — "
                f"excluding from phylogeny"
            )
            continue
        if r.assembly_name == outgroup_asm_name:
            og_leaves, og_status = build_outgroup_leaves(
                assembly_name=r.assembly_name,
                classifications=r.classifications,
                compleasm_result=compleasm_for_phylo,
                ref_lengths_norm=ref_ctx.ref_lengths_norm,
                min_ref_assignment=args.phylo_outgroup_min_ref_assignment,
            )
            if og_status == "haploid_pooled":
                logger.info(
                    f"Outgroup {r.assembly_name!r}: detected c=1; pooled all "
                    f"single-copy BUSCOs into one leaf regardless of chromosome assignment"
                )
            elif og_status == "polyploid_filtered":
                kept_sgs = sorted({lb.leaf.ref_subgenome for lb in og_leaves if lb.leaf.ref_subgenome})
                logger.info(
                    f"Outgroup {r.assembly_name!r}: c>=2; retained subgenome(s) "
                    f"{kept_sgs} with ref-chromosome-assignment quality >= "
                    f"{args.phylo_outgroup_min_ref_assignment:.0%}; other subgenomes dropped"
                )
            elif og_status == "polyploid_unusable":
                logger.warning(
                    f"Outgroup {r.assembly_name!r}: c>=2 but no reference subgenome reached "
                    f"the assignment-quality threshold ({args.phylo_outgroup_min_ref_assignment:.0%}). "
                    f"This usually means the outgroup is too divergent for subgenome "
                    f"separation against this reference.  Excluding it from the tree; "
                    f"phylogeny will proceed unrooted unless --phylo-outgroup matches another taxon."
                )
            leaves.extend(og_leaves)
        else:
            leaves.extend(build_assembly_leaves(r.assembly_name, r.classifications, compleasm_for_phylo))

    ref_labels: List[str] = []
    if ref_compleasm_future is not None:
        try:
            ref_compleasm: Optional[CompleasmResult] = ref_compleasm_future.result()
        except Exception as e:
            logger.warning(f"Reference compleasm failed: {e}")
            ref_compleasm = None
        if ref_compleasm is not None:
            ref_leaves = build_reference_leaves(reference_name, ref_compleasm)
            if len(ref_leaves) > 1:
                logger.info(
                    f"Reference is polyploid: split into {len(ref_leaves)} per-subgenome "
                    f"leaves ({[lb.leaf.label for lb in ref_leaves]})"
                )
            leaves.extend(ref_leaves)
            ref_labels = [lb.leaf.label for lb in ref_leaves]

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
        ref_labels,
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

    # Sidecar so refresh_reports.py (and any downstream consumer) can recover
    # the outgroup label without re-deriving it from the run.
    outgroup_sidecar = phylo_dir / "species_tree.outgroup.txt"
    outgroup_sidecar.write_text(f"{outgroup_label or ''}\n")

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
