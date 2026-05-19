"""Build a concatenated protein supermatrix from per-leaf single-copy BUSCOs.

For each BUSCO that is single-copy in every retained leaf, extract the leaf's
protein sequence, align with MAFFT L-INS-i, trim columns with trimAl
-gappyout, then horizontally concatenate per-gene alignments into one
supermatrix with consistent taxon order.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from dnadis.alignment.external_tools import run_mafft, run_trimal
from dnadis.phylogeny.busco_extraction import BuscoHit, LeafBuscos, LeafId
from dnadis.utils.logging_config import get_logger

logger = get_logger("phylogeny")


@dataclass
class SupermatrixResult:
    supermatrix_fasta: Path
    n_genes_used: int
    n_columns: int
    leaf_order: List[str]
    per_gene_dir: Path


def _parse_translated_protein_fasta(path: Path) -> Dict[Tuple[str, str, int, int], str]:
    """Index compleasm's translated_protein.fasta by (busco_id, contig, start, end).

    Header format written by compleasm::

        >{Protein_name}|{Contig_name}:{Contig_Start}-{Contig_Stop}

    where ``Protein_name`` starts with ``{busco_id}_``.
    """
    out: Dict[Tuple[str, str, int, int], str] = {}
    if path is None or not path.exists():
        return out

    cur_key: Optional[Tuple[str, str, int, int]] = None
    cur_seq: List[str] = []

    def flush() -> None:
        if cur_key is not None:
            out[cur_key] = "".join(cur_seq)

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                cur_seq = []
                header = line[1:]
                try:
                    name_part, coord_part = header.split("|", 1)
                    contig, span = coord_part.rsplit(":", 1)
                    start_s, end_s = span.split("-", 1)
                    busco_id = name_part.split("_", 1)[0]
                    cur_key = (busco_id, contig, int(start_s), int(end_s))
                except ValueError:
                    cur_key = None
            else:
                cur_seq.append(line.strip())
        flush()

    return out


def _select_single_copy_hit(lb: LeafBuscos, busco_id: str) -> Optional[BuscoHit]:
    """Return the lone hit for ``busco_id`` in ``lb`` (regardless of compleasm's
    global Complete/Duplicated label), or ``None`` if there is ≠1 hit."""
    hits = lb.hits_by_busco.get(busco_id, [])
    if len(hits) != 1:
        return None
    h = hits[0]
    if h.contig is None or h.start is None or h.end is None:
        return None
    return h


def _resolve_sequence(
    protein_index: Dict[Tuple[str, str, int, int], str],
    busco_id: str,
    hit: BuscoHit,
) -> Optional[str]:
    """Look up a hit's protein, tolerating off-by-one differences in coords.

    Compleasm's full_table coordinates are usually identical to the headers
    in ``translated_protein.fasta``, but we try a few neighbours to be safe.
    """
    contig = hit.contig
    assert contig is not None and hit.start is not None and hit.end is not None
    key = (busco_id, contig, hit.start, hit.end)
    if key in protein_index:
        return protein_index[key]
    for ds in (-1, 1, -2, 2):
        for de in (0, -1, 1):
            k = (busco_id, contig, hit.start + ds, hit.end + de)
            if k in protein_index:
                return protein_index[k]
    return None


def write_per_gene_fastas(
    leaves: List[LeafBuscos],
    shared_genes: List[str],
    out_dir: Path,
) -> Tuple[Dict[str, Path], List[str]]:
    """Write one per-gene unaligned FASTA per shared single-copy BUSCO.

    Returns ``(gene_to_fasta, leaf_labels)`` where the per-gene FASTA lists
    one record per leaf in the same order as ``leaf_labels``.  Genes for
    which any leaf's protein sequence cannot be resolved are skipped.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    leaf_labels = [lb.leaf.label for lb in leaves]

    # Per-leaf protein index (parse each translated_protein.fasta once).
    by_path: Dict[Path, Dict[Tuple[str, str, int, int], str]] = {}
    for lb in leaves:
        if lb.translated_protein_path is None:
            continue
        if lb.translated_protein_path not in by_path:
            by_path[lb.translated_protein_path] = _parse_translated_protein_fasta(
                lb.translated_protein_path
            )

    gene_to_fasta: Dict[str, Path] = {}
    for bid in shared_genes:
        records: List[Tuple[str, str]] = []
        ok = True
        for lb in leaves:
            hit = _select_single_copy_hit(lb, bid)
            idx = by_path.get(lb.translated_protein_path, {}) if lb.translated_protein_path else {}
            if hit is None:
                ok = False
                break
            seq = _resolve_sequence(idx, bid, hit)
            if seq is None or not seq:
                ok = False
                break
            records.append((lb.leaf.label, seq))
        if not ok:
            logger.debug(f"Skipping {bid}: missing protein for at least one leaf")
            continue
        fasta = out_dir / f"{bid}.faa"
        with fasta.open("w") as out:
            for name, seq in records:
                out.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
        gene_to_fasta[bid] = fasta

    return gene_to_fasta, leaf_labels


def _align_and_trim_one(
    bid: str,
    fasta: Path,
    aln_path: Path,
    trim_path: Path,
    mafft_err: Path,
    trimal_err: Path,
) -> Optional[Tuple[str, Path]]:
    """Run MAFFT L-INS-i then trimAl -gappyout on one gene.

    Module-level (not nested) so it can be dispatched to a remote worker.
    """
    if not run_mafft(
        fasta, aln_path,
        threads=1,
        err_path=mafft_err,
        algorithm="linsi",
        quiet=True,
    ):
        return None
    if not run_trimal(
        aln_path, trim_path,
        mode="gappyout",
        err_path=trimal_err,
        quiet=True,
    ):
        return None
    return bid, trim_path


def _align_and_trim_chunk(
    chunk: List[Tuple[str, Path, Path, Path, Path, Path]],
    inner_threads: int = 4,
) -> List[Tuple[str, Path]]:
    """Process a chunk of genes inside a single SLURM job.

    ``chunk`` is a list of ``(bid, fasta, aln_path, trim_path, mafft_err, trimal_err)``
    tuples.  Genes within the chunk run concurrently on a per-job
    ``ThreadPoolExecutor`` with ``inner_threads`` workers.

    Returns the list of successfully completed ``(bid, trim_path)`` pairs;
    failures are silently dropped (per-gene errors are written to the
    per-gene .err files already).
    """
    out: List[Tuple[str, Path]] = []
    if inner_threads <= 1 or len(chunk) <= 1:
        for args in chunk:
            r = _align_and_trim_one(*args)
            if r is not None:
                out.append(r)
        return out

    with ThreadPoolExecutor(max_workers=inner_threads) as pool:
        futures = [pool.submit(_align_and_trim_one, *args) for args in chunk]
        for fut in as_completed(futures):
            r = fut.result()
            if r is not None:
                out.append(r)
    return out


def align_and_trim_genes(
    gene_to_fasta: Dict[str, Path],
    msa_dir: Path,
    threads: int = 1,
    executor=None,
    use_cluster: bool = False,
    cluster_config=None,
    chunk_size: int = 256,
    inner_threads: int = 4,
) -> Dict[str, Path]:
    """Run MAFFT L-INS-i then trimAl -gappyout on each per-gene FASTA.

    Returns a dict mapping busco_id → trimmed MSA path.  Genes that fail
    either step are omitted.

    Scheduling:
        * ``use_cluster=True`` (and ``executor`` set): genes are grouped into
          chunks of ``chunk_size`` and each chunk is submitted as one SLURM
          job via the executor.  Each chunk job uses ``inner_threads`` cores
          internally so ``inner_threads`` genes within the chunk run
          concurrently.  Total SLURM jobs = ``ceil(N_genes / chunk_size)``.
        * Otherwise: a local ``ThreadPoolExecutor`` with ``threads`` workers
          runs the jobs on the coordinator node.
    """
    msa_dir.mkdir(parents=True, exist_ok=True)
    aln_dir = msa_dir / "mafft_linsi"
    trim_dir = msa_dir / "trimal_gappyout"
    aln_dir.mkdir(parents=True, exist_ok=True)
    trim_dir.mkdir(parents=True, exist_ok=True)

    n_genes = len(gene_to_fasta)

    def _paths(bid: str) -> Tuple[Path, Path, Path, Path]:
        return (
            aln_dir / f"{bid}.aln.faa",
            trim_dir / f"{bid}.trim.faa",
            aln_dir / f"{bid}.mafft.err",
            trim_dir / f"{bid}.trimal.err",
        )

    results: Dict[str, Path] = {}

    if use_cluster and executor is not None:
        from dnadis.utils.distributed import ResourceSpec, clamp_resources
        items = list(gene_to_fasta.items())
        chunks: List[List[Tuple[str, Path, Path, Path, Path, Path]]] = []
        for i in range(0, len(items), chunk_size):
            sub = items[i:i + chunk_size]
            chunks.append([
                (bid, fasta, *_paths(bid)) for bid, fasta in sub
            ])
        spec = ResourceSpec(
            cores=inner_threads,
            memory_gb=8.0,
            time_minutes=240,
            job_name="dnadis_phylo_msa",
        )
        if cluster_config is not None:
            spec = clamp_resources(spec, cluster_config)
        logger.info(
            f"Aligning {n_genes} per-gene FASTAs via SLURM: "
            f"{len(chunks)} chunked job(s) × ≤{chunk_size} genes "
            f"× {inner_threads} inner thread(s)..."
        )
        futures = [
            executor.submit(
                _align_and_trim_chunk,
                chunk=chunk,
                inner_threads=inner_threads,
                resource_spec=spec,
            )
            for chunk in chunks
        ]
        for fut in futures:
            try:
                chunk_results = fut.result()
            except Exception as e:
                logger.warning(f"Per-gene MSA chunk failed: {e}")
                continue
            for bid, trim_path in chunk_results:
                results[bid] = trim_path
    elif threads <= 1:
        logger.info(
            f"Aligning {n_genes} per-gene FASTAs (MAFFT L-INS-i + trimAl -gappyout, "
            f"serial)..."
        )
        for bid, fasta in gene_to_fasta.items():
            aln, trimmed, mafft_err, trimal_err = _paths(bid)
            r = _align_and_trim_one(bid, fasta, aln, trimmed, mafft_err, trimal_err)
            if r is not None:
                results[r[0]] = r[1]
    else:
        logger.info(
            f"Aligning {n_genes} per-gene FASTAs (MAFFT L-INS-i + trimAl -gappyout, "
            f"{threads} worker thread(s))..."
        )
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = {}
            for bid, fasta in gene_to_fasta.items():
                aln, trimmed, mafft_err, trimal_err = _paths(bid)
                futures[pool.submit(
                    _align_and_trim_one, bid, fasta, aln, trimmed, mafft_err, trimal_err
                )] = bid
            for fut in as_completed(futures):
                r = fut.result()
                if r is not None:
                    results[r[0]] = r[1]

    logger.info(f"Per-gene MSAs: {len(results)}/{n_genes} succeeded")
    return results


def _read_aligned_fasta(path: Path) -> Dict[str, str]:
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


def concatenate_supermatrix(
    trimmed_msas: Dict[str, Path],
    leaf_order: List[str],
    out_fasta: Path,
) -> SupermatrixResult:
    """Horizontally concatenate per-gene trimmed MSAs into one supermatrix.

    Per-leaf rows are written in ``leaf_order``.  A leaf missing entirely
    from a gene MSA is padded with gap characters of the alignment width
    (this should not happen since per-gene FASTAs are built from shared
    single-copy genes, but the safety net keeps the matrix rectangular).
    """
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    concat: Dict[str, List[str]] = {label: [] for label in leaf_order}
    n_cols = 0
    n_genes = 0

    for bid in sorted(trimmed_msas.keys()):
        msa = _read_aligned_fasta(trimmed_msas[bid])
        if not msa:
            continue
        width = len(next(iter(msa.values())))
        if any(len(s) != width for s in msa.values()):
            logger.warning(f"Skipping {bid}: ragged MSA")
            continue
        for label in leaf_order:
            concat[label].append(msa.get(label, "-" * width))
        n_cols += width
        n_genes += 1

    with out_fasta.open("w") as out:
        for label in leaf_order:
            seq = "".join(concat[label])
            out.write(f">{label}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")

    return SupermatrixResult(
        supermatrix_fasta=out_fasta,
        n_genes_used=n_genes,
        n_columns=n_cols,
        leaf_order=list(leaf_order),
        per_gene_dir=out_fasta.parent,
    )
