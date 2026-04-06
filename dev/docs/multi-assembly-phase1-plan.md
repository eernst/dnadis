# Multi-Assembly Support — Phase 1: Core Infrastructure

## Context

The pipeline currently processes exactly one query assembly per invocation. The monolithic `main()` in `cli.py` (~1400 lines) handles reference preparation, per-assembly analysis, and output generation in a single function. To support analyzing multiple assemblies against a common reference, we need to: (1) refactor `main()` to separate shared reference work from per-assembly work, (2) add input mechanisms for specifying multiple assemblies, and (3) route per-assembly output into subfolders.

This plan covers **infrastructure only** — the comparative visualizations (hybrid subgenome configuration summary, chromosome synteny plot) will be designed and implemented separately after this foundation is in place.

---

## 1. New data structures in `models.py`

Add a `ReferenceContext` dataclass to hold all shared reference state:

```python
@dataclass
class ReferenceContext:
    ref: Path
    ref_lengths_norm: Dict[str, int]
    ref_orig_to_norm: Dict[str, str]
    ref_norm_to_orig: Dict[str, str]
    ref_ids_raw: Set[str]
    ref_gc_all: Dict[str, float]
    ref_gc_mean: Optional[float]
    ref_gc_std: Optional[float]
    ref_gff3: Optional[Path]          # filtered GFF3 (or None)
    tx2loc: Dict                       # transcript -> (chrom, start, end, strand)
    tx2gene: Dict                      # transcript -> gene_id
    proteins_faa: Optional[Path]       # extracted proteins (protein mode only)
    chrC_ref: Optional[Path]
    chrM_ref: Optional[Path]
    rdna_ref: Optional[Path]
    chr_like_minlen: int
    ref_lengths_tsv: Path
```

No `AssemblyResult` dataclass needed yet (that's for Phase 2 comparative analyses).

## 2. Refactor `main()` into three functions

### `prepare_reference(args, ref_outprefix) -> ReferenceContext`

Extract reference reading, ID maps, GC baseline, ref_lengths.tsv, GFF3 processing, protein extraction, organelle/rDNA reference preparation, and chr_like_minlen computation.

**Key detail**: Reference-only outputs (ref_lengths.tsv, proteins_faa, filtered GFF3) use `ref_outprefix`. In single-assembly mode, `ref_outprefix == outprefix` (backward-compatible). In multi-assembly mode, `ref_outprefix = output_dir / "reference" / "reference"`.

### `run_assembly(args, ref_ctx, qry, assembly_name, reads, outprefix)`

Extract the per-assembly body: read query FASTA, query GC, synteny analysis, gate filtering, detection phases, classification, output, and plotting.

Takes `ReferenceContext` instead of recomputing reference data. Receives `qry`, `assembly_name`, `reads`, and `outprefix` as parameters.

### `main()` becomes a thin dispatcher

```python
def main():
    args = parse_args()
    # ... logging, config ...

    if has_multi_assembly_input(args):
        assemblies = resolve_assemblies(args)
        output_dir = Path(args.output_dir)
        ref_outprefix = output_dir / "reference" / "reference"
        ref_ctx = prepare_reference(args, ref_outprefix)
        for asm_path, asm_name, asm_reads in assemblies:
            asm_outprefix = output_dir / asm_name / asm_name
            run_assembly(args, ref_ctx, asm_path, asm_name, asm_reads, asm_outprefix)
    else:
        outprefix = Path(args.outprefix)
        ref_ctx = prepare_reference(args, outprefix)
        run_assembly(args, ref_ctx, args.query, args.assembly_name, args.reads, outprefix)
```

## 3. CLI argument changes

Make `-q`/`--query` and `-o`/`--outprefix` no longer `required=True`. Instead, validate post-parse:

- **Single-assembly mode**: `-q` and `-o` both required (error if either is missing).
- **Multi-assembly mode**: `--fofn` or `--assembly-dir` required, plus `--output-dir`. `-q`/`-o` must NOT be set.

New arguments in a "Multi-assembly" group:

```python
multi = p.add_argument_group("Multi-assembly mode")
multi.add_argument("--fofn", type=_validate_input_path,
    help="File-of-filenames: TSV with columns [path] and optional [name], [reads]")
multi.add_argument("--assembly-dir", type=_validate_input_dir,
    help="Directory of FASTA files to process as multiple assemblies")
multi.add_argument("--output-dir", type=str,
    help="Output directory for multi-assembly mode (required with --fofn/--assembly-dir)")
```

## 4. Input resolution functions

**New file**: `dnadis/utils/multi_assembly.py`

- `parse_fofn()`: Read TSV with header, required column `path`, optional `name`/`reads`
- `scan_assembly_dir()`: Glob for FASTA files, derive names by stripping extensions
- `resolve_assemblies()`: Dispatch to parse_fofn or scan_assembly_dir

## 5. Output directory structure

Multi-assembly mode:
```
{output_dir}/
  reference/
    reference.ref_lengths.tsv
    reference.ref_proteins.fa
    reference.ref_gff3.filtered.gff3
  {assembly_name_1}/
    {assembly_name_1}.contig_summary.tsv
    ...
  {assembly_name_2}/
    ...
```

Single-assembly mode: Unchanged.

## 6. Error handling

- If one assembly fails, log the error and continue processing remaining assemblies
- At the end, summarize: "N/M assemblies completed successfully"
- Single-assembly mode: errors propagate as before

## 7. Logging

- Each assembly starts with a phase banner: `"=== Assembly: {name} ({i}/{n}) ==="`

---

## Implementation Status

**Completed** (2026-02-14):
- `models.py`: Added `ReferenceContext` dataclass
- `utils/multi_assembly.py`: Created with `parse_fofn()`, `scan_assembly_dir()`, `resolve_assemblies()`
- `cli.py`: Extracted `prepare_reference()` and `run_assembly()`, added multi-assembly CLI args, validation, dispatch loop
- `utils/config.py`: Added `multi_assembly` section to schema and template
- `CLAUDE.md`: Documented multi-assembly mode
- `tests/test_multi_assembly.py`: 26 unit tests for input parsing
- All 87 tests pass (54 existing + 26 new + 7 read_depth)
