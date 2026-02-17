# Multi-Assembly Support — Phase 2 & 3: Comparative Analysis and Visualization

## Context

Phase 1 (complete, committed on `feat/multi-assembly`) refactored `main()` into `prepare_reference()` → `ReferenceContext` and `run_assembly()`, added multi-assembly CLI args (`--fofn`, `--assembly-dir`, `--output-dir`), and a dispatch loop in `main()`. Currently `run_assembly()` returns `None` and no cross-assembly comparison is produced.

Phase 2 adds a data layer: `run_assembly()` returns an `AssemblyResult` summary, and after all assemblies complete, comparison TSVs are written. Phase 3 adds a comparison HTML report with interactive plots rendered via an Rmd template.

---

## Phase 2: Comparative Analysis Data Layer

### 2.1 New dataclasses in `models.py`

- `ChromRefSummary`: Per-reference-chromosome summary for one assembly (ref_id, ref_length, n_contigs, total_assigned_bp, best_ref_coverage, is_full_length, has_both_telomeres, mean_identity)
- `AssemblyResult`: Full summary metrics from one assembly's finalization (contiguity, classification counts/bp, chromosome completeness, quality metrics, contamination, rDNA, organelles, read depth, per-ref detail, output file paths)

### 2.2 New file: `final_finalizer/output/comparison.py`

- `_compute_n50_l50()`: N50/L50 from contig lengths
- `build_assembly_result()`: Computes all metrics from classifications, qry_lengths, ev, etc.
- `write_comparison_summary_tsv()`: One row per assembly with 29 comparison columns
- `write_chromosome_completeness_tsv()`: Long-format (ref_id × assembly) matrix

### 2.3 Changes to `cli.py`

- `run_assembly()` return type changed from `None` to `AssemblyResult`
- `build_assembly_result()` call inserted after classification summary logging
- Multi-assembly dispatch loop collects results and writes comparison TSVs

## Phase 3: Comparative Visualization

### 3.1 New Rmd template: `comparison_report.tmpl.Rmd`

Report sections:
1. **Overview**: Summary table, classification distribution stacked bar, size comparison
2. **Chromosome Completeness**: Heatmap (ref × assembly, fill = coverage), completeness table
3. **Quality Comparison** (conditional): Identity, collinearity, GC deviation violin+box plots
4. **Contamination Comparison** (conditional): Bar chart and summary table
5. **Depth Comparison** (conditional): Depth distribution violin+box per assembly

### 3.2 New function: `plotting.py:run_comparison_report()`

Same pattern as `run_unified_report()`: check Rscript/rmarkdown, load template, substitute placeholders, render.

### 3.3 Wired up in `cli.py main()`

After writing comparison TSVs, if `--plot` and results exist, calls `run_comparison_report()`.

## Files created/modified

| File | Change |
|------|--------|
| `final_finalizer/models.py` | Added `ChromRefSummary` and `AssemblyResult` dataclasses |
| `final_finalizer/output/comparison.py` | **New**: `build_assembly_result()`, TSV writers |
| `final_finalizer/cli.py` | Changed `run_assembly()` return type; added result collection and comparison outputs |
| `final_finalizer/output/plotting.py` | Added `run_comparison_report()` |
| `final_finalizer/output/comparison_report.tmpl.Rmd` | **New**: Comparison report Rmd template |
| `tests/test_comparison.py` | **New**: 24 unit tests for Phase 2 functions |

## Verification

- All 111 non-integration tests pass
- 24 new tests cover N50/L50, build_assembly_result, comparison TSV writers, chromosome completeness matrix
