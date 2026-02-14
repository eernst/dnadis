# Unified HTML Report — Implementation Plan

## Context

The pipeline currently generates individual standalone outputs: chromosome overview (PDF/HTML), classification summary bar (PDF/HTML), depth overview (PDF/HTML), and contaminant table (HTML). There's no single document that aggregates all results with summary statistics and classification tables. A unified HTML report will provide a comprehensive, human-readable summary of the analysis — combining all existing visualizations with new gt-based summary tables for assembly quality assessment and genomic content understanding.

**Key decisions** (from previous discussion):
- **R Markdown** template with `__PLACEHOLDER__` tokens (matches existing `.tmpl.R` pattern)
- **Always generated** when `--plot` is specified (no separate flag)
- **Keep R** for all plotting (no D3/JS rewrite)
- **gt package** for tables with inline gradient bars and sparklines

## Design: Interactive Widget Embedding + New Tables

The Rmd template **aggregates** existing interactive outputs plus **new gt summary tables**. Strategy for embedding existing plots as interactive widgets:

- Each existing `.tmpl.R` template saves its ggiraph/gt widget object as an `.rds` file (~2 extra lines per template: `saveRDS(widget, path)`)
- The Rmd loads `.rds` files and prints the widgets — R Markdown natively renders htmlwidgets (ggiraph) and gt objects inline
- `self_contained: true` in the Rmd YAML ensures all JavaScript, CSS, and assets are embedded in the single HTML file

**Plot embedding strategy**:
- **Chromosome overview**: Load saved ggiraph widget from `.chromosome_overview.rds` → interactive with hover tooltips
- **Classification summary bar**: Load saved ggiraph widget from `.classification_summary_bar.rds` → interactive
- **Depth overview**: Load saved ggiraph widget from `.depth_overview.rds` → interactive
- **Contaminant table**: Load saved gt object from `.contaminant_table.rds` → renders as styled HTML table

This avoids duplicating any plot code while preserving full interactivity. Standalone templates are changed minimally (just add `saveRDS()`). The Rmd acts as a report shell that loads pre-built widgets alongside new gt summary tables.

**`--plot-html` dependency**: The ggiraph widgets are only constructed inside `if (plot_html)` blocks in the existing templates (the interactive versions use separate ggplot construction with `geom_*_interactive()` layers — ~400 extra lines for chromosome overview). The `.rds` files are saved alongside the standalone HTML widgets, only when `--plot-html` is active. The `saveRDS()` calls go inside the existing `if (plot_html)` blocks, right next to `saveWidget()`.

**Fallback for `--plot` without `--plot-html`**: The unified report still generates with all gt summary tables (the main new content). For sections where the widget `.rds` doesn't exist, the chunk renders a brief note: "Run with --plot-html for the interactive version of this plot." The standalone PDF is always generated regardless. This keeps the change to existing templates minimal (~2 lines each) without restructuring the `if (plot_html)` blocks.

## Files to create

### 1. `final_finalizer/output/unified_report.tmpl.Rmd`

The main R Markdown template with all report sections.

## Files to modify

### 2. `final_finalizer/output/plotting.py` — Add `run_unified_report()`
### 3. `final_finalizer/output/chromosome_overview.tmpl.R` — Add `saveRDS()`
### 4. `final_finalizer/output/classification_summary_bar.tmpl.R` — Add `saveRDS()`
### 5. `final_finalizer/output/depth_overview.tmpl.R` — Add `saveRDS()`
### 6. `final_finalizer/output/contaminant_table.tmpl.R` — Save gt object as `.rds`
### 7. `final_finalizer/cli.py` — Call `run_unified_report()`
### 8. `refresh_plots.py` — Add Rmd support
### 9. `docs/output_formats.md` — Document new output

## Implementation order

1. Archive this plan to `dev/docs/unified-html-report-plan.md`
2. Add `saveRDS()` calls to all 4 existing templates
3. Create `unified_report.tmpl.Rmd`
4. Add `run_unified_report()` to `plotting.py`
5. Add CLI integration in `cli.py`
6. Iterate on remaining Rmd sections
7. Update `refresh_plots.py`
8. Update `docs/output_formats.md`
9. Run tests to verify no regressions
