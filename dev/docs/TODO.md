# TODO

## Consider snakemake

* [ ] executorlib is pretty flaky - moving to snakemake could simplify installation (as it handles conda dependencies in a more user-transparent fashion), simplify SLURM cluster execution, possibly make the pipeline more robust, and support more platforms. Need to look at pros/cons, but pretty sure this needs to happen. If it does, we'll want to use final_finalizer.py as a convenience wrapper around snakemake (and retain our current flags). This will be a big refactor.

## Implement subcommands

* [ ] There are an overwhelming number of command line options available, some of which only apply to one or the other synteny mode. Look into the potential benefit of splitting these two modes of operation into subcommands, opening up the possibility of partitioning other use cases in the same way.

* [ ] We should implement a "clean" subcommand that can clean up residual executorlib_cache directories and any other intermediates (introduce a --keep-intermediates flag for the other subcommands to retain them in the first place). Consider whether we need to retain anything in the *_classification output subdirectory by default.

## NUMT/NUPT detection

* [ ] Add NUMT and NUPT detection using determined organelle reference contigs for alignment and accepted thresholds of sequence retention to call. Add stats tables to assembly report and include avg # and length in the aggregate homologous chromsome stats in the comparison report.

## Plotting improvements

* [ ] Add branding to the top of the report index/nav at left

* [ ] Comparison report: show reference-relative aggregate homologous chromosome stats, e.g. a table with box plots of lengths for each chromosome and min/med/max text next to it, rDNA presence/absence/mixed as full/empty/half-filled circle

### 1. Assembly Overview

#### Assembly Comparison Summary table

* [x] Combine "Chr assigned", "Chr unassigned", and "Chimeric" columns a single composite "assigned/unassigned/chimeric" column (splitting the header into three lines to fit a narrower width) and join the values with a " / ".

* [x] Ensure that we are only showing chrC for organisms that have a plastid genome, i.e. plants, algae, etc.

### 2. Chromosome Assignment

#### Overview

* [x] Use identity rather than reference coverage for fill intensity.

* [x] Use dark reference color for the small numbers for greater visibility.

#### Completeness Table

* [x] Introduce separate tabs for each reference (sub)genome, for example, for the A+P+T reference, we'd have Completeness Tbl. A, Completeness Tbl. P, Completeness Tbl. T as separate tabs and tables.

* [x] Naturally sort the chromosome order in the tables.

### 3. Synteny

* [x] Our current layout is cramped. Need to make a few modifications:
  * [x] Split the separate references into their own tabs titled by the reference, e.g. riparian plots for "Subgenome A", "Subgenome P", and "Subgenome T" would each be displayed in separate tabs, and separate PDFs should be produced for each.
  * [x] Transpose the plots so that the chromosomes are on the x-axis and the assemblies are on the y-axis. The chrs should be naturally ordered from e.g. chr1...chr21, left to right.

* [x] Currently, no ribbons are drawn for chromosome fragments (those with the _f# suffix). We need to handle these cases as well. — Resolved: fragments are `chrom_assigned` contigs and flow through the same pill + ribbon pipeline. They have macro_blocks, get pills, and get ribbons.

* [x] Multiple query subgenomes mapping to the same reference (e.g. chr2A_c1, chr2A_c2) are currently laid out side-by-side in the plot (good), but overlap (bad) making it hard to distinguish the multiple pills, and making the ribbons hard to read. Can we detect this situation (where any of the query assemblies have multiple subgenomes mapped to a single reference) and introduce extra horizontal space between the chromosome columns in that specific plot, so that the multiple query chromosomes won't overlap? We should also horizontally shift the centers of them so that the center of the composite (chr2A_c1, chr2A_c2) would align with another query assembly above or below with just one homologous query chromosome. In actuality, we should be horizontally centering all pills on each reference column, but if there are multiple pills in the same ref column, the composition of the multiple pills should be centered.

* [ ] We might actually need to do all vs all alignments within the set of contigs across all assemblies assigned to each reference subgenome. This would allow for more complete syntenic relationship tracking and plotting over rearrangements, potentially.

### 8. Contamination Comparison

* [x] Top Taxa tables: show binomial (genus + species) in the species table, add percentage labels in bars, hover tooltips showing full name and per-taxon assembly lists, ellipsis clipping for long names.

### Rename sections

* [x] Drop the "Comparison" suffix from all sections titles

* [x] Chromosome Synteny -> Macro Synteny

* [x] Organelle -> Organelles

* [x] Contamination -> Contaminants

## Better handling of unidentified contaminants

* [ ] I'm seeing frequent cases of what look like bacterial contaminant genomes in plant assemblies failing to be identified as contaminants with some potential knock-on effects in reporting and analysis - are they being bundled into non-chromosomes for the compleasm evals? We should consider whether we can flag these as unknown contaminants, provided they have no alignment-based similarity to the chromosome-assigned contigs.

## Query subgenome segmentation

* [x] Implement and document Gaussian mixture model segmentation of query subgenomes on the basis of reference alignment identity

## compleasm evals
[compleasm](https://github.com/huangnengCSU/compleasm) support added in phase 17 (`--compleasm-lineage`). Runs on `chrs.fasta` and `non_chrs.fasta` (debris + unclassified + contaminants). Results included in `comparison_summary.tsv`.

* [x] Add compleasm runs on segregated datasets - the chromosome-assigned contigs and all other debris

* [x] Aggregate compleasm results for the multiassembly comparison report as a gt table

  * [x] Show classification category output numerically in the table row for each assembly

  * [x] Add compleasm S/D/F/I/M columns to `comparison_summary.tsv` (numeric counts + percentages)

  * [x] Plot in the comparison Rmd as a horizontal 100% stacked bar with the default BUSCO/compleasm color scheme

## Report infrastructure

* [x] Factor shared R setup (packages, fonts, theming, OI palette, classification colors, helpers) into `output/reports/common.R` and shared CSS into `output/reports/common.css`, sourced by both report templates.

* [x] Rename `unified_report` to `assembly_report` for clarity; move templates into `output/reports/` subdirectory.

* [x] Replace all hardcoded bar colors and CSS gradients with `bar_cell()` helper using Okabe-Ito palette.

* [x] Add `--comparison-name` CLI argument and register all `--*-name` args in TOML config schema.

## Translocation-aware chromosome assignment

* [x] Replace raw-score-based reference assignment with reference span fraction (`qr_ref_span_bp / ref_length`) to avoid bias toward larger reference chromosomes in translocation cases. Implemented as a reassignment pass in `classify_all_contigs` using proportionate coverage.

## Reference chromosome filtering

* [ ] Filter reference sequences used for synteny assignment to exclude small scaffolds that can soak up spurious mappings and add noise. Precedence:

  1. **User-specified list** (`--ref-chr-list`): Only use the listed reference sequence IDs for chromosome assignment. Highest priority, explicit control.

  2. **Naming convention detection**: If no list provided, test reference sequence names against widely used chromosome naming conventions (chr1, Chr1, chromosome_1, etc. — extend `is_nuclear_chromosome()`). Keep only sequences matching recognized chromosome patterns.

  3. **Automatic gap detection**: If naming conventions don't match (e.g., scaffold_001, contig_123), sort reference sequences by length, detect the largest length gap, and retain only sequences above the gap. This handles references with a clear size discontinuity between chromosomes and unplaced scaffolds.

  4. **Keep-all flag** (`--ref-keep-all-contigs`): Opt-in to retain all reference sequences for assignment, overriding the above filters.

  Filtered-out reference sequences should still be present in the FASTA (for organelle/rDNA detection) but excluded from synteny-based chromosome assignment and the chromosome overview plot. This affects alignment target selection, assignment scoring, and all downstream outputs — needs careful integration across the pipeline.

## Incremental multi-assembly re-runs

* [ ] When re-running a multi-assembly job with a modified FOFN (added/removed/reordered assemblies), detect changes relative to existing cached results and compute a minimal set of jobs to update. This includes:
  * Per-assembly phases: skip assemblies whose inputs haven't changed and whose outputs are valid
  * Pairwise synteny: recompute only pairs affected by new/removed/reordered assemblies; invalidate pairs where the adjacent assembly changed
  * Rescue pairwise: regenerate based on updated assembly chain membership
  * Comparison report: always regenerate since row ordering and rescue logic depend on the full FOFN
  * Consider storing a manifest (FOFN hash, assembly list, per-assembly input checksums) alongside outputs to detect staleness

## Performance

* [ ] Add bgzip and indexing by default of pipeline outputs (FASTA, GFF3).

## Testing

* [ ] Need to develop test cases for classification, especially for complicated karyotypes and assembly issues involving multiple query subgenomes mapping to the same reference, segmentation of those query subgenomes, fragments vs complete, and the handling of contigs pre- and post-scaffolding.

## Mammalian, other non-plant genome support

* [ ] The pipeline could be particularly useful in long-read survey sequencing of cancer genomes, where polyaneuploidy and large-scale rearrangements can be common.

* [ ] Ensure we have no hard organismal assumptions - non-plants won't have a plastid organellar genome. We should either bundle a few specific rDNA references or maybe extract Infernal/Rfam eukaryotic references and use those.
