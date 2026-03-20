# TODO

## Consider snakemake

* [ ] executorlib is pretty flaky - moving to snakemake could simplify installation (as it handles conda dependencies in a more user-transparent fashion), simplify SLURM cluster execution, and support more platforms.

## Implement subcommands

* [ ] There are an overwhelming number of command line options available, some of which only apply to one or the other synteny mode. Look into the potential benefit of splitting these two modes of operation into subcommands, opening up the possibility of partitioning other use cases in the same way.

* [ ] We should implement a "clean" subcommand that can clean up residual executorlib_cache directories and any other intermediates (introduce a --keep-intermediates flag for the other subcommands to retain them in the first place). Consider whether we need to retain anything in the *_classification output subdirectory by default.

## Plotting improvements

* [ ] Add branding to the top of the report index/nav at left

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

* [ ] Our current layout is cramped. Need to make a few modifications:
  * [ ] Split the separate references into their own tabs titled by the reference, e.g. riparian plots for "Subgenome A", "Subgenome P", and "Subgenome T" would each be displayed in separate tabs, and separate PDFs should be produced for each.
  * [ ] Transpose the plots so that the chromosomes are on the x-axis and the assemblies are on the y-axis. The chrs should be naturally ordered from e.g. chr1...chr21, left to right.

* [ ] Currently, no ribbons are drawn for chromosome fragments (those with the _f# suffix). We need to handle these cases as well.

### 8. Contamination Comparison

* [x] Top Taxa tables: show binomial (genus + species) in the species table, add percentage labels in bars, hover tooltips showing full name and per-taxon assembly lists, ellipsis clipping for long names.

### Rename sections

* [ ] Drop the "Comparison" suffix from all sections titles

* [ ] Chromosome Synteny -> Macro Synteny

* [ ] Organelle -> Organelles

* [ ] Contamination -> Contaminants

## Better handling of unidentified contaminants

* [ ] I'm seeing frequent cases of what look like bacterial contaminant genomes in plant assemblies failing to be identified as contaminants with some potential knock-on effects in reporting and analysis - are they being bundled into non-chromosomes for the compleasm evals? We should consider whether we can flag these as unknown contaminants, provided they have no alignment-based similarity to the chromosome-assigned contigs.

## compleasm evals
[compleasm](https://github.com/huangnengCSU/compleasm) support added in phase 17 (`--compleasm-lineage`). Runs on `chrs.fasta` and `non_chrs.fasta` (debris + unclassified + contaminants). Results included in `comparison_summary.tsv`.

* [x] Add compleasm runs on segregated datasets - the chromosome-assigned contigs and all other debris

* [ ] Aggregate compleasm results for the multiassembly comparison report as a gt table

  * [ ] Show classification category output numerically in the table row for each assembly

  * [x] Add compleasm S/D/F/I/M columns to `comparison_summary.tsv` (numeric counts + percentages)

  * [x] Plot in the comparison Rmd as a horizontal 100% stacked bar with the default BUSCO/compleasm color scheme

## Report infrastructure

* [x] Factor shared R setup (packages, fonts, theming, OI palette, classification colors, helpers) into `output/reports/common.R` and shared CSS into `output/reports/common.css`, sourced by both report templates.

* [x] Rename `unified_report` to `assembly_report` for clarity; move templates into `output/reports/` subdirectory.

* [x] Replace all hardcoded bar colors and CSS gradients with `bar_cell()` helper using Okabe-Ito palette.

* [x] Add `--comparison-name` CLI argument and register all `--*-name` args in TOML config schema.

## Performance

* [ ] Add bgzip and indexing by default of pipeline outputs (FASTA, GFF3).

## Mammalian, other non-plant genome support

* [ ] The pipeline could be particularly useful in long-read survey sequencing of cancer genomes, where polyaneuploidy and large-scale rearrangements can be common.

* [ ] Ensure we have no hard organismal assumptions - non-plants won't have a plastid organellar genome. We should either bundle a few specific rDNA references or maybe extract Infernal/Rfam eukaryotic references and use those.
