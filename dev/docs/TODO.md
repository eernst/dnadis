# TODO

## Implement subcommands

* [ ] There are an overwhelming number of command line options available, some of which only apply to one or the other synteny mode. Look into the potential benefit of splitting these two modes of operation into subcommands, opening up the possibility of partitioning other use cases in the same way.

* [ ] We should implement a "clean" subcommand that can clean up residual executorlib_cache directories and any other intermediates (introduce a --keep-intermediates flag for the other subcommands to retain them in the first place). Consider whether we need to retain anything in the *_classification output subdirectory by default.

## Plotting improvements

### 1. Assembly Overview

#### Assembly Comparison Summary table

* [ ] Combine "Chr assigned", "Chr unassigned", and "Chimeric" columns a single composite "assigned/unassigned/chimeric" column (splitting the header into three lines to fit a narrower width) and join the values with a " / ".

* [ ] Ensure that we are only showing chrC for organisms that have a plastid genome, i.e. plants, algae, etc.

* [ ] Boldface the max value for each column. If there is a tie, bold all entries with the max value. For any composite column, apply this rule to the individual component pieces. Factor out any generalizable code for reuse with other tables.

### 2. Chromosome Assignment 

#### Overview

* [ ] Use identity rather than reference coverage for fill intensity.

* [ ] Use dark reference color for the small numbers for greater visibility.

#### Completeness Table

* [ ] Introduce separate tabs for each reference (sub)genome, for example, for the A+P+T reference, we'd have Completeness Tbl. A, Completeness Tbl. P, Completeness Tbl. T as separate tabs and tables. 

* [ ] Naturally sort the chromosome order in the tables.

### 3. Synteny

* [ ] Ribbons are not drawn between all chrom-assigned contigs and their nearest neighbor to the left. There seem to be two potential cases/causes:
  1. In some cases, the assembly neighbor to the immediate left does not have an assigned contig for a particular chromosome, and thus there is no underlying synteny information to draw ribbons from. This happens because we don't perform all-vs-all alignments or "rescue" these cases by doing a supplementary alignment of that individual chromosome with no partner to the immediate left to the nearest assembly to the left that does possess the query chromosome. Example: Ref.T chr10. Let's consider implementing one of these resolutions.
  2. Almost no ribbons are drawn between la0028 and its neighbor la0077. We need to investigate this case.

### 8. Contamination Comparison

* [ ] Add a top contaminants pie chart and top 5 table beneath it

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

* [ ] Add BUSCO runs on segregated datasets - the chromosome-assigned contigs and all other debris

* [ ] Aggregate compleasm and BUSCO results for the multiassembly comparison report as a gt table

  * [ ] Show classification category output numerically in the table row for each assembly

  * [x] Add compleasm S/D/F/I/M columns to `comparison_summary.tsv` (numeric counts + percentages)

  * [x] Plot in the comparison Rmd as a horizontal 100% stacked bar with the default BUSCO/compleasm color scheme

## Performance

* [ ] Add bgzip and indexing by default of pipeline outputs (FASTA, GFF3).

## Mammalian, other non-plant genome support

* [ ] The pipeline could be particularly useful in long-read survey sequencing of cancer genomes, where polyaneuploidy and large-scale rearrangements can be common.

* [ ] Ensure we have no hard organismal assumptions - non-plants won't have a plastid organellar genome. We should either bundle a few specific rDNA references or maybe extract Infernal/Rfam eukaryotic references and use those.