# TODO

## Improve Chromosome Synteny plot

* [ ] Ribbons are not drawn between all chrom-assigned contigs and their nearest neighbor to the left. There seem to be two potential cases/causes:
  1. In some cases, the assembly neighbor to the immediate left does not have an assigned contig for a particular chromosome, and thus there is no underlying synteny information to draw ribbons from. This happens because we don't perform all-vs-all alignments or "rescue" these cases by doing a supplementary alignment of that individual chromosome with no partner to the immediate left to the nearest assembly to the left that does possess the query chromosome. Example: Ref.T chr10. Let's consider implementing one of these resolutions.
  2. Almost no ribbons are drawn between la0028 and its neighbor la0077. We need to investigate this case.

## compleasm evals
[compleasm](https://github.com/huangnengCSU/compleasm) support added in phase 17 (`--compleasm-lineage`). Runs on `chrs.fasta` and `non_chrs.fasta` (debris + unclassified + contaminants). Results included in `comparison_summary.tsv`.

* [x] Add compleasm runs on segregated datasets - the chromosome-assigned contigs and all other debris

* [ ] Add BUSCO runs on segregated datasets - the chromosome-assigned contigs and all other debris

* [ ] Aggregate compleasm and BUSCO results for the multiassembly comparison report as a gt table

  * [ ] Show classification category output numerically in the table row for each assembly

  * [x] Add compleasm S/D/F/I/M columns to `comparison_summary.tsv` (numeric counts + percentages)

  * [x] Plot in the comparison Rmd as a horizontal 100% stacked bar with the default BUSCO/compleasm color scheme

## Mammalian, other non-plant genome support

* [ ] The pipeline could be particularly useful in long-read survey sequencing of cancer genomes, where polyaneuploidy and large-scale rearrangements can be common.

* [ ] Ensure we have no hard organismal assumptions - non-plants won't have a plastid organellar genome. We should either bundle a few specific rDNA references or maybe extract Infernal/Rfam eukaryotic references and use those.