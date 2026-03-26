# TODO

## Current PR: Reporting Improvements

* [ ] Add branding to the top of the report index/nav at left

* [ ] Add a legend to the Chromosome Assignments compound plot (contig composition + identity panel). Should explain: blue/colored blocks = on-target synteny, red blocks = off-target synteny (potential rearrangement), pill background = contig extent, identity values = matches/alignment_length from best chain to assigned reference.

* [ ] Comparison report: show reference-relative aggregate homologous chromosome stats, e.g. a table with box plots of lengths for each chromosome and min/med/max text next to it, rDNA array presence/absence/mixed as full/empty/half-filled circle, rearrangements detected.

* [ ] Consider adding a flag to allow different sort orders for assemblies in the report. Current should be using fixed sort order by fofn input order, but users may prefer sorting by global identity to the reference, or other metrics.

* [ ] We might actually need to do all vs all alignments within the set of contigs across all assemblies assigned to each reference subgenome. This would allow for more complete syntenic relationship tracking and plotting over rearrangements, potentially.

---

## Future PRs

### Pipeline architecture

* [ ] **Snakemake migration**: executorlib is pretty flaky - moving to snakemake could simplify installation (as it handles conda dependencies in a more user-transparent fashion), simplify SLURM cluster execution, possibly make the pipeline more robust, and support more platforms. Need to look at pros/cons, but pretty sure this needs to happen. If it does, we'll want to use final_finalizer.py as a convenience wrapper around snakemake (and retain our current flags). This will be a big refactor.

* [ ] **Subcommands**: There are an overwhelming number of command line options available, some of which only apply to one or the other synteny mode. Look into the potential benefit of splitting these two modes of operation into subcommands, opening up the possibility of partitioning other use cases in the same way.

* [ ] **Clean subcommand**: Implement a "clean" subcommand that can clean up residual executorlib_cache directories and any other intermediates (introduce a --keep-intermediates flag for the other subcommands to retain them in the first place). Consider whether we need to retain anything in the *_classification output subdirectory by default.

* [ ] **CLI flag naming**: Review CLI flag names for clarity. `--min-span-bp` controls the minimum synteny block span for chromosome assignment gates, but the name doesn't convey this. Consider renaming to something like `--min-synteny-span` or `--gate-min-span`. This flag also implicitly controls the minimum detectable rearrangement size (since rearrangement detection operates on macro_blocks that passed this threshold).

* [ ] **Incremental multi-assembly re-runs**: When re-running a multi-assembly job with a modified FOFN (added/removed/reordered assemblies), detect changes relative to existing cached results and compute a minimal set of jobs to update. This includes:
  * Per-assembly phases: skip assemblies whose inputs haven't changed and whose outputs are valid
  * Pairwise synteny: recompute only pairs affected by new/removed/reordered assemblies; invalidate pairs where the adjacent assembly changed
  * Rescue pairwise: regenerate based on updated assembly chain membership
  * Comparison report: always regenerate since row ordering and rescue logic depend on the full FOFN
  * Consider storing a manifest (FOFN hash, assembly list, per-assembly input checksums) alongside outputs to detect staleness

### Classification and assignment

* [ ] **Fragmented reciprocal translocations**: The current reciprocal detection requires the partner ref to have zero assigned contigs. This fails when the assembly is fragmented and one fragment has already been independently assigned to the partner ref.

  **Possible approach**: Use the coordinate layout of synteny blocks on the main translocation contig to determine which fragments belong with which arm. For each fragment assigned to the same ref as the main contig, check whether its synteny blocks overlap the main contig's R1 region or R2 region. If they overlap the R2 region (the translocated arm), reassign the fragment to R2. This requires:
  - Detecting that the main contig spans two reference chromosomes (already visible in the per-assembly macro_blocks)
  - Partitioning the main contig's query coordinates into R1 and R2 zones
  - Checking each fragment's macro_block coordinates against these zones
  - Only reassigning when the overlap with the R2 zone significantly exceeds overlap with the R1 zone

* [ ] **Reference chromosome filtering**: Filter reference sequences used for synteny assignment to exclude small scaffolds that can soak up spurious mappings and add noise. Precedence:

  1. **User-specified list** (`--ref-chr-list`): Only use the listed reference sequence IDs for chromosome assignment. Highest priority, explicit control.

  2. **Naming convention detection**: If no list provided, test reference sequence names against widely used chromosome naming conventions (chr1, Chr1, chromosome_1, etc. — extend `is_nuclear_chromosome()`). Keep only sequences matching recognized chromosome patterns.

  3. **Automatic gap detection**: If naming conventions don't match (e.g., scaffold_001, contig_123), sort reference sequences by length, detect the largest length gap, and retain only sequences above the gap. This handles references with a clear size discontinuity between chromosomes and unplaced scaffolds.

  4. **Keep-all flag** (`--ref-keep-all-contigs`): Opt-in to retain all reference sequences for assignment, overriding the above filters.

  Filtered-out reference sequences should still be present in the FASTA (for organelle/rDNA detection) but excluded from synteny-based chromosome assignment and the chromosome overview plot. This affects alignment target selection, assignment scoring, and all downstream outputs — needs careful integration across the pipeline.

### Detection

* [ ] **NUMT/NUPT detection**: Add NUMT and NUPT detection using determined organelle reference contigs for alignment and accepted thresholds of sequence retention to call. Add stats tables to assembly report and include avg # and length in the aggregate homologous chromsome stats in the comparison report.

* [ ] **Read depth evidence for rearrangement confidence**: Use read depth uniformity at predicted rearrangement breakpoints to refine confidence. A real rearrangement should show uniform depth across the breakpoint; an assembly misjoin often shows a depth drop or spike at the junction.

  **Phase ordering**: Swap Phase 12 (rearrangement detection) and Phase 13 (read depth) so that depth data is available when rearrangement confidence is scored. Single-pass — no need for a two-pass approach.

  **Depth resolution**: Consider reducing mosdepth window size from 500bp to 100bp for finer breakpoint resolution. Estimated file size increase ~2-3× (not 5×) due to gzip compression efficiency on autocorrelated depth values.

  **Implementation**:
  - After depth calculation, for each rearrangement breakpoint candidate, extract depth in a ±5kb window around the junction from the mosdepth bed.gz
  - Compute depth uniformity metric (e.g., coefficient of variation in the breakpoint window vs flanking regions)
  - Uniform depth across breakpoint → upgrade confidence; depth anomaly → add caveat or downgrade
  - No dependency on `--keep-bam` for depth-based evidence

  **Opportunistic split-read analysis**: If `--keep-bam` is in use, additionally check for split/discordant reads spanning predicted breakpoints. This provides direct confirmation but is optional — depth uniformity alone is sufficient for confidence adjustment.

* [ ] **Unidentified contaminants**: Frequent cases of what look like bacterial contaminant genomes in plant assemblies failing to be identified as contaminants with some potential knock-on effects in reporting and analysis - are they being bundled into non-chromosomes for the compleasm evals? We should consider whether we can flag these as unknown contaminants, provided they have no alignment-based similarity to the chromosome-assigned contigs.

### Testing

* [ ] Develop test cases for classification, especially for complicated karyotypes and assembly issues involving multiple query subgenomes mapping to the same reference, segmentation of those query subgenomes, fragments vs complete, and the handling of contigs pre- and post-scaffolding.

* [ ] Develop synthetic test cases for chromosome assignment with known rearrangements: reciprocal translocations (balanced and unbalanced), Robertsonian translocations, inversions, whole-arm translocations, fusions, and fissions. Include fragmented assembly variants of each scenario. Verify assignment outcomes match biological expectations across scenarios.

### Performance

* [ ] Add bgzip and indexing by default of pipeline outputs (FASTA, GFF3).

### Broader scope

* [ ] **Mammalian/non-plant support**: The pipeline could be particularly useful in long-read survey sequencing of cancer genomes, where polyaneuploidy and large-scale rearrangements can be common. Ensure we have no hard organismal assumptions - non-plants won't have a plastid organellar genome. We should either bundle a few specific rDNA references or maybe extract Infernal/Rfam eukaryotic references and use those.
