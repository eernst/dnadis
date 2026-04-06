# rDNA annotation
* Add a routine to annotate rDNA loci, and cache the results as with other potentially long-running analyses.
* We should consider how to derive a consensus 45S nuclear rDNA repeat. 
* Display annotated regions in the contig composition panel to aid in understanding assembly quality, completeness and genome structure.
* Output a properly-formatted bed12 file with rDNA annotations, sorted, bgzipped and tabix-indexed.
* Refer to /home/eernst/scriptorium/annotation/infernal.sh

# Chromosome synteny plot
Consider using integrating GENESPACE, or maybe directly from the synteny data we have using ggalluvial. This will be a separate full-page plot, 7.2 in. wide and up to 10 in. tall.

# Reference-based scaffolding
We currently detect 5' and 3' telomere sequence, and there are frequently cases of split-arm assemblies that are made obvious by our visualization. We also have the relative positions of each query contig on the nearest reference genome, so we could consider producing an edited version of the assembly FASTA that includes not just new names for contigs (ref. chromosome assigned, organelles), but also scaffolded contigs joined by a fixed length of Ns where we have good evidence. Need to consider whether to include smaller contigs in the scaffolding process (in which case we might want to use ragtag with just the target reference chromosome and contigs with alignments to it), or to only join chromosome-assigned contigs.

# Unified HTML report
* We should generate a single, unified HTML report with all current plots and additional basic statistics and classification tables. For the tables, we should look into the continued use of the gt package to embed plots where helpful to provide quick graphical feedback. The report should be comprehensive, but not riddled with unnecessary technical detail. 

We should try to maintain consistent styling throughout where possible, i.e. using the same set of typefaces and semantic use of text size throughout. Use the r-dataviz-expert subagent for design and implementation assistance.
