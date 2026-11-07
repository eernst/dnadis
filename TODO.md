# TODO

## Classification for extra copies of a chromosomal region

On a chromosome split into query subgenome sets, a fragment whose reference
footprint is already covered by every full-length copy, but which is not
near-identical to any of them (so phase-5 chromosome debris misses it), is
flagged `subgenome_basis = extra_copy` and left without a set. It stays
`chrom_fragment`. These contigs mix several origins that read depth partly
separates:

- collapsed duplicates or amplified regions (depth ratio ~2×)
- divergent haplotigs of one homeolog (depth ratio ~0.5×)
- segmental duplications or an extra homeolog of a region (depth ratio ~1×)

Consider a dedicated classification (e.g. `chrom_extra_copy`) with its own
FASTA, report colour and comparison columns, annotated with depth ratio.
Related: the footprint-containment route to `chrom_debris` gates on identity
to the *reference* (`--chr-debris-min-identity`), which never passes against
a divergent reference; a containment rule based on coverage by the assembled
copies would catch these contigs regardless of reference distance.

## True contaminant (adapter / vector) detection phase

The centrifuger screen (phase 6) classifies contigs by taxonomy and is best
understood as **cobiont** detection — it flags sequence from organisms that
co-occur with the target (symbionts, commensals, gut/endophytic microbiota,
co-purified organisms). These are biologically real and not necessarily
artifacts.

Genuine **contaminants** — sequencing adapters, cloning/expression vectors,
spike-ins, and other engineered or technical sequence — are a different problem
and call for a different method (e.g. matching against UniVec / adapter
databases, exact/near-exact local alignment rather than k-mer taxonomic
classification).

Add a dedicated contaminant-detection phase that:

- Runs **before** the centrifuger cobiont screen, so adapter/vector sequence is
  removed (or flagged) before taxonomic classification.
- Uses an alignment-based methodology against vector/adapter references
  (e.g. UniVec) rather than taxonomic k-mer matching.
- Emits its own classification category (`contaminant`) and output FASTA/TSV,
  kept distinct from the `cobiont` category introduced when the centrifuger
  screen was renamed.

## Scaffolding results in the reports

With `--scaffold`, the comparison report's Scaffolding tab summarizes joins per
assembly, and both reports draw scaffold brackets or pills in their synteny
plots. Most other report sections still describe contigs:

- Nx curves, the contiguity heatmap and chromosome completeness use contig
  lengths. Scaffold-level versions (or a contig/scaffold toggle) would show
  what scaffolding changed.
- The per-assembly report has no scaffolding summary: joins per scaffold,
  RagTag location/orientation confidence, and the fragments left on their own
  with the reason (no reference placement, rDNA, sister-placed, no anchor).
- `comparison_summary.tsv` has no scaffolding columns (joined Mb, joins,
  scaffold N50, chromosomes in one scaffold).
- Joins between copy sets or against sister-assembly evidence are not flagged
  in either report.
