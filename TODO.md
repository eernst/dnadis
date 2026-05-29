# TODO

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
