# Design Decisions

This document records significant design decisions, the alternatives considered,
and the rationale for the chosen approach.

---

## Scaffold naming preserves contig suffixes (2026-03-23)

**Decision**: Trivial scaffolds (single-contig groups and T2T contigs) use
`clf.new_name` as the scaffold/FASTA header name, preserving suffixes like
`_f1`, `_c2`, `_B`.

**Context**: The scaffolded FASTA (`chrs.fasta`, `chrs.A.fasta`) is used for
pairwise assembly-vs-assembly alignment. The comparison report's ribbon plot
matches pairwise macro_block contig names against pill names from the contig
summary. If the scaffold name differs from the contig summary name, ribbons
fail to draw.

**Alternatives considered**:
- Fix the comparison report's name mapping to handle the mismatch. Rejected
  because it would require maintaining a fragile scaffold→contig lookup, and
  the mismatch propagates to other downstream consumers of the FASTA.

**Consequence**: Pairwise alignment files now use the full renamed contig names
(e.g., `chr14A_f1` instead of `chr14A`), which means existing pairwise caches
must be cleared when re-running after this change.

---

## Non-T2T chromosome copies retained in scaffolded output (2026-03-23)

**Decision**: When a T2T contig exists in a scaffolding group, other
`chrom_assigned` contigs in the same group are emitted as their own trivial
scaffolds. Only debris contigs are dropped.

**Context**: Previously, the T2T path dropped all non-T2T contigs with the
rationale "the T2T contig already represents the complete chromosome." This
is wrong when the group contains multiple chromosome copies from different
query subgenomes (e.g., `chr19A_c1` and `chr19A_c2`) that ended up in the
same scaffolding group because subgenome segmentation didn't split them.

**Consequence**: `chrs.fasta` may now contain multiple contigs for the same
reference chromosome. This is correct — they represent distinct chromosomal
copies, not fragments of the same chromosome.

---

## Full-length classification uses length ratio rescue (2026-03-23)

**Decision**: `classify_full_length()` accepts an optional `length_ratio`
parameter (contig_len / ref_len). Contigs >= 80% of the reference chromosome
length are classified as full-length regardless of alignment-based
`ref_coverage`.

**Context**: In cross-species nucleotide mode comparisons (~60% identity),
alignment-based `ref_coverage` can be as low as 0.33 for a contig that is
92% of the reference length. The `ref_coverage` metric measures what fraction
of the reference is spanned by alignments — not how much of the chromosome
the contig represents. At low identity, many regions simply don't align.

**Alternatives considered**:
- Lower the `full_length_threshold` globally. Rejected because it would
  weaken the classification for within-species comparisons where ref_coverage
  is a reliable signal.
- Use a different coverage metric. The length ratio is complementary to
  ref_coverage — it captures physical completeness while ref_coverage
  captures sequence-level similarity.

**Threshold choice**: 0.80 was chosen as a conservative default. A contig
that is 80% of the reference length is unlikely to be a "fragment" in any
biological sense. The threshold could be made configurable via CLI if needed.

---

## GMM-based query subgenome inference with paired validation (2026-03-23)

**Decision**: Query subgenome inference uses a 1D Gaussian Mixture Model
fitted to alignment identities across all multi-copy chromosomes, with
model selection validated by per-chromosome pairing consistency.

**Previous approach**: Per-chromosome gap detection with an adaptive
threshold based on `k * std_dev` of primary-copy identities, floored by
`0.5 * (1 - mean_identity)`. This failed for cross-species comparisons
because the divergence-scaled floor (e.g., 0.16 at 68% mean identity) was
higher than the actual inter-subgenome gaps (0.06-0.09).

**New approach (3 stages)**:

1. **GMM fitting**: Fit 1D GMMs for k=1..3 using pure-Python EM. Select
   initial k by BIC (Bayesian Information Criterion).

2. **Paired validation**: BIC can be too conservative for small sample
   sizes (n~36). If BIC picks k=1, also try k=2 and validate using
   per-chromosome pairing: for chromosomes with exactly 2 copies, each
   copy should land in a different cluster. Accept the split if >= 70% of
   testable chromosomes are consistent AND cluster means are separated by
   >= 3% identity.

3. **Rank rescue**: For chromosomes where all copies land in the same
   cluster (borderline cases), rank-assign by identity (highest -> primary,
   next -> B, etc.) when the global split is already validated.

**Why pure Python**: The GMM operates on ~40-60 identity values per
assembly. EM converges in a handful of iterations on this scale. Adding
scipy or sklearn as a dependency would be disproportionate.

**Why not Otsu's method**: Otsu assumes equal variance and is inherently
binary (k=2 only). Extending to k>1 requires exhaustive search. GMM
naturally handles unequal variances and extends to k=3+ with BIC-based
selection. On this data, both find the same clusters, but GMM is more
principled and extensible.

**Why paired validation over BIC alone**: BIC's penalty term
(`num_params * log(n)`) is strong for small n. With 36 data points, BIC
favours k=1 even when the k=2 GMM finds well-separated clusters
(means 0.686 vs 0.614). The paired structure across chromosomes provides
strong biological evidence: if 17/18 chromosomes have their two copies
landing in different clusters, the split is real.

**Minimum separation threshold (3%)**: Prevents splitting when clusters
are very close (e.g., allelic variation within a single subgenome). BIC
can prefer k>1 for marginally better fit even when cluster means differ
by only 1-2%. The 3% floor ensures only biologically meaningful
divergence triggers a split.

---

## Permissive minimap2 chaining in nucleotide mode

**Decision**: Nucleotide mode always uses permissive minimap2 chaining
parameters (`--max-chain-skip=300 -z 10000,1000 -r 50000`).

**Context**: Even for byte-for-byte identical sequences, minimap2's
conservative default parameters split alignments at complex repetitive
regions (long homopolymer runs, tandem repeats). For chromosome
classification and structural composition visualisation, these kb-scale
gaps are noise that obscures megabase-scale patterns.

**Consequence**: Produces clean chromosome-scale blocks suitable for
compositional analysis. Downstream filtering (identity thresholds,
minimum alignment length, gate filtering) prevents spurious assignments.

---

## Gate-based chromosome assignment (AND logic)

**Decision**: Contigs must satisfy ALL criteria for chromosome assignment.
In protein mode: >= 5 segments, >= 3 genes, >= 50kb span, >= 20% span
fraction. In nucleotide mode: >= 1 segment, >= 50kb span, >= 20% span
fraction.

**Context**: Single conserved genes or repetitive sequences can produce
isolated hits that don't indicate chromosome-level homology. AND logic
ensures multiple independent lines of evidence converge.

**Why nucleotide mode has lower segment threshold**: Perfect full-length
alignments often produce very few segments (even 1) because they are
continuous matches. Protein hits are inherently fragmented.

---

## Interval tree for overlap detection in chain parsing

**Decision**: Uses `intervaltree` package for O(n log n) overlap
detection instead of O(n^2) pairwise comparison.

**Context**: Large assemblies with many protein alignments can have
tens of thousands of hits. The interval tree makes chain parsing
tractable for these cases.

**Consequence**: `intervaltree` is a hard dependency (no fallback).

---

## Two-gate contaminant filtering

**Decision**: Contaminant classification requires BOTH a centrifuger
score >= 1000 AND coverage >= 0.50.

**Context**: Single-gate filtering produces false positives. Conserved
genes (e.g., histones, ribosomal proteins) match non-target taxa with
high score but low coverage. Requiring coverage >= 50% ensures the
match represents genuine contamination, not just a few conserved genes.

---

## Depth caching with metadata validation

**Decision**: Read depth results are cached with metadata (query MD5,
window size, etc.) and automatically invalidated if inputs change.

**Context**: Read alignment and depth calculation are the most
expensive optional steps. Caching avoids redundant computation across
re-runs, but stale caches from changed inputs would produce incorrect
results.

---

## Atomic gzipped file writes via temp-and-rename

**Decision**: Gzipped output files (PAF, etc.) are written to a `.tmp`
file first, then atomically renamed on success. The temp file is cleaned
up on any exception.

**Location**: `utils/io_utils.py`

**Context**: In HPC/SLURM environments, jobs can be killed at any time
(timeout, OOM). A truncated `.paf.gz` file left on disk would be treated
as a valid cached result on the next run, silently corrupting downstream
analysis.

---

## Collinearity score: bp-weighted consecutive concordance

**Decision**: Compute a 0.0–1.0 collinearity score measuring what
fraction of aligned bp maintains reference coordinate order when chains
are sorted by query midpoint. Strand-aware: forward chains should have
increasing ref coordinates, reverse chains decreasing.

**Location**: `alignment/chain_parsing.py`

**Context**: Collinearity provides a confidence signal orthogonal to
identity or coverage. Highly ordered evidence strengthens confidence
(collinearity >= 0.9 with >= 5 chains can upgrade medium to high
confidence), while disordered evidence flags potential chimeras or
misassemblies.

---

## Top-K chain scoring for reference assignment

**Decision**: Score-based reference assignment uses the sum of top-K
chains (default K=3, configurable via `--assign-chain-topk`) rather
than all chains.

**Location**: `alignment/chain_parsing.py`

**Context**: Low-quality chains (small, low-identity) can dilute the
assignment score. In polyploid genomes, secondary homeologs produce many
small chains that could overwhelm the primary assignment if all chains
were summed equally. Top-K scoring focuses on the strongest evidence.

---

## Off-target cluster detection with gene density calibration

**Decision**: Rearrangement candidates are detected by identifying
off-target synteny clusters whose gene density exceeds a threshold
calibrated to the genome's own gene density (default: 10% of median
on-target density).

**Location**: `classification/classifier.py`

**Context**: A fixed gene density threshold would fail across genomes
with vastly different gene densities. Auto-calibrating to
`density_frac × median_on_target_density` ensures that gene-poor
regions (pericentromeres) are still detected, since even gene-poor
chromosomal content has higher density than scattered paralogous noise.

---

## Preference for mm2plus over minimap2

**Decision**: The tool checks for `mm2plus` first and falls back to
`minimap2`. The choice is cached globally for the process lifetime.

**Location**: `alignment/external_tools.py`

**Context**: mm2plus is a drop-in replacement with additional features.
Preferring it when available allows users to benefit from improvements
without any code or configuration changes.

---

## Early chromosome exclusion from BLAST phases

**Decision**: After synteny analysis (phase 2), contigs identified as
chromosome-assigned are excluded from the BLAST-based detection phases
(organelle, rDNA). A filtered FASTA is written containing only
non-chromosome contigs, and BLAST runs against this smaller input.

**Location**: `cli.py` (between phase 2 and phases 3-4)

**Context**: Chromosome contigs are large and numerous. Running BLAST
on the full assembly would be slow and produce irrelevant hits (e.g.,
organelle-derived insertions in nuclear chromosomes). Excluding them
significantly speeds up detection phases.

---

## Command injection prevention via regex validation

**Decision**: User-provided extra arguments (e.g., `--miniprot-extra-args`,
`--minimap2-extra-args`) are validated against a regex that rejects
shell metacharacters: `` [;&|`$(){}[\]<>!\\'"#] ``.

**Location**: `alignment/external_tools.py`

**Context**: Defence-in-depth measure. Arguments are already passed via
`shlex.split()` and `subprocess` list mode (not shell=True), but the
regex provides an additional layer of protection and clearer error
messages if a user accidentally includes shell syntax.

---

## TOCTOU fix: single stat() call for file validation

**Decision**: `file_exists_and_valid()` uses a single `os.stat()` call
to check both existence and non-zero size, rather than separate
`exists()` + `getsize()` calls.

**Location**: `utils/io_utils.py`

**Context**: In concurrent/distributed environments (SLURM with
multiple assemblies), a file could be deleted or truncated between an
`exists()` check and a subsequent `getsize()` check. A single `stat()`
call is atomic with respect to the filesystem.

---

## Stderr tee-writing for HPC diagnostics

**Decision**: stderr is duplicated to both the console and a structured
log file via a custom `_TeeWriter` class.

**Location**: `cli.py`

**Context**: In HPC environments, subprocesses (minimap2, BLAST, etc.)
write diagnostics to stderr, and unhandled Python exceptions produce
tracebacks on stderr. Without tee-writing, this output is lost when
jobs run non-interactively. The tee captures everything alongside the
structured `[info]`/`[warn]`/`[error]` log.

---

## Subgenome tie-breaking uses identity rank only (2026-03-23)

**Decision**: When per-chromosome subgenome assignment is ambiguous (both
copies land in the same GMM cluster), rank-assign by identity rather
than incorporating additional metrics (collinearity, synteny score, GC).

**Alternatives considered**:
- Multivariate GMM using identity + GC + synteny score. Rejected because
  GC differences between subgenomes are tiny (0.2-0.6%) and collinearity
  gaps are inconsistent in sign across chromosomes — neither provides a
  reliable directional signal.
- K-mer distance or tandem repeat fingerprinting for assembly-intrinsic
  subgenome signals. Scoped as future work (see
  `dev/docs/kmer-subgenome-inference-plan.md`) — would add dependencies
  and complexity for marginal improvement in the current pipeline.

**Rationale**: Identity rank is the simplest approach that's consistent
with the global pattern. For genuinely ambiguous cases where the
per-chromosome identity gap is very small (e.g., 2%), every
reference-based metric is near-tied, so additional signals wouldn't
help. The rank-rescue only fires when the global split is already
validated (≥70% of chromosomes consistent), so the direction of
assignment is well-established.

---

## BAM/CRAM format auto-detection via header inspection

**Decision**: Read input format is detected by running `samtools view -H`
and checking for `@SQ` header lines. Aligned BAM/CRAM (with `@SQ`)
skips the alignment step; unaligned BAM (uBAM, no `@SQ`) is treated
like FASTQ.

**Location**: `analysis/read_depth.py`

**Context**: Users may provide pre-aligned BAMs from upstream pipelines
or unaligned uBAMs from PacBio. Auto-detection eliminates the need for
a `--reads-format` flag and avoids re-aligning already-aligned data.
