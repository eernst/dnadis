# K-mer Based Subgenome Inference — Future Work Plan

## Motivation

The current subgenome inference uses reference alignment identity as the
primary signal: contigs from subgenome A align to the reference at ~69%
identity while subgenome B aligns at ~61%. This works well when the
identity gap is clear (>3%), but fails for chromosomes where both copies
have similar divergence from the reference (e.g., both copies at ~0.62
identity with only a 2% gap). Additional reference-based metrics (collinearity, synteny score,
GC content) don't provide a consistent directional signal for
tie-breaking.

Assembly-intrinsic approaches — comparing contigs to *each other* rather
than to a reference — could provide an orthogonal signal. Contigs from
the same subgenome share more recent common ancestry and should be more
similar to each other (in k-mer composition, repeat landscape, etc.)
than to contigs from the other subgenome.

## Proposed Approaches

### 1. K-mer distance matrix + clustering

**Concept**: Compute pairwise k-mer distances between all chromosome-
assigned contigs, then cluster into subgenome groups.

**Implementation sketch**:
- For each chrom_assigned contig, compute a k-mer sketch (e.g., MinHash
  via sourmash, or FracMinHash for containment-based distance)
- Build a pairwise distance matrix across all multi-copy contigs
- Cluster into k groups (where k = number of expected subgenomes, from
  GMM identity-based inference or user-specified)
- Validate: contigs assigned to the same reference chromosome should
  land in different clusters

**Tools/dependencies**:
- sourmash (Python, conda-installable): FracMinHash sketches, pairwise
  compare, already supports genome-scale inputs
- Alternative: mash (C++, fast but less flexible)
- Alternative: pure Python with collections.Counter on short k-mers
  (k=21), but slow for chromosome-scale contigs

**Strengths**:
- No reference required — works for orphan assemblies
- Captures repeat landscape divergence (often the strongest subgenome
  signal in allopolyploids)
- Natural validation/complement to identity-based approach

**Weaknesses**:
- Computationally heavier: O(n²) pairwise comparisons for n contigs
- Requires additional dependency (sourmash or mash)
- May not work well for very recent allopolyploids where subgenomes
  haven't diverged in k-mer composition

**Estimated effort**: Medium. Sourmash has a Python API that makes
sketch computation and comparison straightforward. The clustering step
reuses the existing GMM infrastructure.

### 2. Tandem repeat fingerprinting via TRF

**Concept**: Run Tandem Repeat Finder on each contig, extract a
repeat-unit spectrum (lengths, GC, copy numbers), and cluster contigs
by repeat landscape similarity.

**Implementation sketch**:
- Run TRF on each chrom_assigned contig
- Extract features: distribution of repeat unit lengths, total tandem
  repeat fraction, dominant satellite unit sequences
- Compute pairwise distances based on repeat feature vectors
- Cluster as above

**Strengths**:
- Satellite repeats often diverge faster than coding sequences between
  subgenomes — potentially the strongest signal for recent polyploids
- Could detect subgenome-specific centromeric satellites

**Weaknesses**:
- TRF is slow on chromosome-scale contigs (minutes per contig)
- New external dependency
- Signal is concentrated in repeat-rich regions (centromeres,
  subtelomeres) — sparse for gene-rich chromosome arms
- Feature engineering required (what aspects of the repeat landscape
  to compare)

**Estimated effort**: High. TRF output parsing, feature extraction,
and distance computation are non-trivial. Less mature Python ecosystem
than k-mer tools.

### 3. Hybrid: k-mer validation of identity-based segmentation

**Concept**: Use k-mer distances not as the primary segmentation method,
but as a validation/rescue step after identity-based GMM inference.

**Implementation sketch**:
- Run identity-based GMM segmentation (current approach)
- For chromosomes where segmentation is ambiguous (both copies in same
  cluster, or low confidence), compute k-mer distances between the
  ambiguous contig and representative contigs from each subgenome
  cluster
- Assign the ambiguous contig to the cluster whose representatives are
  more similar by k-mer distance

**Strengths**:
- Only runs k-mer analysis on ambiguous cases (fast)
- Doesn't require changing the primary segmentation logic
- Minimal dependency: could use a lightweight pure-Python k-mer counter
  on just the ambiguous contigs + a few representatives

**Weaknesses**:
- Only helps with tie-breaking, doesn't improve the primary signal
- May not help if the ambiguous case is genuinely ambiguous (both
  subgenomes equidistant)

**Estimated effort**: Low-medium. Could be implemented as an optional
post-processing step.

## Recommendation

Start with **approach 3** (hybrid validation) as it has the lowest risk
and effort. If the identity-based GMM correctly segments most
chromosomes (typically 17/18 or better), k-mer validation only needs
to handle the residual cases. This avoids adding a heavy dependency for
a marginal improvement.

If there's demand for reference-free subgenome inference (e.g., for
orphan species), pursue **approach 1** (full k-mer clustering) as a
standalone mode, likely using sourmash as an optional dependency
(similar to how compleasm and centrifuger are optional).

**Approach 2** (TRF) is the most speculative and should only be pursued
if satellite repeat divergence proves to be a uniquely strong signal
that k-mers miss.

## Integration Points

- `classification/classifier.py:infer_query_subgenomes()` — add an
  optional k-mer validation step after GMM assignment
- New module: `detection/kmer_subgenome.py` or similar
- CLI flags: `--subgenome-kmer-validation` (enable hybrid mode),
  `--subgenome-kmer-only` (reference-free mode)
- Optional dependency: sourmash (for approach 1/3 with external tool)
  or pure Python (for approach 3 lightweight version)
