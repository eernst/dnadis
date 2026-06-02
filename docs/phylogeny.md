# Phylogeny: per-leaf single-copy BUSCO accounting

The species-tree pipeline in `dnadis/phylogeny/` builds a leaf for every
`(assembly, reference_subgenome, query_subgenome)` combination and chooses
the BUSCOs that are single-copy *within each leaf*.  Compleasm's global
`Complete`/`Duplicated` label is **not** what governs that decision — it is
recomputed per leaf from the rows in `full_table_busco_format.tsv`.

This document explains how the same BUSCO appears differently in each
ploidy combination, so the per-leaf accounting matches biological intent.

## Leaf identity

A leaf is identified by:

- **assembly**: name of the query (or `reference`)
- **ref_subgenome**: suffix of the assigned reference chromosome name
  (e.g. `chr3A` → `A`, `chr3P` → `P`); `None` if the reference is monoploid
- **query_subgenome**: within-query duplication group inferred by
  `classifier.infer_query_subgenomes()` (primary haplotype = `None`,
  secondary = `B`, etc.)

The leaf label is `assembly[_ref_sg][_query_sg]`.  Examples:

| assembly | ref_sg | query_sg | label        |
|----------|--------|----------|--------------|
| asm1     | None   | None     | `asm1`       |
| asm1     | `A`    | None     | `asm1_A`     |
| asm1     | None   | `B`      | `asm1_B`     |
| asm1     | `A`    | `B`      | `asm1_A_B`   |

## How compleasm labels rows under each ploidy combination

For one BUSCO gene, the `full_table_busco_format.tsv` produced by
compleasm-on-chrs.fasta has one row per detected hit (status set per
*global* count: `Complete` = 1 hit, `Duplicated` = ≥2 hits, `Fragmented` =
partial, `Missing` = none).  The table below shows what dnadis sees for a
single BUSCO under each input topology, and which leaves count it as
single-copy.

| Input topology | Hits in compleasm output | Compleasm global status | Per-leaf hit counts | Single-copy in which leaves? |
|---|---|---|---|---|
| Collapsed haploid / diploid query, monoploid reference | 1 hit on the chrom set | `Complete` | `(asm, None, None)`: 1 | ✓ in the single leaf |
| Polyploid reference, haploid/collapsed query (chr1A *and* chr1P present) | 2 hits — one on the A-mapped contig, one on the P-mapped contig | `Duplicated` | `(asm, A, None)`: 1, `(asm, P, None)`: 1 | ✓ in **both** leaves |
| Monoploid reference, polyploid query (primary + alt haplotypes on the same chr) | 2 hits — one on each haplotype contig | `Duplicated` | `(asm, None, None)`: 1, `(asm, None, B)`: 1 | ✓ in **both** leaves |
| Polyploid reference *and* polyploid query | 4 hits — one per (ref_sg × query_sg) combination | `Duplicated` | `(asm, A, None)`: 1, `(asm, A, B)`: 1, `(asm, P, None)`: 1, `(asm, P, B)`: 1 | ✓ in **all four** leaves |
| Genuine duplication within one subgenome (e.g. tandem duplicate on chr3A) | ≥2 hits, all on contigs assigned to the same subgenome | `Duplicated` | `(asm, A, None)`: 2+ | ✗ in that leaf (excluded from supermatrix for *that* gene) |
| Subgenome-specific loss (present in A, missing from P) | 1 hit on the A-mapped contig | `Complete` | `(asm, A, None)`: 1, `(asm, P, None)`: 0 | ✓ in A only |
| Fragmented (partial hit, low score) | hit row with `Fragmented` status | `Fragmented` | excluded at parse time | ✗ everywhere — translated protein is unreliable |
| Missing | no hit row beyond the bare `Missing` entry | `Missing` | excluded at parse time | ✗ everywhere |

## Why this matters

If `single_copy_genes()` filtered on `status == "Complete"` alone, every
BUSCO that persists as one homeolog per ancestral subgenome would be
discarded — because compleasm marks the whole row set `Duplicated` when
≥ 2 copies are found, regardless of where they sit.  In an allopolyploid
with *n* well-retained ancestral subgenomes, most conserved single-copy
ancestral genes appear *n* times in compleasm's output and are therefore
all global-Duplicated; the share of `Duplicated` BUSCOs can exceed 90 %.
Filtering them out leaves single-digit per-leaf completeness and the
threshold (`--phylo-min-busco-completeness`, default 50 %) drops nearly
every leaf.  Recomputing per-leaf counts fixes this without changing the
threshold's semantics: it still means "≥ 50 % of the lineage's BUSCOs
present exactly once in this leaf".

## Supermatrix construction

A BUSCO contributes a row to the supermatrix only when it is single-copy
in **every retained** leaf.  Under the matrix above:

- Globally `Complete` BUSCOs are usable whenever the one hit lands on a
  retained leaf and that leaf has no other copies elsewhere.
- Globally `Duplicated` BUSCOs are usable when every retained leaf has
  exactly one copy.
- Globally `Fragmented`/`Missing` BUSCOs are never usable.

The protein sequence for each single-copy hit is looked up in compleasm's
`translated_protein.fasta` by `(busco_id, contig, start, end)`, aligned
per gene with MAFFT L-INS-i, trimmed with trimAl `-gappyout`, and the
trimmed alignments are concatenated horizontally into `supermatrix.faa`
with a consistent leaf order across genes.  IQ-TREE then runs on the
supermatrix; SH-aLRT/UFBoot branch supports appear in that order in the
output Newick (`species_tree.treefile`).

## Rooting

The tree is rooted by `--phylo-outgroup`: `none` (unrooted, default),
`reference`, `auto` (most-divergent taxon — not biologically conclusive),
or a query assembly name.  When a polyploid leaf set shares one outgroup
name, all of its subgenome leaves are passed to IQ-TREE as a comma-joined
outgroup clade.

## Phylogeny-only outgroups

A distant taxon is often useful *only* to root the tree: it diverged too
far to align to the reference, so its composition, pairwise synteny, and
per-chromosome assignments are noise rather than signal.  Designate such
an assembly with `--phylo-outgroup-only NAME` (repeatable; each value may
be a comma-separated list).  A phylogeny-only outgroup:

- Runs a **minimal pipeline** — compleasm only.  No synteny, detection,
  classification, read depth, scaffolding, or per-assembly outputs are
  produced beyond its `compleasm/` directory.
- Is **excluded** from `comparison_summary.tsv`,
  `chromosome_completeness.tsv`, pairwise synteny, and the HTML report.
  Exclusion is structural — the assembly never enters the `results` list
  those outputs iterate.
- **Contributes leaves** to the species tree and **auto-roots** it (so a
  separate `--phylo-outgroup` is unnecessary; if both are given,
  `--phylo-outgroup-only` supersedes it).

Because the minimal path has no classifications, subgenome identity for a
phylogeny-only outgroup is read from its **own contig names** rather than
from reference assignment: contigs sharing an alpha suffix (chr1A / chr2A
→ `A`, chr1B → `B`) form one leaf per suffix
(`build_outgroup_only_leaves`).  A contig with no recognizable suffix —
the common case for a divergent outgroup — pools into a single leaf
labeled with the bare assembly name.  This differs from
`--phylo-outgroup NAME` on a *regular* assembly
(`build_outgroup_leaves`), which still runs the full pipeline and infers
subgenomes from ploidy and reference-assignment quality.

Requires multi-assembly mode and `--compleasm-lineage`; at least one
non-outgroup assembly must remain.
