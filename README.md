# dnadis

**Genome assembly post-processing tool for contig classification, homology assignment, and quality control** - the *de novo* assembly disambiguator

`dnadis` classifies contigs from a *de novo* genome assembly into biological categories (chromosomes, organelles, rDNA, contaminants, debris) using nucleotide or protein-anchored synteny with a reference, organelle/rDNA alignments, and taxonomic classification. Beyond classification, it evaluates assembly quality through BUSCO completeness scoring (via compleasm), syntenic block coverage, and alignment identity metrics, and produces rich interactive HTML reports for both individual assemblies and multi-assembly comparisons. Multi-assembly mode aggregates individual assessments for easy comparisons between various assemblies of the same individual, between different individuals of the same species, or even between multiple species.

## Features

- **Dual synteny modes**:
  - **Nucleotide mode** (default): Whole-genome nucleotide alignment (minimap2) for structural composition analysis
  - **Protein mode**: Protein-anchored synteny blocks (miniprot) for gene-level classification
- **Chromosome assignment** with gate-based filtering to prevent spurious assignments
- **Subgenome resolution** for polyploid genomes with multiple chromosome sets
- **Organelle detection** (chloroplast/chrC, mitochondrion/chrM) via BLAST
- **rDNA contig identification** via BLAST against reference rDNA
- **rDNA consensus building and annotation** (optional) with structure-based sub-feature detection (18S, 5.8S, 25S, ITS1, ITS2) using Infernal/Rfam
- **Contaminant screening** via centrifuger taxonomic classification
- **Debris detection** for assembly fragments (chromosome debris, organelle debris)
- **Chimera flagging** for contigs with evidence from multiple chromosomes
- **Orientation determination** for chromosome-assigned contigs
- **Contig renaming** to reference-based names following the scheme `chr<ref>(_<subgenome>)?(_c<copy>|_f<frag>)?` — for example: `chr1A` (full-length, single copy), `chr1A_B` (full-length, query subgenome B), `chr1A_f1`/`chr1A_f2` (chromosome fragments by descending length), `chr1A_c1`/`chr1A_c2` (rare duplicate full-length copies), `chr1A_B_f1` (fragment from subgenome B); non-chromosome contigs are named `contig_1`, `contig_2`, etc. by descending length (see [Contig Naming Scheme](#contig-naming-scheme))
- **BUSCO completeness evaluation** (optional) via compleasm, run separately on chromosome-assigned and non-chromosome contigs
- **Reference-guided scaffolding** (optional) producing chromosome-scale pseudomolecules with AGP output (uses RagTag if available, otherwise built-in scaffolder)
- **Read depth analysis** (optional) with automated downsampling and caching
- **Multi-assembly mode**: analyze multiple assemblies concurrently against a shared reference via `--fofn` (file-of-filenames TSV) or `--assembly-dir` (directory scan)
- **Cross-assembly comparison reports**: interactive HTML tables (via gt) aggregating classification and BUSCO completeness results across all assemblies
- **Interactive HTML reports** per assembly with chromosome overview, classification summary, read depth, and contaminant table (enabled by default; requires rmarkdown + pandoc). Reports are self-contained by default (single portable file); use `--no-self-contained-html` for HTML + companion `_files/` directory (faster rendering)
- **Distributed SLURM execution** (optional) via executorlib: compute-intensive phases submitted as individual SLURM jobs with automatic resource estimation

## Contig Naming Scheme

Chromosome-assigned contigs are renamed using the pattern `chr<ref>(_<subgenome>)?(_c<copy>|_f<frag>)?`, where `<ref>` is the reference chromosome identifier (e.g., `1A`, `3B`). Non-chromosome contigs are named `contig_1`, `contig_2`, etc., ordered by descending length.

### Chromosome name suffixes

| Suffix | Meaning | Example |
|--------|---------|---------|
| _(none)_ | Primary contig assigned to this reference chromosome — either the only copy, or the highest-identity copy when multiple query subgenomes are detected | `chr1A` |
| `_B`, `_C`, … | Secondary query subgenome label — the query assembly carries multiple homeologous copies that map to the same reference chromosome, resolved into distinct subgenomes by identity clustering. The primary (highest-identity) copy is always unsuffixed. | `chr1A_B` |
| `_c1`, `_c2`, … | Multiple full-length copies that could not be resolved into distinct subgenomes — the copies align to the reference at nearly identical identity levels, so the GMM clustering cannot distinguish them. This typically arises when both haplotypes of a phased assembly are provided in a single query FASTA, or when the assembler produces near-identical duplicate copies. Ordered by descending alignment identity; `_c1` is always the highest-identity copy. See [Limitations](#limitations-and-expectations) for details. | `chr1A_c1`, `chr1A_c2` |
| `_f1`, `_f2`, … | Chromosome fragments (contigs not classified as full-length), ordered by descending alignment identity to the reference | `chr1A_f1`, `chr1A_f2` |

Suffixes compose left-to-right: subgenome first, then copy/fragment. For example, `chr1A_B_f1` is the longest fragment of chr1A from query subgenome B.

A contig is classified as **full-length** when its syntenic coverage of the reference chromosome meets `--full-length-ref-coverage` (default: 0.70) and/or it has telomeres at both ends. The `_f` suffix therefore flags genuine assembly fragmentation rather than short contigs of arbitrary origin.

Subgenome labels (`_B`, `_C`, …) are inferred by clustering query contigs that all map to the same reference chromosome by sequence-identity similarity. The label letters are assigned in order of cluster size (largest cluster = A, or inherits the reference subgenome letter when the reference already carries one, e.g., `chr1A` → `chr1A_B` means the new subgenome differs from subgenome A).

The authoritative implementation is `dnadis/classification/classifier.py:generate_contig_names()`. See [Output Formats](docs/output_formats.md#contig_summarytsv) for the corresponding `contig` TSV column.

## Query Subgenome Segmentation

When multiple contigs from the query assembly map to the same reference chromosome — as is expected for polyploid genomes where two or more homeologous chromosome sets are present — `dnadis` automatically attempts to assign each contig to a distinct query subgenome. The result appears in the contig name as a subgenome suffix (`_B`, `_C`, …) or, when segmentation is not possible, as a copy suffix (`_c1`, `_c2`, …).

### How it works

Subgenome assignment is based on the observation that homeologous chromosomes from different ancestral subgenomes typically align to the same reference chromosome at detectably different sequence identity levels. The tool clusters contigs using a 1-D Gaussian Mixture Model (GMM) fitted to alignment identities collected across all multi-copy reference chromosomes within each reference subgenome:

1. **Global identity clustering**: Alignment identities from all multi-copy reference chromosomes are pooled and fitted with a GMM. The number of components is selected by BIC (Bayesian Information Criterion). Fitting across multiple chromosomes simultaneously provides more signal than fitting each chromosome independently.
2. **Paired validation**: For chromosomes with exactly _k_ copies, each copy should land in a different cluster. The model is accepted if at least 70% of testable chromosomes satisfy this pairing requirement. If BIC initially selects a model that fails paired validation, the tool attempts rescue by testing alternative component counts.
3. **Per-chromosome assignment**: Once global cluster means are established, each contig is assigned to the nearest cluster. Chromosomes where all copies land in the same cluster despite rescue attempts are rank-assigned by identity (highest identity = primary subgenome).

No user configuration is required. The feature runs automatically whenever multiple contigs map to the same reference chromosome.

### Limitations and expectations

**Identity divergence requirement**: The method relies on a detectable difference in alignment identity between subgenomes. Separation is validated by per-chromosome pairing consistency (each chromosome's copies should land in different clusters), not by a fixed identity threshold. Even small identity gaps can be detected if the pattern is consistent across chromosomes. Subgenomes that are nearly equidistant from the reference (e.g., two closely related diploid progenitors) may not produce a consistent pairing signal, and contigs will be labeled `_c1`/`_c2` instead.

**Nucleotide mode with divergent references**: In nucleotide mode (`--synteny-mode nucleotide`, the default) using a cross-species reference at ~60% identity, alignment coverage is lower and identity estimates are noisier than in within-species comparisons. Subgenome segmentation is less reliable in this setting; protein mode (`--synteny-mode protein`) is generally more robust for distantly related references.

**Haplotype-collapsed assemblies**: Assemblers that produce a single primary assembly (e.g., in primary/alternate mode) collapse both haplotypes of each chromosome into one contig. There is nothing to segment — the tool correctly reports a single copy with no subgenome suffix. This is expected behavior.

**Haplotype-phased assemblies**: Assemblers like hifiasm (Hi-C mode) produce separate FASTA files per haplotype (e.g., `h1tg`/`h2tg` contigs). If both haplotypes are combined into a single query FASTA, they typically align to the reference at nearly identical identity levels and may not be separable by the GMM. The tool will label them `_c1`/`_c2` (unsegmented copies) rather than `_B` (distinct subgenome). This is the correct outcome: the two copies represent haplotypes of the same subgenome, not distinct homeologous subgenomes. Note that dnadis does not parse assembler-specific contig name conventions (e.g., `h1tg`/`h2tg`) to infer haplotype groupings — its segmentation is based solely on alignment identity to the reference, so the resulting `_c1`/`_c2` ordering may not correspond to the assembler's `h1`/`h2` labels.

**Recommended workflow for phased assemblies**: Rather than combining haplotypes into a single FASTA, run each haplotype as a separate assembly using `--fofn` or `--assembly-dir`. Each haplotype gets clean chromosome names without `_c` suffixes, and the multi-assembly comparison report shows how the haplotypes differ.

**Autopolyploids**: Autopolyploid genomes (e.g., autotetraploids) carry multiple chromosome sets derived from the same ancestral species. Unlike allopolyploids, the subgenomes have identical evolutionary distance from the reference, so identity-based clustering cannot distinguish them. These will also be labeled `_c1`/`_c2`.

**Minimum data requirement**: Clustering requires at least 4 multi-copy reference chromosomes to attempt model fitting, and at least 3 chromosomes with exactly _k_ copies for paired validation. Assemblies with very few assigned chromosomes may not trigger subgenome segmentation.

**Most effective use case**: Subgenome segmentation is most informative for allopolyploids, where the query carries two or more homeologous chromosome sets derived from different ancestral species, each at a distinct evolutionary distance from the reference. In these cases, the identity gap between subgenomes is driven by real evolutionary divergence and is usually large enough for reliable clustering.

## Installation

### Dependencies
[![CI](https://github.com/eernst/dnadis/actions/workflows/ci.yml/badge.svg)](https://github.com/eernst/dnadis/actions/workflows/ci.yml)

**Required:**
- Python 3.11+
- [intervaltree](https://github.com/chaimleib/intervaltree) - efficient overlap detection
- [miniprot](https://github.com/lh3/miniprot) - protein-to-genome alignment
- [gffread](https://github.com/gpertea/gffread) - GFF3/FASTA processing
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - organelle/rDNA detection

**Optional:**
- [minimap2](https://github.com/lh3/minimap2) or [mm2plus](https://github.com/lh3/mm2plus) - nucleotide synteny QA and read alignment for depth analysis
- [samtools](https://github.com/samtools/samtools) - BAM/CRAM handling for depth analysis
- [mosdepth](https://github.com/brentp/mosdepth) - efficient depth calculation
- [rasusa](https://github.com/mbhall88/rasusa) - FASTQ downsampling for depth analysis
- [centrifuger](https://github.com/mourisl/centrifuger) - contaminant detection
- [taxonkit](https://github.com/shenwei356/taxonkit) + NCBI taxonomy database - taxonomic lineage for contaminant table (see below)
- [RagTag](https://github.com/malonge/RagTag) - improved reference-guided scaffolding (for `--scaffold`; built-in scaffolder used as fallback)
- [executorlib](https://github.com/pyiron/executorlib) + [pysqa](https://github.com/pyiron/pysqa) + [h5py](https://github.com/h5py/h5py) - distributed SLURM job submission (required for `--cluster`; pysqa and h5py are optional executorlib dependencies not pulled in by default)
- [infernal](http://eddylab.org/infernal/) - structure-based rRNA annotation with Rfam covariance models (for rDNA consensus building; enabled by default, skip with `--skip-rdna-consensus`; bundled Rfam database)
- [compleasm](https://github.com/huangnengCSU/compleasm) - BUSCO completeness evaluation (requires `--compleasm-lineage`; install in a **separate conda environment** due to dependency conflicts — see below)
- R with ggplot2, dplyr, readr, stringr, tibble, tidyr, patchwork, ggnewscale, pacman - visualization (enabled by default; skip with `--skip-plot`)
- rmarkdown + pandoc - unified HTML report generation (enabled by default; skip with `--skip-plot`)

### Conda environment

An `environment.yml` is provided for the full installation:

```bash
conda env create -f environment.yml
conda activate dnadis
```

Alternatively, create the environment manually:

**Minimal** — core classification pipeline with interactive HTML reports:

```bash
conda create -n dnadis -c conda-forge -c bioconda \
    python=3.14 intervaltree \
    miniprot gffread blast mm2plus \
    r-base r-ggplot2 r-dplyr r-readr r-stringr r-tibble r-tidyr \
    r-patchwork r-ggnewscale r-pacman r-ggiraph r-htmlwidgets r-scales \
    r-gt r-gtextras r-svglite r-xml2 r-rmarkdown r-ggridges \
    r-colorspace r-ggokabeito r-ggrepel r-showtext r-sysfonts \
    r-systemfonts freetype libxml2 xz pandoc
conda activate dnadis
```

**Full** — all features including reports, contaminant screening, read depth, rDNA annotation, reference-guided scaffolding, SLURM distribution, and taxonomic lineage:

```bash
conda create -n dnadis -c conda-forge -c bioconda \
    python=3.14 intervaltree \
    miniprot gffread blast mm2plus \
    ragtag \
    samtools mosdepth rasusa \
    centrifuger taxonkit infernal \
    executorlib pysqa h5py \
    r-base r-ggplot2 r-dplyr r-readr r-stringr r-tibble r-tidyr \
    r-patchwork r-ggnewscale r-pacman r-ggiraph r-htmlwidgets r-scales \
    r-gt r-gtextras r-svglite r-xml2 r-rmarkdown r-ggridges \
    r-colorspace r-ggokabeito r-ggrepel r-showtext r-sysfonts \
    r-systemfonts freetype libxml2 xz pandoc
conda activate dnadis
```

**Compleasm** (BUSCO completeness evaluation) — must be in a **separate conda environment** due to dependency conflicts (dendropy version clash):

```bash
conda create -n compleasm -c bioconda -c conda-forge compleasm
```

The pipeline auto-detects the `compleasm` conda environment, or you can pass `--compleasm-path /path/to/compleasm` explicitly. Use `--compleasm-lineage <lineage>` (e.g., `embryophyta`, `liliopsida`, `eukaryota`) to enable BUSCO evaluation.

**Taxonkit** — for full taxonomic lineage (Domain → Family → Genus → Species) in the contaminant table, download the NCBI taxonomy database after installing taxonkit:

```bash
taxonkit download --data-dir ~/.taxonkit
```

Without taxonkit or the database, contaminant tables show species names only.

### Development

To run the test suite:
```bash
conda install -n dnadis -c conda-forge pytest pytest-cov
conda run -n dnadis pytest -q
```

Latest tested conda package versions (CI):
<!-- conda-versions-start -->
- python: 3.14.4
- miniprot: 0.18
- gffread: 0.12.9
- blast: 2.17.0
- mm2plus: 1.2
- ragtag: 2.1.0
- centrifuger: 1.1.1
<!-- conda-versions-end -->

## Quick Start

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o output/ \
    --ref-gff3 reference.gff3
```

## Usage

### Required arguments

| Argument | Description |
|----------|-------------|
| `-r, --ref` | Reference genome FASTA (can be gzipped) |
| `-q, --query` | Query assembly FASTA to classify (single-assembly mode) |
| `-o, --output-dir` | Output directory (reference/ and per-assembly subdirectories created inside) |
| `--ref-gff3` | Reference GFF3 with protein-coding gene annotations. Required for `--synteny-mode protein`. |

### Multi-assembly mode

`-q/--query` is for single-assembly mode. To process multiple assemblies against the same reference, use one of:

| Argument | Description |
|----------|-------------|
| `--fofn` | Tab-separated file-of-filenames with column `path` and optional columns `name` and `reads`. One assembly per row. |
| `--assembly-dir` | Directory to scan for FASTA files (`.fasta`, `.fa`, `.fna`, `.fasta.gz`, `.fa.gz`, `.fna.gz`). Assembly names are derived from filenames. |
| `--assembly-name` | Override the assembly name for the single-assembly (`-q`) case (default: derived from query filename stem). |

`--fofn` and `--assembly-dir` are mutually exclusive and cannot be combined with `-q`. When ≥2 assemblies complete, cross-assembly comparison outputs (summary TSV, chromosome completeness TSV, and an interactive comparison HTML report) are written to the top-level output directory (see [Multi-assembly outputs](#multi-assembly-outputs)).

### Common options

| Argument | Description | Default |
|----------|-------------|---------|
| `-t, --threads` | Number of threads | 8 |
| `--skip-plot` | Skip unified HTML report generation | off |
| `--no-self-contained-html` | Produce non-self-contained HTML reports (HTML + companion `_files/` directory). Faster rendering and smaller files. Default: self-contained (all resources embedded in a single portable file). | off |
| `-v, --verbose` | Enable verbose (DEBUG level) logging | off |
| `--quiet` | Suppress INFO messages (only warnings and errors) | off |
| `--log-file` | Write logs to file (in addition to stderr) | none |
| `--config` | Load configuration from TOML file (CLI args override) | none |
| `--dump-config` | Print TOML config template and exit | off |
| `-C, --chr-like-minlen` | Min contig length for chromosome classification | 25% of smallest nuclear ref chromosome |
| `--add-subgenome-suffix` | Suffix for non-polyploid references (e.g., 'A') | none |

### Classification references

| Argument | Description |
|----------|-------------|
| `--chrC-ref` | Chloroplast reference FASTA (default: extract from --ref) |
| `--chrM-ref` | Mitochondrion reference FASTA (default: extract from --ref) |
| `--rdna-ref` | rDNA reference FASTA, or 'default' for bundled Arabidopsis 45S |
| `--skip-rdna-consensus` | Skip building consensus 45S rDNA and annotating rDNA loci (enabled by default) |
| `--centrifuger-idx` | Centrifuger index prefix for contaminant screening |

### Read depth analysis

| Argument | Description | Default |
|----------|-------------|---------|
| `--reads` | Reads for depth analysis (FASTQ/BAM/CRAM). Auto-detects format. | none |
| `--reads-type` | Read type: `lrhq`, `r9`, or `sr` | lrhq |
| `--skip-depth` | Skip depth analysis even if `--reads` provided | off |
| `--depth-window-size` | Window size for mosdepth | 1000 |
| `--depth-target-coverage` | Target coverage for downsampling before alignment (0 to disable) | 0 |
| `--keep-depth-bam` | Keep aligned BAM file after depth analysis | off |

The `--reads` option accepts:
- FASTQ files (`.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`)
- Unaligned BAM/CRAM files (will be aligned with minimap2)
- Pre-aligned BAM/CRAM files (used directly)

Read type to minimap2 preset mapping:
- `lrhq` → `-x lr:hqae` (PacBio HiFi, ONT R10; high-quality long reads with error rate < 1%)
- `r9` → `-x map-ont` (ONT R9 reads, standard accuracy)
- `sr` → `-x sr` (Illumina short reads)

**Read downsampling**: If `--depth-target-coverage` is set to a non-zero value (e.g., 20 for 20X coverage), reads are automatically downsampled before alignment using rasusa (for FASTQ) or samtools (for BAM/CRAM). This reduces computational time while maintaining sufficient coverage for depth-based quality assessment. Pre-aligned BAM/CRAM files are not downsampled. Default: 0 (disabled).

**Depth caching**: Alignment and depth results are cached and automatically reused on subsequent runs if inputs match. Use `--keep-depth-bam` to retain the aligned BAM file after depth calculation (default: deleted to save space). Cached depth results can be reused even if the BAM is deleted.

### Synteny mode selection

| Argument | Description | Default |
|----------|-------------|---------|
| `--synteny-mode` | Synteny evidence source: `protein` (miniprot) or `nucleotide` (minimap2) | nucleotide |

**Nucleotide mode** (default): Uses minimap2 whole-genome nucleotide alignment with permissive chaining for chromosome-scale structural composition analysis. Creates megabase-scale synteny blocks by chaining through repetitive regions. Suitable for both within-species and cross-species comparisons. Ideal for identifying structural features like chromosomal fusions, homeologous recombination, or introgression events.

**Protein mode**: Uses miniprot protein-anchored synteny. Requires `--ref-gff3`. Ideal for detecting conserved gene content across distantly related species.

### Pipeline toggles

| Argument | Description |
|----------|-------------|
| `--skip-organelles` | Skip organelle detection |
| `--skip-rdna` | Skip rDNA detection |
| `--skip-contaminants` | Skip contaminant detection |

### BUSCO completeness options

| Argument | Description | Default |
|----------|-------------|---------|
| `--compleasm-lineage` | BUSCO lineage for compleasm evaluation (e.g., `eukaryota`, `viridiplantae`, `embryophyta`). Required for compleasm to run. | none |
| `--compleasm-library` | Path to pre-downloaded compleasm lineage files (avoids runtime download) | auto-download |
| `--compleasm-path` | Path to compleasm executable. If unset, auto-detects from a `compleasm` conda environment or `PATH`. | auto-detect |
| `--skip-compleasm` | Skip compleasm even if `--compleasm-lineage` is specified | off |

When `--compleasm-lineage` is set, phase 18 runs compleasm on two FASTA subsets: chromosome-assigned contigs (`*.chrs.fasta`) and non-chromosome contigs (`*.non_chrs.fasta`, combining debris + unclassified + contaminants). Both runs are submitted in parallel. Results are included in the per-assembly unified HTML report and in the multi-assembly `comparison_summary.tsv`.

### Scaffolding options

| Argument | Description | Default |
|----------|-------------|---------|
| `--scaffold` | Produce reference-guided scaffolded chromosome sequences (uses RagTag if available, otherwise built-in scaffolder) | off |
| `--scaffold-gap-size` | Number of Ns between contigs in scaffolded output | 100 |

When `--scaffold` is enabled, chromosome-assigned contigs are grouped by reference chromosome and ordered into pseudomolecules. The scaffolder handles haplotype-aware grouping for polyploid assemblies (e.g., contigs assigned to chr1A are scaffolded separately from chr1B). Single T2T contigs that span a full chromosome produce trivial (single-component) AGP entries. Multi-contig chromosomes are ordered by reference position, either via RagTag (if installed) or the built-in scaffolder.

### Distributed computing (SLURM cluster)

| Argument | Description | Default |
|----------|-------------|---------|
| `--cluster` | Submit compute phases as SLURM jobs via executorlib | off |
| `--partition` | SLURM partition for distributed jobs | cpuq |
| `--max-threads-dist` | Max threads per distributed job | 64 |
| `--max-mem-dist` | Max memory (GB) per distributed job | 128 |
| `--max-time-dist` | Max wall time (minutes) per distributed job | 720 |

**Requires [executorlib](https://github.com/pyiron/executorlib), [pysqa](https://github.com/pyiron/pysqa), and [h5py](https://github.com/h5py/h5py):**
```bash
conda install -n dnadis -c conda-forge executorlib pysqa h5py
```

When `--cluster` is enabled, compute-intensive phases (synteny alignment, BLAST detection, debris detection, contaminant screening, read depth, compleasm) are submitted as individual SLURM jobs with per-job resource control. In multi-assembly mode, assemblies run concurrently with each submitting its own SLURM jobs. Without `--cluster`, all phases run locally.

## Output Files

### FASTA outputs

| File | Description |
|------|-------------|
| `*.chrs.fasta` | Chromosome-assigned contigs (renamed and reoriented; scaffolded pseudomolecules if `--scaffold`) |
| `*.organelles.fasta` | Organelle contigs (chrC, chrM) |
| `*.rdna.fasta` | rDNA-containing contigs |
| `*.contaminants.fasta` | Contaminant contigs |
| `*.debris.fasta` | Assembly debris (fragments, duplicates) |
| `*.unclassified.fasta` | Contigs that couldn't be classified |
| `*.non_chrs.fasta` | Combined non-chromosome contigs (debris + unclassified + contaminants; produced when `--compleasm-lineage` is set) |
| `*.scaffolded.fasta` | Scaffolded chromosome pseudomolecules (if `--scaffold`) |
| `*.scaffolded.agp` | AGP 2.0 file describing scaffold structure (if `--scaffold`) |

### Summary tables

| File | Description |
|------|-------------|
| `*.contig_summary.tsv` | Per-contig classification with metrics |
| `*.evidence_summary.tsv` | Synteny evidence per contig × reference |
| `*.segments.tsv` | Individual synteny segments |
| `*.macro_blocks.tsv` | Aggregated synteny macro-blocks |
| `*.ref_lengths.tsv` | Reference chromosome lengths |
| `*.contaminants.tsv` | Detailed contaminant summary with taxonomic lineage (if contaminants detected) |
| `*.rdna_annotations.tsv` | rRNA sub-feature annotations in TSV format (produced by default; skip with `--skip-rdna-consensus`) |
| `*.rdna_arrays.tsv` | rDNA array locations per contig (produced when arrays are detected; skip with `--skip-rdna-consensus`) |

### GFF3 outputs

| File | Description |
|------|-------------|
| `*.rdna_annotations.gff3` | Hierarchical rRNA gene annotations with 18S, 5.8S, 25S, ITS1, ITS2 sub-features (produced by default; skip with `--skip-rdna-consensus`) |

### Multi-assembly outputs

Produced in the top-level output directory when ≥2 assemblies complete. File names are prefixed with `--comparison-name` (default: `comparison`):

| File | Description |
|------|-------------|
| `{prefix}_summary.tsv` | Per-assembly classification counts and BUSCO completeness scores (S/D/F/I/M) for chromosome and non-chromosome contig sets |
| `{prefix}_chromosome_completeness.tsv` | Per-reference-chromosome coverage and completeness across all assemblies |
| `{prefix}.comparison_report.html` | Interactive HTML comparison report with gt tables summarizing classification and BUSCO results across assemblies |

### Key columns in `contig_summary.tsv`

| Column | Description |
|--------|-------------|
| `contig` | New contig name (with chromosome assignment) |
| `original_name` | Original contig name from input FASTA |
| `classification` | Category: chrom_assigned, chrom_unassigned, organelle_complete, organelle_debris, rDNA, contaminant, chrom_debris, debris, unclassified |
| `classification_confidence` | Confidence level: high, medium, or low |
| `reversed` | Whether contig was reverse-complemented |
| `contaminant_taxid` | NCBI taxonomy ID (for contaminants) |
| `contaminant_sci` | Scientific name (for contaminants) |
| `assigned_ref_id` | Best-matching reference chromosome |
| `ref_gene_proportion` | Fraction of reference chromosome genes aligned to this contig (0.0-1.0) |
| `genes_per_Mbp` | Gene density (aligned reference genes per Mbp of query sequence) |
| `depth_mean` | Mean read depth (if `--reads` provided) |
| `depth_median` | Median read depth (if `--reads` provided) |
| `depth_std` | Standard deviation of read depth (if `--reads` provided) |
| `depth_breadth_1x` | Fraction of bases with ≥1x coverage (if `--reads` provided) |
| `depth_breadth_10x` | Fraction of bases with ≥10x coverage (if `--reads` provided) |

For complete column documentation for all TSV files, see [docs/output_formats.md](docs/output_formats.md).

### Visualization

Enabled by default; skip with `--skip-plot` (requires R with ggplot2 and related packages, plus rmarkdown + pandoc).

| File | Description |
|------|-------------|
| `*.assembly_report.html` | HTML report per assembly with all plots (chromosome overview, classification bar, read depth overview, contaminant table). Self-contained by default (all resources embedded). Use `--no-self-contained-html` to produce HTML + companion `*_files/` directory for faster rendering. |
| `*.chromosome_overview.pdf` | Multi-panel plot showing contig composition, subgenome support, and alignment identity (exported from the assembly report) |
| `*.depth_overview.pdf` | Read depth visualization by classification and chromosome (if `--reads` provided; exported from the assembly report) |

## Classification Pipeline

Reference preparation runs first: read reference genome, compute GC statistics, prepare organelle/rDNA references, and (in protein mode) extract proteins from GFF3. Then the per-assembly pipeline runs these phases:

| Phase | Description |
|-------|-------------|
| 1 | **Read query assembly** — parse FASTA, compute contig lengths and GC |
| 2 | **Synteny analysis** — protein mode (miniprot) or nucleotide mode (minimap2); chain alignments into synteny blocks and identify chromosome candidates |
| 3 | **Organelle detection** — BLAST against chrC/chrM references (skip with `--skip-organelles`) |
| 4 | **rDNA detection** — BLAST against rDNA reference (skip with `--skip-rdna`) |
| 5 | **Chromosome debris detection** — high-coverage, high-identity matches to assembled chromosomes (minimap2) |
| 6 | **Contaminant detection** — centrifuger taxonomic classification with two-gate filtering (score ≥1000, coverage ≥0.50; requires `--centrifuger-idx`) |
| 7 | **Debris classification** — reference-based debris detection for remaining contigs |
| 8 | **Gene count statistics** — compute gene proportion metrics (if GFF3 provided) |
| 9 | **Orientation determination** — determine strand for chromosome contigs based on synteny votes |
| 10 | **Telomere detection** — scan contig ends for telomeric repeats (skip with `--skip-telomeres`) |
| 11 | **Classification** — assign all contigs to categories with confidence levels; rename contigs to reference-based names (e.g., `chr1A`, `chr1A_f1`, `contig_1`) |
| 12 | **Rearrangement detection** — identify contigs with evidence from multiple reference chromosomes |
| 13 | **Read depth analysis** — align reads and compute per-contig depth metrics (optional, requires `--reads`) |
| 14 | **rDNA consensus building** — build species-specific 45S rDNA consensus and annotate rRNA sub-features (skip with `--skip-rdna-consensus`) |
| 15 | **Reference-guided scaffolding** — order and orient contigs into chromosome-scale pseudomolecules with AGP output (optional, `--scaffold`) |
| 16 | **Write FASTA outputs** — classified FASTA files (chromosomes, organelles, rDNA, etc.) |
| 17 | **Write summary TSV** — per-contig classification table, evidence summaries, and visualizations |
| 18 | **BUSCO completeness evaluation** — compleasm run on `*.chrs.fasta` and `*.non_chrs.fasta` in parallel (optional, requires `--compleasm-lineage`; skip with `--skip-compleasm`) |

## rDNA Consensus and Annotation

By default, the tool builds a species-specific consensus 45S rDNA sequence from the query assembly and uses it to annotate ribosomal RNA genes with accurate sub-feature boundaries. This step runs automatically unless `--skip-rdna-consensus` is set.

**Pipeline:**
1. Extract rDNA-containing regions from the assembly using BLAST against a reference 45S sequence
2. Self-align extracted regions to detect repeat periodicity
3. Cluster individual copies and select a consensus representative
4. Annotate rRNA sub-features (18S, 5.8S, 25S/28S) and internal transcribed spacers (ITS1, ITS2) using Infernal/cmscan with Rfam covariance models
5. Write GFF3 file with hierarchical feature structure (rRNA_gene parent with rRNA and ITS children)

**Sub-feature annotation:**
Uses Infernal covariance models from Rfam 15.0 for structure-based rRNA boundary detection. Provides accurate gene boundaries based on conserved secondary structure. Bundled models (5S, 5.8S, 18S, 28S) are stored in `dnadis/data/rfam/euk-rrna.cm` and automatically pressed (indexed) on first use.

**Output files:**
- `*.rdna_annotations.gff3`: Hierarchical GFF3 with proper Sequence Ontology terms (SO:0001637 for rRNA_gene, SO:0000252 for rRNA, SO:0000635 for ITS)
- Each rRNA locus includes properly nested child features for 18S, ITS1, 5.8S, ITS2, and 25S (or 28S)

**Dependencies:**
- Required: BLAST+ (makeblastdb, blastn) for locus detection
- Required: [Infernal](http://eddylab.org/infernal/) for sub-feature annotation (included in the recommended conda environment)

**Use cases:**
- Locate nucleolar organizer regions (NORs) on chromosomes
- Quantify rDNA copy number across the genome
- Generate species-specific probe for improved rDNA detection
- Annotate rRNA gene boundaries for downstream analysis

## Contig Classification Categories

| Category | Description |
|----------|-------------|
| `chrom_assigned` | Chromosome-length contig assigned to a reference chromosome via synteny |
| `chrom_unassigned` | Chromosome-length contig without reference assignment (novel or failed synteny gates) |
| `organelle_complete` | Complete organelle genome (chrC or chrM) |
| `organelle_debris` | Partial organelle sequence |
| `rDNA` | Ribosomal DNA repeat unit |
| `contaminant` | Sequence from contaminating organism |
| `chrom_debris` | High-coverage (≥80%), high-identity (≥90%) duplicate of an assembled chromosome contig |
| `debris` | Assembly debris with reference nucleotide coverage (≥50%) or protein homology (≥2 miniprot hits) |
| `unclassified` | Could not be classified |

## Classification Confidence Levels

Each contig is assigned a confidence level (`high`, `medium`, or `low`) indicating the reliability of its classification based on multiple lines of evidence.

**Mode-dependent criteria**: Confidence scoring for `chrom_assigned` and `debris` categories differs between protein and nucleotide synteny modes to account for the different types of evidence available.

### Confidence Criteria by Category

| Classification | High | Medium | Low |
|----------------|------|--------|-----|
| **chrom_assigned** (protein mode) | Gene proportion ≥20% AND GC deviation <2σ | Gene proportion 10-20% OR GC deviation 2-3σ | Gene proportion <10% OR GC deviation >3σ |
| **chrom_assigned** (nucleotide mode) | Ref coverage ≥30% AND identity ≥50% AND GC deviation <2σ | Ref coverage 10-30% OR GC deviation 2-3σ | Ref coverage <10% OR identity <50% OR GC deviation >3σ |
| **chrom_unassigned** | — | GC deviation <2σ | GC deviation ≥2σ |
| **organelle_complete** | Coverage ≥90% | Coverage 80-90% | — |
| **organelle_debris** | — | Coverage ≥60% | Coverage <60% |
| **rDNA** | Coverage ≥80% AND identity ≥95% | Coverage 50-80% OR identity <95% | Coverage <60% |
| **contaminant** | Coverage ≥80% OR GC deviation >2σ | Coverage 50-80% | Coverage <50% |
| **chrom_debris** | Coverage ≥90%, identity ≥95%, GC deviation <3σ | Coverage 80-90% OR identity 90-95% OR GC deviation ≥3σ | — |
| **debris** (protein mode) | Coverage ≥80% OR ≥5 protein hits | Coverage 50-80% OR 2-5 protein hits | GC deviation >3σ AND no synteny |
| **debris** (nucleotide mode) | Coverage ≥80% | Coverage 50-80% | Coverage <55% OR GC deviation >3σ |
| **unclassified** | — | — | Always (no evidence) |

### Evidence Factors

**Gene proportion** (`ref_gene_proportion`) [protein mode only]: Fraction of reference chromosome genes aligned to the query contig (e.g., 0.85 = 85% of reference genes found on this contig). Higher values indicate stronger synteny support. This metric assesses completeness relative to the reference chromosome.

**Reference coverage** (`ref_coverage`) [both modes]: Fraction of reference chromosome spanned by aligned segments (0.0-1.0). In nucleotide mode, this is the primary metric for assessing chromosome assignment quality.

**Synteny score** (`synteny_score`) [mode-dependent]: Summary score for assignment strength, computed differently by mode:
- Protein mode: `min(1.0, gene_proportion × 2)` — normalized gene proportion
- Nucleotide mode: `min(1.0, ref_coverage)` — directly uses reference coverage

**GC deviation**: How many standard deviations (σ) the contig's GC content differs from a GC baseline. Large deviations may indicate contamination or unusual sequences.
- For `chrom_assigned`: compared to reference nuclear genome GC (validates synteny-based assignment)
- For all other categories: compared to assembly chromosome GC (more appropriate for divergent genomes where the assembly may differ significantly from the reference)

**Coverage**: Fraction of contig covered by alignments (0.0-1.0). Higher coverage indicates more complete matches.

**Identity**: Alignment identity (0.0-1.0). Values ≥0.95 indicate very high sequence similarity. In nucleotide mode, identity <0.5 triggers a confidence downgrade for `chrom_assigned` contigs.

**Protein hits** [protein mode only]: Number of unique reference protein-coding genes with miniprot alignments to the contig.

## Thresholds

### Chromosome assignment

**Common thresholds (both modes):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assign-min-frac` | 0.10 | Min synteny coverage of contig |
| `--assign-min-ratio` | 1.25 | Min best/second score ratio |
| `--min-span-frac` | 0.20 | Min span fraction of contig |
| `--min-span-bp` | 50000 | Min absolute span in bp |

**Protein mode specific:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assign-minlen-protein` | 150 | Min target span (bp) for protein-anchor synteny block building |
| `--miniprot-min-genes` | 3 | Min unique genes for assignment |
| `--miniprot-min-segments` | 5 | Min synteny segments for protein mode |

**Nucleotide mode specific:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assign-minlen` | 10000 | Min alignment span for nucleotide synteny block building |
| Min segments | 1 | Hardcoded (nucleotide mode requires ≥1 segment) |

All gate criteria must be satisfied for chromosome assignment (AND logic). The segment count threshold differs: protein mode requires ≥5 segments (configurable), while nucleotide mode requires ≥1 segment (hardcoded) because perfect full-length alignments produce fewer segments.

### Organelle detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--organelle-min-cov` | 0.80 | Min query coverage |
| `--chrC-len-tolerance` | 0.05 | Length tolerance for chrC |
| `--chrM-len-tolerance` | 0.20 | Length tolerance for chrM |

### Contaminant detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--contaminant-min-score` | 1000 | Min centrifuger score (~1kb matching sequence with k=31) |
| `--contaminant-min-coverage` | 0.50 | Min query coverage (low coverage may indicate conserved genes, not contamination) |

### Debris detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--chr-debris-min-cov` | 0.80 | Min coverage vs chromosomes |
| `--chr-debris-min-identity` | 0.90 | Min identity vs chromosomes |
| `--debris-min-cov` | 0.50 | Min coverage vs reference |

## Configuration Files

`dnadis` supports TOML configuration files for managing complex parameter sets and reproducible analysis workflows.

### Generating a config template

```bash
./dnadis.py --dump-config > my_config.toml
```

This creates a complete configuration file with all parameters and their current default values. Edit this file to customize your analysis.

### Using a config file

```bash
./dnadis.py --config my_config.toml
```

**Important:** Command-line arguments override config file values. This allows you to have a base configuration but adjust specific parameters via CLI:

```bash
# Use config but override threads
./dnadis.py --config my_config.toml -t 64
```

### Example config file

```toml
# Required
[required]
ref = "/path/to/reference.fasta"
query = "/path/to/assembly.fasta"
output_dir = "/path/to/output_directory"
ref_gff3 = "/path/to/reference.gff3"

# Common Options
[common]
threads = 32
plot = true
chr_like_minlen = 1000000

# Read Depth Analysis
[read_depth]
reads = "/path/to/reads.fastq.gz"
reads_type = "lrhq"
depth_target_coverage = 20.0
keep_depth_bam = false

# Reference Inputs For Classification
[reference_inputs]
centrifuger_idx = "/path/to/centrifuger/nt"
rdna_ref = "default"

# Thresholds can be customized per section
[thresholds_chromosome]
assign_min_frac = 0.10
assign_min_ratio = 1.25
```

## Example Workflows

### Basic chromosome assignment (protein mode)

```bash
./dnadis.py \
    -r TAIR10.fasta \
    -q my_assembly.fasta \
    -o results/ \
    --synteny-mode protein \
    --ref-gff3 TAIR10.gff3 \
    -t 32

```

### Nucleotide mode for structural composition analysis

Use nucleotide mode when you want to detect chromosome-scale structural features like fusions, translocations, or homeologous exchanges. Nucleotide mode is the default and uses whole-genome nucleotide alignment; it is suitable for both within-species and cross-species comparisons.

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o results/ \
    -t 32

```

Note: Nucleotide mode does not require `--ref-gff3`, but if provided, gene count statistics will be included in the output.

### Polyploid genome with contaminant screening

```bash
./dnadis.py \
    -r wheat_ref.fasta \
    -q wheat_assembly.fasta \
    -o results/ \
    --ref-gff3 wheat_ref.gff3 \
    --centrifuger-idx /path/to/centrifuger/nt \
    --rdna-ref default \
    -t 64

```

### Non-polyploid with subgenome suffix

```bash
./dnadis.py \
    -r rice_ref.fasta \
    -q rice_assembly.fasta \
    -o results/ \
    --ref-gff3 rice_ref.gff3 \
    --add-subgenome-suffix A

```

### With reference-guided scaffolding

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o results/ \
    --ref-gff3 reference.gff3 \
    --scaffold \
    -t 32
```

This produces `*.scaffolded.fasta` (chromosome pseudomolecules) and `*.scaffolded.agp` (AGP 2.0 describing scaffold structure). If RagTag is installed, it is used for ordering; otherwise, the built-in scaffolder orders contigs by reference position.

### With read depth analysis (HiFi reads)

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o results/ \
    --ref-gff3 reference.gff3 \
    --reads hifi_reads.fastq.gz \
    --reads-type lrhq \
    -t 32
```

### With pre-aligned BAM for depth analysis

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o results/ \
    --ref-gff3 reference.gff3 \
    --reads aligned_reads.bam \
    -t 32
```

### Using configuration file with verbose logging

```bash
# Generate config template
./dnadis.py --dump-config > wheat_config.toml

# Edit wheat_config.toml to set paths and parameters

# Run with config and verbose logging
./dnadis.py \
    --config wheat_config.toml \
    --verbose \
    --log-file wheat_analysis.log
```

### With downsampled reads for faster depth analysis

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o results/ \
    --ref-gff3 reference.gff3 \
    --reads hifi_reads.fastq.gz \
    --reads-type lrhq \
    --depth-target-coverage 20 \
    --keep-depth-bam \
    -t 32
```

### With BUSCO completeness evaluation

```bash
./dnadis.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o results/ \
    --compleasm-lineage embryophyta \
    -t 32
```

Compleasm runs on chromosome-assigned contigs and non-chromosome contigs in parallel. Results appear in `*.assembly_report.html` and (in multi-assembly runs) in `comparison_summary.tsv`. The pipeline auto-detects compleasm from a `compleasm` conda environment, or you can pass `--compleasm-path` explicitly.

### Multi-assembly mode (file-of-filenames)

Analyze multiple assemblies against a shared reference. The TSV file has columns `path`, `name`, and optionally `reads` (one assembly per row):

```bash
# assemblies.tsv: path <tab> name [<tab> reads]
./dnadis.py \
    -r reference.fasta \
    --ref-gff3 reference.gff3 \
    --fofn assemblies.tsv \
    -o multi_output/ \
    --compleasm-lineage viridiplantae \
    -t 32
```

Cross-assembly outputs (`comparison_summary.tsv`, `comparison_chromosome_completeness.tsv`, `comparison.comparison_report.html`) are written to `multi_output/` once all assemblies complete.

### Multi-assembly mode with SLURM cluster

```bash
./dnadis.py \
    -r reference.fasta \
    --ref-gff3 reference.gff3 \
    --fofn assemblies.tsv \
    -o multi_output/ \
    --cluster \
    --partition cpuq \
    --max-threads-dist 32 \
    --max-mem-dist 128
```

Assemblies run concurrently; each submits its own SLURM jobs for compute-intensive phases.

## Algorithm Details

### Synteny Mode Details

See [Synteny mode selection](#synteny-mode-selection) for an overview of the two modes.

#### Nucleotide mode chaining

Minimap2 alignments use permissive chaining parameters to create continuous chromosome-scale blocks:
- `--max-chain-skip=300`: Chain through repetitive regions
- `-z 10000,1000`: Tolerate long gaps in alignment chains
- `-r 50000`: Large bandwidth (50kb) to accommodate structural variations

This permissive approach chains through homopolymer runs, tandem repeats, and ambiguous regions that would otherwise fragment alignments. The resulting megabase-scale blocks are suitable for chromosome classification and compositional analysis, balanced by downstream filtering (identity thresholds, minimum alignment length, gate filtering) to prevent spurious assignments.

#### Protein mode pipeline

Miniprot alignments are:
1. Filtered by alignment quality and identity
2. Mapped to reference genomic coordinates via GFF3
3. Overlapping hits filtered by identity (keeps highest-identity hit at each position)
4. Chained into collinear synteny blocks
5. Aggregated per contig × reference chromosome
6. Scored and ranked for assignment

### Gate-based assignment

Contigs must satisfy ALL gate criteria to be assigned to a chromosome (AND logic). This prevents spurious assignments from single conserved genes, repetitive sequences, or low-complexity regions. See [Chromosome assignment thresholds](#chromosome-assignment) for the full parameter list and mode-specific defaults.

### Reference assignment scoring

In nucleotide mode, each contig is assigned to the reference chromosome with the highest reference span fraction: `ref_span_bp / ref_length`. This metric answers the biologically meaningful question — what fraction of the reference chromosome does this contig represent? — and is size-normalized, preventing the raw-score bias toward larger chromosomes that would otherwise skew assignments in translocation cases. The metric is computed directly in chain parsing from reference lengths extracted from the PAF file.

After initial assignment, a dedicated pass (`_resolve_reciprocal_translocations()`) handles a specific edge case: when exactly two contigs are both assigned to the same reference chromosome and both share the same second-best reference (which has no contigs assigned), this is the signature of a reciprocal translocation. The weaker contig is reassigned to the partner reference.

In protein mode, span fraction falls back to raw synteny score (sum of chain scores; see `--assign-ref-score`) because miniprot PAF does not carry reference chromosome lengths.

Each contig is scored independently; this is not a conflict-aware or globally optimal assignment. Multiple contigs can be assigned to the same reference chromosome, which is expected and correct for polyploid assemblies.

### Chain scoring modes

The `--assign-chain-score` parameter controls how synteny chains are weighted for chromosome assignment:

| Mode | Formula | Description |
|------|---------|-------------|
| `matches` (default) | `matches` | Most permissive. Favors chains with high absolute matching base count, regardless of alignment quality. Best for initial assignments. |
| `qbp_ident` | `qbp × identity` | Balances coverage and identity. Favors chains that cover large query regions with good identity. Useful for discriminating between homeologs. |
| `matches_ident` | `matches × identity` | Most stringent. Heavily weights high-identity alignments. Use when reference and query are closely related and you want to penalize divergent matches. |

Where:
- `matches`: Total matching bases across chain
- `qbp`: Query base pairs (union of alignment intervals on query contig)
- `identity`: matches / aln_len

The chain score determines ranking for `--assign-chain-topk` (default: top 3 chains per contig-reference pair contribute to assignment score).

### Subgenome discrimination

For polyploid genomes (e.g., wheat with A, B, D subgenomes), dnadis discriminates between chromosome sets using alignment identity distributions.

**How it works:**
1. Reference chromosomes are labeled with subgenome suffixes (chr1A, chr1B, chr1D)
2. Query contigs are assigned to reference chromosomes based on synteny evidence
3. Identity distributions distinguish homeologous chromosome sets
4. Visualization shows identity grouped by assigned subgenome with mean ± SD error bars

**Expected pattern:**
- Within-subgenome assignments show higher identity (e.g., query A-subgenome → reference chr*A: ~98-99% identity)
- Cross-subgenome assignments show lower identity (e.g., query A-subgenome → reference chr*B: ~94-96% identity)

**Identity calculation:** Identity is computed per (contig, reference) pair as the aggregate across all chains: `sum(matches) / sum(aln_len)`. For polyploid genomes, this metric effectively captures the evolutionary distance between homeologous chromosome sets.

**For non-polyploid genomes:** Use `--add-subgenome-suffix A` to add a single subgenome label to the reference. This can be used to "bootstrap" a synthetic hybrid reference chromosome set.

### Debris detection algorithm

Two complementary approaches:

1. **Chromosome debris** - Contigs with high-coverage (≥80%) and high-identity (≥90%) matches to already-assembled chromosome contigs (detected via minimap2 asm5 alignment)

2. **Reference-based debris** - Contigs with moderate reference coverage (>50%) but insufficient synteny support for chromosome assignment

### Contaminant visualization

When contaminants are detected (with `--centrifuger-idx`) and plotting is enabled (the default; disable with `--skip-plot`), dnadis includes a **contaminant table** in the assembly report (`*.assembly_report.html`): an interactive HTML table showing top contaminants ranked by abundance (depth × length if depth data available, otherwise total length). Features include species-level aggregation with binomial names, inline gradient bars for Total Mb and Depth columns, colored domain badges, family grouping, and min-max spread values for multi-contig entries. HTML-only output (CSS gradients don't render to PDF).

The table filters to high-confidence contaminants (coverage ≥ `--contaminant-min-coverage`, default 0.50) to reduce noise from conserved gene matches. See [docs/output_formats.md](docs/output_formats.md#contaminant-table-visualization-html) for detailed documentation on the table format and interpretation.

## Glossary

| Term | Definition |
|------|------------|
| **Alignment identity** | Fraction of matching bases in alignment (`matches / aln_len`). Differs from BLAST "percent identity" which includes query length normalization. |
| **Synteny block/chain** | Group of collinear alignments chained together based on query position, reference position, and diagonal offset. Represents a contiguous region of conserved gene order. |
| **Union bp** | Total non-overlapping base pairs covered by alignments, computed by merging overlapping intervals. Measures breadth of alignment coverage. |
| **Gate-based filtering** | Chromosome assignment requires passing ALL criteria (AND logic): minimum segments, genes, span bp, and span fraction. Prevents spurious assignments from single conserved genes or repetitive elements. |
| **Homeologs** | Corresponding chromosomes from different subgenomes in a polyploid (e.g., chr1A vs chr1B in wheat). Show sequence similarity (85-95% identity) but are distinct from homologs (between species) or paralogs (within genome). |
| **Subgenome** | One of multiple ancestral genome copies in a polyploid. Identified by single-letter suffix (A, B, C, etc.) on chromosome names. |

## Citation

If you use dnadis in your work, please cite it. The repository includes a
[`CITATION.cff`](CITATION.cff) file; GitHub displays a "Cite this repository"
button that exports BibTeX or APA formats.

A suggested citation:

> Ernst, E. (2026). *dnadis: genome assembly comprehension and curation tool for contig classification, homology assignment, and quality control* (Version 0.1.0) [Computer software]. https://github.com/eernst/dnadis

Please also cite the underlying tools that dnadis invokes:
- [miniprot](https://github.com/lh3/miniprot)
- [minimap2](https://github.com/lh3/minimap2)
- [gffread](https://github.com/gpertea/gffread)
- [centrifuger](https://github.com/mourisl/centrifuger) (if used)
- [taxonkit](https://github.com/shenwei356/taxonkit) (if used for contaminant visualization)
- [Infernal](http://eddylab.org/infernal/) and [Rfam](https://rfam.org/) (if rDNA consensus building was used for rRNA annotation):
  - Nawrocki EP, Eddy SR (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29(22):2933-2935. <https://doi.org/10.1093/bioinformatics/btt509>
  - Kalvari I, et al. (2021). Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Research, 49(D1):D192-D200. <https://doi.org/10.1093/nar/gkaa1047>

## License

BSD 3-Clause License. Copyright 2026 HHMI. See [LICENSE](LICENSE) for the full text.

## Author

Evan Ernst (eernst@cshl.edu)
