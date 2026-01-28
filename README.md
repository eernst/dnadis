# final_finalizer

**Genome assembly finalization tool for contig classification and quality control**

`final_finalizer` classifies contigs from a *de novo* genome assembly into biological categories (chromosomes, organelles, rDNA, contaminants, debris, unclassified) using protein-anchored synteny evidence, organelle/rDNA alignments, and taxonomic classification.

## Features

- **Dual synteny modes**:
  - **Protein mode** (default): Protein-anchored synteny blocks (miniprot) for gene-level classification
  - **Nucleotide mode**: Whole-genome nucleotide alignment (minimap2) for structural composition analysis
- **Chromosome assignment** with gate-based filtering to prevent spurious assignments
- **Subgenome resolution** for polyploid genomes with multiple chromosome sets
- **Organelle detection** (chloroplast/chrC, mitochondrion/chrM) via BLAST
- **rDNA contig identification** via BLAST against reference rDNA
- **Contaminant screening** via centrifuger taxonomic classification
- **Debris detection** for assembly fragments (chromosome debris, organelle debris)
- **Chimera flagging** for contigs with evidence from multiple chromosomes
- **Orientation determination** for chromosome-assigned contigs
- **Read depth analysis** (optional) with automated downsampling and caching
- **Publication-ready visualizations** (PDF plots via R/ggplot2)

## Installation

### Dependencies
[![CI](https://github.com/eernst/final_finalizer/actions/workflows/ci.yml/badge.svg)](https://github.com/eernst/final_finalizer/actions/workflows/ci.yml)

**Required:**
- Python 3.9+
- Python packages:
  - [tomli](https://pypi.org/project/tomli/) - TOML configuration file support (Python <3.11 only)
  - [intervaltree](https://pypi.org/project/intervaltree/) - efficient overlap detection (recommended for performance)
- [miniprot](https://github.com/lh3/miniprot) - protein-to-genome alignment
- [gffread](https://github.com/gpertea/gffread) - GFF3/FASTA processing
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - organelle/rDNA detection

**Optional:**
- [minimap2](https://github.com/lh3/minimap2) or [mm2plus](https://github.com/lh3/mm2plus) - nucleotide synteny QA and read alignment for depth analysis
- [samtools](https://github.com/samtools/samtools) - BAM/CRAM handling for depth analysis
- [mosdepth](https://github.com/brentp/mosdepth) - efficient depth calculation
- [rasusa](https://github.com/mbhall88/rasusa) - FASTQ downsampling for depth analysis
- [centrifuger](https://github.com/mourisl/centrifuger) - contaminant detection
- [taxonkit](https://github.com/shenwei356/taxonkit) + NCBI taxonomy database - taxonomic lineage for contaminant visualization (see below)
- R with ggplot2, dplyr, readr, stringr, tibble, tidyr, patchwork, ggnewscale, pacman - visualization (`--plot`)
- R with treemapify - contaminant phylogenetic breakdown visualization
- R with ggiraph, htmlwidgets, pandoc - interactive visualization (`--plot-html`)

### Conda environment

```bash
conda create -n final_finalizer python=3.10 miniprot gffread blast mm2plus centrifuger -c bioconda -c conda-forge
conda activate final_finalizer

# Install Python dependencies
pip install -r requirements.txt
```

Optional (read depth analysis):
```bash
conda install -n final_finalizer -c bioconda samtools mosdepth rasusa
```

Optional (plotting):
```bash
conda install -n final_finalizer -c conda-forge r-base r-ggplot2 r-dplyr r-readr r-stringr r-tibble r-tidyr r-patchwork r-ggnewscale r-pacman r-ggiraph r-htmlwidgets r-scales libxml2 pandoc
```

Note: The contaminant treemap visualization requires the R `treemapify` package, which depends on `libxml2`. After installing `libxml2` via conda, you may need to create a symlink for R to find it:
```bash
ln -sf libxml2.so.16 $CONDA_PREFIX/lib/libxml2.so
```

Optional (taxonkit for contaminant phylogenetic visualization):

For full Domain → Family → Genus → Species breakdown in the contaminant treemap, install taxonkit and the NCBI taxonomy database:

```bash
# Install taxonkit
conda install -n final_finalizer -c bioconda taxonkit

# Download and set up NCBI taxonomy database
# Option 1: Let taxonkit download it (requires ~1.5GB)
taxonkit download --data-dir ~/.taxonkit

# Option 2: Manual download
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p ~/.taxonkit
tar -xzf taxdump.tar.gz -C ~/.taxonkit
```

Without taxonkit, the contaminant plot will show a simpler treemap by genus (parsed from scientific names). With taxonkit and the taxonomy database, you get a full hierarchical treemap showing the phylogenetic breakdown of detected contaminants (Domain > Family > Genus > Species).

Latest tested conda package versions (CI):
<!-- conda-versions-start -->
- python: 3.10
- miniprot: unknown
- gffread: unknown
- blast: unknown
- minimap2: unknown
- centrifuger: unknown
<!-- conda-versions-end -->

## Quick Start

```bash
./final_finalizer.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o output_prefix \
    --ref-gff3 reference.gff3 \
    --plot \
    --plot-html
```

## Usage

### Required arguments

| Argument | Description |
|----------|-------------|
| `-r, --ref` | Reference genome FASTA (can be gzipped) |
| `-q, --query` | Query assembly FASTA to classify |
| `-o, --outprefix` | Output file prefix |
| `--ref-gff3` | Reference GFF3 with protein-coding gene annotations. Required for `--synteny-mode protein`. |

### Common options

| Argument | Description | Default |
|----------|-------------|---------|
| `-t, --threads` | Number of threads | 8 |
| `--plot` | Generate PDF visualization | off |
| `--plot-html` | Also generate interactive HTML visualization | off |
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
- `lrhq` → `-x lr:hqae` (PacBio HiFi/Duplex, ONT Q20+; high-quality long reads with error rate < 1%)
- `r9` → `-x map-ont` (ONT reads, standard accuracy)
- `sr` → `-x sr` (Illumina short reads)

**Read downsampling**: If `--depth-target-coverage` is set to a non-zero value (e.g., 20 for 20X coverage), reads are automatically downsampled before alignment using rasusa (for FASTQ) or samtools (for BAM/CRAM). This reduces computational time while maintaining sufficient coverage for depth-based quality assessment. Pre-aligned BAM/CRAM files are not downsampled. Default: 0 (disabled).

**Depth caching**: Alignment and depth results are cached and automatically reused on subsequent runs if inputs match. Use `--keep-depth-bam` to retain the aligned BAM file after depth calculation (default: deleted to save space). Cached depth results can be reused even if the BAM is deleted.

### Synteny mode selection

| Argument | Description | Default |
|----------|-------------|---------|
| `--synteny-mode` | Synteny evidence source: `protein` (miniprot) or `nucleotide` (minimap2) | protein |

**Protein mode** (default): Uses miniprot protein-anchored synteny. Requires `--ref-gff3`. Ideal for detecting conserved gene content across distantly related species.

**Nucleotide mode**: Uses minimap2 whole-genome nucleotide alignment with permissive chaining for chromosome-scale structural composition analysis. Creates megabase-scale synteny blocks by chaining through repetitive regions. Suitable for both within-species and cross-species comparisons. Ideal for identifying structural features like chromosomal fusions, homeologous recombination, or introgression events.

### Pipeline toggles

| Argument | Description |
|----------|-------------|
| `--skip-organelles` | Skip organelle detection |
| `--skip-rdna` | Skip rDNA detection |
| `--skip-contaminants` | Skip contaminant detection |

## Output Files

### FASTA outputs

| File | Description |
|------|-------------|
| `*.chrs.fasta` | Chromosome-assigned contigs (reoriented if needed) |
| `*.organelles.fasta` | Organelle contigs (chrC, chrM) |
| `*.rdna.fasta` | rDNA-containing contigs |
| `*.contaminants.fasta` | Contaminant contigs |
| `*.debris.fasta` | Assembly debris (fragments, duplicates) |
| `*.unclassified.fasta` | Contigs that couldn't be classified |

### Summary tables

| File | Description |
|------|-------------|
| `*.contig_summary.tsv` | Per-contig classification with metrics |
| `*.evidence_summary.tsv` | Synteny evidence per contig × reference |
| `*.segments.tsv` | Individual synteny segments |
| `*.macro_blocks.tsv` | Aggregated synteny macro-blocks |
| `*.ref_lengths.tsv` | Reference chromosome lengths |

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

| File | Description |
|------|-------------|
| `*.chromosome_overview.pdf` | Multi-panel plot showing contig composition, subgenome support, and alignment identity |
| `*.depth_overview.pdf` | Read depth visualization by classification and chromosome (if `--reads` provided) |
| `*.depth_overview.html` | Interactive version with tooltips (if `--plot-html` and `--reads` provided) |
| `*.contaminant_treemap.pdf` | Phylogenetic breakdown of contaminants as hierarchical treemap (if contaminants detected) |
| `*.contaminant_bandage.pdf` | Individual contaminant contigs visualized as Bandage-style shapes colored by taxonomy (if contaminants detected) |
| `*.contaminants.tsv` | Detailed contaminant summary with taxonomic lineage |

## Classification Pipeline

The tool runs these phases in order:

1. **Reference preparation** - Read reference genome and compute GC statistics
2. **Synteny analysis** (mode-dependent):
   - **Protein mode**: Extract proteins from GFF3 (gffread), align to query assembly (miniprot)
   - **Nucleotide mode**: Whole-genome alignment with permissive chaining (minimap2)
3. **Synteny block building** - Chain alignments into synteny blocks; identify chromosome candidates
4. **Organelle detection** - BLAST non-chromosome contigs against chrC/chrM references
5. **rDNA detection** - BLAST against rDNA reference
6. **Chromosome debris detection** - High-coverage, high-identity matches to assembled chromosomes
7. **Contaminant detection** - Centrifuger taxonomic classification with two-gate filtering (score ≥1000, coverage ≥0.50)
8. **Debris classification** - Reference-based debris detection for remaining contigs
9. **Gene count statistics** (if GFF3 provided) - Compute gene proportion metrics
10. **Orientation determination** - Determine strand for chromosome contigs based on synteny votes
11. **Final classification** - Assign all contigs to categories with confidence levels
12. **Read depth analysis** (optional) - Align reads and compute per-contig depth metrics
13. **Output generation** - Write classified FASTAs, summary tables, and visualizations

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

### Confidence Criteria by Category

| Classification | High | Medium | Low |
|----------------|------|--------|-----|
| **chrom_assigned** | Gene proportion ≥20% AND GC deviation <2σ | Gene proportion 10-20% OR GC deviation 2-3σ | Gene proportion <10% OR GC deviation >3σ |
| **chrom_unassigned** | — | GC deviation <2σ | GC deviation ≥2σ |
| **organelle_complete** | Coverage ≥90% | Coverage 80-90% | — |
| **organelle_debris** | — | Coverage ≥60% | Coverage <60% |
| **rDNA** | Coverage ≥80% AND identity ≥95% | Coverage 50-80% OR identity <95% | Coverage <60% |
| **contaminant** | Coverage ≥80% OR GC deviation >2σ | Coverage 50-80% | Coverage <50% |
| **chrom_debris** | Coverage ≥90%, identity ≥95%, GC deviation <3σ | Coverage 80-90% OR identity 90-95% OR GC deviation ≥3σ | — |
| **debris** | Coverage ≥80% OR ≥5 protein hits | Coverage 50-80% OR 2-5 protein hits | GC deviation >3σ AND no synteny |
| **unclassified** | — | — | Always (no evidence) |

### Evidence Factors

**Gene proportion** (`ref_gene_proportion`): Fraction of reference chromosome genes aligned to the query contig (e.g., 0.85 = 85% of reference genes found on this contig). Higher values indicate stronger synteny support. This metric assesses completeness relative to the reference chromosome.

**GC deviation**: How many standard deviations (σ) the contig's GC content differs from a GC baseline. Large deviations may indicate contamination or unusual sequences.
- For `chrom_assigned`: compared to reference nuclear genome GC (validates synteny-based assignment)
- For all other categories: compared to assembly chromosome GC (more appropriate for divergent genomes where the assembly may differ significantly from the reference)

**Coverage**: Fraction of contig covered by alignments (0.0-1.0). Higher coverage indicates more complete matches.

**Identity**: Alignment identity (0.0-1.0). Values ≥0.95 indicate very high sequence similarity.

**Protein hits**: Number of unique reference protein-coding genes with miniprot alignments to the contig.

## Thresholds

### Chromosome assignment

**Common thresholds (both modes):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assign-min-frac` | 0.10 | Min synteny coverage of contig |
| `--assign-min-ratio` | 1.25 | Min best/second score ratio |
| `--miniprot-min-span-frac` | 0.20 | Min span fraction of contig (applies to both modes) |
| `--miniprot-min-span-bp` | 50000 | Min absolute span in bp (applies to both modes) |

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

All gate criteria must be satisfied for chromosome assignment (AND logic).

**Note**: Despite the `miniprot-` prefix on some parameters, span fraction and span bp thresholds apply to both modes. The segment count threshold differs: protein mode requires ≥5 segments (configurable), while nucleotide mode requires ≥1 segment (hardcoded) because perfect full-length alignments produce fewer segments.

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

`final_finalizer` supports TOML configuration files for managing complex parameter sets and reproducible analysis workflows.

### Generating a config template

```bash
./final_finalizer.py --dump-config > my_config.toml
```

This creates a complete configuration file with all parameters and their current default values. Edit this file to customize your analysis.

### Using a config file

```bash
./final_finalizer.py --config my_config.toml
```

**Important:** Command-line arguments override config file values. This allows you to have a base configuration but adjust specific parameters via CLI:

```bash
# Use config but override threads
./final_finalizer.py --config my_config.toml -t 64
```

### Example config file

```toml
# Required
[required]
ref = "/path/to/reference.fasta"
query = "/path/to/assembly.fasta"
outprefix = "output_prefix"
ref_gff3 = "/path/to/reference.gff3"

# Common Options
[common]
threads = 32
plot = true
plot_html = true
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
./final_finalizer.py \
    -r TAIR10.fasta \
    -q my_assembly.fasta \
    -o my_assembly_classified \
    --ref-gff3 TAIR10.gff3 \
    -t 32 \
    --plot \
    --plot-html
```

### Nucleotide mode for structural composition analysis

Use nucleotide mode when you want to detect chromosome-scale structural features like fusions, translocations, or homeologous exchanges. This mode uses whole-genome nucleotide alignment and is suitable for both within-species and cross-species comparisons.

```bash
./final_finalizer.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o assembly_nucleotide \
    --synteny-mode nucleotide \
    -t 32 \
    --plot \
    --plot-html
```

Note: Nucleotide mode does not require `--ref-gff3`, but if provided, gene count statistics will be included in the output.

### Polyploid genome with contaminant screening

```bash
./final_finalizer.py \
    -r wheat_ref.fasta \
    -q wheat_assembly.fasta \
    -o wheat_classified \
    --ref-gff3 wheat_ref.gff3 \
    --centrifuger-idx /path/to/centrifuger/nt \
    --rdna-ref default \
    -t 64 \
    --plot \
    --plot-html
```

### Non-polyploid with subgenome suffix

```bash
./final_finalizer.py \
    -r rice_ref.fasta \
    -q rice_assembly.fasta \
    -o rice_classified \
    --ref-gff3 rice_ref.gff3 \
    --add-subgenome-suffix A \
    --plot \
    --plot-html
```

### With read depth analysis (HiFi reads)

```bash
./final_finalizer.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o assembly_classified \
    --ref-gff3 reference.gff3 \
    --reads hifi_reads.fastq.gz \
    --reads-type lrhq \
    --plot \
    --plot-html \
    -t 32
```

### With pre-aligned BAM for depth analysis

```bash
./final_finalizer.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o assembly_classified \
    --ref-gff3 reference.gff3 \
    --reads aligned_reads.bam \
    --plot \
    --plot-html \
    -t 32
```

### Using configuration file with verbose logging

```bash
# Generate config template
./final_finalizer.py --dump-config > wheat_config.toml

# Edit wheat_config.toml to set paths and parameters

# Run with config and verbose logging
./final_finalizer.py \
    --config wheat_config.toml \
    --verbose \
    --log-file wheat_analysis.log
```

### With downsampled reads for faster depth analysis

```bash
./final_finalizer.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o assembly_classified \
    --ref-gff3 reference.gff3 \
    --reads hifi_reads.fastq.gz \
    --reads-type lrhq \
    --depth-target-coverage 20 \
    --keep-depth-bam \
    --plot \
    -t 32
```

## Algorithm Details

### Synteny Mode Selection

The tool supports two complementary synteny modes:

#### Protein mode (default)

Uses protein homology as the primary evidence source because:
- Proteins are more conserved than nucleotide sequences
- Works across distantly related species
- Robust to repetitive sequences
- Provides gene-level resolution for functional assessment

Miniprot alignments are:
1. Filtered by alignment quality and identity
2. Mapped to reference genomic coordinates via GFF3
3. Overlapping hits filtered by identity (keeps highest-identity hit at each position)
4. Chained into collinear synteny blocks
5. Aggregated per contig × reference chromosome
6. Scored and ranked for assignment

**Performance optimization**: The pipeline uses interval trees (via the `intervaltree` package) for efficient O(n log n) overlap detection during synteny block filtering, replacing the previous O(n²) algorithm. This significantly improves runtime for large assemblies with many protein alignments.

#### Nucleotide mode

Uses whole-genome nucleotide alignment for structural composition analysis:
- Detects actual sequence-level identity and synteny
- Ideal for identifying structural features (fusions, translocations, homeologous exchanges)
- Works for both within-species and cross-species comparisons
- Creates megabase-scale synteny blocks for chromosome architecture analysis

Minimap2 alignments use permissive chaining parameters to create continuous chromosome-scale blocks:
- `--max-chain-skip=300`: Chain through repetitive regions
- `-z 10000,1000`: Tolerate long gaps in alignment chains
- `-r 50000`: Large bandwidth (50kb) to accommodate structural variations

This permissive approach chains through homopolymer runs, tandem repeats, and ambiguous regions that would otherwise fragment alignments. The resulting megabase-scale blocks are suitable for chromosome classification and compositional analysis, balanced by downstream filtering (identity thresholds, minimum alignment length, gate filtering) to prevent spurious assignments.

**Gate threshold differences**: Nucleotide mode requires ≥1 segment (vs. ≥5 for protein mode) because perfect full-length nucleotide alignments produce fewer segments than fragmented protein hits.

### Gate-based assignment

Contigs must pass "gates" to be assigned:
- Minimum number of synteny segments
- Minimum contig span covered
- Minimum unique genes aligned

This prevents spurious assignments from:
- Single conserved genes
- Repetitive sequences
- Low-complexity regions

### Debris detection

Two complementary approaches:

1. **Chromosome debris** - Contigs with high-coverage (≥80%) and high-identity (≥90%) matches to already-assembled chromosome contigs (detected via minimap2 asm5 alignment)

2. **Reference-based debris** - Contigs with moderate reference coverage (>50%) but insufficient synteny support for chromosome assignment

### Contaminant visualization

When contaminants are detected (with `--centrifuger-idx`) and plotting is enabled (`--plot`), final_finalizer generates two complementary visualizations:

1. **Phylogenetic treemap** (`*.contaminant_treemap.pdf`): Hierarchical treemap showing taxonomic breakdown (Kingdom → Family → Genus → Species) with area proportional to total contamination span. Requires taxonkit for full taxonomic lineage; falls back to genus-level grouping if unavailable.

2. **Bandage-style plot** (`*.contaminant_bandage.pdf`): Individual contaminant contigs visualized as geometric shapes:
   - Circular contigs (names ending in "c"): Drawn as rings/donuts
   - Linear contigs (names ending in "l"): Drawn as pills (elongated ellipses)
   - Size uses log-scale areas so small contigs remain visible
   - Colored by taxonomic family
   - Ordered by decreasing read depth (if `--reads` provided)

Both visualizations filter to high-confidence contaminants (coverage ≥ `--contaminant-min-coverage`, default 0.50) to reduce noise from conserved gene matches.

## Citation

If you use this tool, please cite:
- [miniprot](https://github.com/lh3/miniprot)
- [minimap2](https://github.com/lh3/minimap2)
- [gffread](https://github.com/gpertea/gffread)
- [centrifuger](https://github.com/mourisl/centrifuger) (if used)
- [taxonkit](https://github.com/shenwei356/taxonkit) (if used for contaminant visualization)

## License

MIT License

## Author

Evan Ernst (eernst@cshl.edu)
