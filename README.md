# final_finalizer

**Genome assembly finalization tool for contig classification and quality control**

`final_finalizer` classifies contigs from a *de novo* genome assembly into biological categories (chromosomes, organelles, rDNA, contaminants, debris, unclassified) using protein-anchored synteny evidence, organelle/rDNA alignments, and taxonomic classification.

## Features

- **Chromosome assignment** using protein-anchored synteny blocks (miniprot)
- **Subgenome resolution** for polyploid genomes with multiple chromosome sets
- **Organelle detection** (chloroplast/chrC, mitochondrion/chrM) via BLAST
- **rDNA contig identification** via BLAST against reference rDNA
- **Contaminant screening** via centrifuger taxonomic classification
- **Debris detection** for assembly fragments (chromosome debris, organelle debris)
- **Chimera flagging** for contigs with evidence from multiple chromosomes
- **Orientation determination** for chromosome-assigned contigs
- **Publication-ready visualizations** (PDF plots via R/ggplot2)

## Installation

### Dependencies
[![CI](https://github.com/eernst/final_finalizer/actions/workflows/ci.yml/badge.svg)](https://github.com/eernst/final_finalizer/actions/workflows/ci.yml)

**Required:**
- Python 3.9+
- [miniprot](https://github.com/lh3/miniprot) - protein-to-genome alignment
- [gffread](https://github.com/gpertea/gffread) - GFF3/FASTA processing
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - organelle/rDNA detection

**Optional:**
- [minimap2](https://github.com/lh3/minimap2) or [mm2plus](https://github.com/lh3/mm2plus) - nucleotide synteny QA
- [centrifuger](https://github.com/mourisl/centrifuger) - contaminant detection
- R with ggplot2, dplyr, readr, stringr, tibble, tidyr, patchwork, ggnewscale, pacman - visualization (`--plot`)
- R with ggiraph, htmlwidgets, pandoc - interactive visualization (`--plot-html`)

### Conda environment

```bash
conda create -n final_finalizer python=3.10 miniprot gffread blast mm2plus centrifuger -c bioconda -c conda-forge
conda activate final_finalizer
```

Optional (plotting):
```bash
conda install -n final_finalizer -c conda-forge r-base r-ggplot2 r-dplyr r-readr r-stringr r-tibble r-tidyr r-patchwork r-ggnewscale r-pacman r-ggiraph r-htmlwidgets pandoc
```

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
| `--ref-gff3` | Reference GFF3 with protein-coding gene annotations |

### Common options

| Argument | Description | Default |
|----------|-------------|---------|
| `-t, --threads` | Number of threads | 8 |
| `--plot` | Generate PDF visualization | off |
| `--plot-html` | Also generate interactive HTML visualization | off |
| `-C, --chr-like-minlen` | Min contig length for chromosome classification | 80% of smallest nuclear ref chromosome |
| `--add-subgenome-suffix` | Suffix for non-polyploid references (e.g., 'A') | none |

### Classification references

| Argument | Description |
|----------|-------------|
| `--chrC-ref` | Chloroplast reference FASTA (default: extract from --ref) |
| `--chrM-ref` | Mitochondrion reference FASTA (default: extract from --ref) |
| `--rdna-ref` | rDNA reference FASTA, or 'default' for bundled Arabidopsis 45S |
| `--centrifuger-idx` | Centrifuger index prefix for contaminant screening |

### Pipeline toggles

| Argument | Description |
|----------|-------------|
| `--nt-synteny` | Run nucleotide synteny alignments (minimap2) for QA |
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
| `genes_per_Mbp` | Gene density (aligned reference genes per Mbp of query sequence) |

For complete column documentation for all TSV files, see [docs/output_formats.md](docs/output_formats.md).

### Visualization

| File | Description |
|------|-------------|
| `*.chromosome_overview.pdf` | Multi-panel plot showing contig composition, subgenome support, and alignment identity |

## Classification Pipeline

The tool runs these phases in order:

1. **Reference protein extraction** - Extract proteins from GFF3 using gffread
2. **Protein-anchor synteny** - Align proteins to query assembly (miniprot)
3. **Synteny block building** - Chain alignments into synteny blocks; identify chromosome candidates
4. **Organelle detection** - BLAST non-chromosome contigs against chrC/chrM references
5. **rDNA detection** - BLAST against rDNA reference
6. **Chromosome debris detection** - High-coverage, high-identity matches to assembled chromosomes
7. **Contaminant detection** - Centrifuger taxonomic classification
8. **Debris classification** - Reference/protein-based debris detection for remaining contigs
9. **Orientation determination** - Determine strand for chromosome contigs based on synteny votes
10. **Final classification** - Assign all contigs to categories
11. **Output generation** - Write classified FASTAs and summary tables

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

**Gene proportion**: Fraction of reference chromosome genes aligned to the query contig. Higher values indicate stronger synteny support.

**GC deviation**: How many standard deviations (σ) the contig's GC content differs from a GC baseline. Large deviations may indicate contamination or unusual sequences.
- For `chrom_assigned`: compared to reference nuclear genome GC (validates synteny-based assignment)
- For all other categories: compared to assembly chromosome GC (more appropriate for divergent genomes where the assembly may differ significantly from the reference)

**Coverage**: Fraction of contig covered by alignments (0.0-1.0). Higher coverage indicates more complete matches.

**Identity**: Alignment identity (0.0-1.0). Values ≥0.95 indicate very high sequence similarity.

**Protein hits**: Number of unique reference protein-coding genes with miniprot alignments to the contig.

## Thresholds

### Chromosome assignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assign-min-frac` | 0.10 | Min synteny coverage of contig |
| `--assign-min-ratio` | 1.25 | Min best/second score ratio |
| `--miniprot-min-genes` | 3 | Min unique genes for assignment |
| `--miniprot-min-segments` | 5 | Min synteny segments |
| `--miniprot-min-span-frac` | 0.20 | Min span fraction of contig |

All gate criteria must be satisfied for chromosome assignment (AND logic).

### Organelle detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--organelle-min-cov` | 0.80 | Min query coverage |
| `--chrC-len-tolerance` | 0.05 | Length tolerance for chrC |
| `--chrM-len-tolerance` | 0.20 | Length tolerance for chrM |

### Contaminant detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--contaminant-min-score` | 150 | Min centrifuger score |

### Debris detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--chr-debris-min-cov` | 0.80 | Min coverage vs chromosomes |
| `--chr-debris-min-identity` | 0.90 | Min identity vs chromosomes |
| `--debris-min-cov` | 0.50 | Min coverage vs reference |

## Example Workflows

### Basic chromosome assignment

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

## Algorithm Details

### Protein-anchor synteny

The tool uses protein homology as the primary evidence source because:
- Proteins are more conserved than nucleotide sequences
- Works across distantly related species
- Robust to repetitive sequences
- Provides gene-level resolution

Miniprot alignments are:
1. Filtered by alignment quality
2. Mapped to reference genomic coordinates via GFF3
3. Chained into collinear synteny blocks
4. Aggregated per contig × reference chromosome
5. Scored and ranked for assignment

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

## Citation

If you use this tool, please cite:
- [miniprot](https://github.com/lh3/miniprot)
- [minimap2](https://github.com/lh3/minimap2)
- [gffread](https://github.com/gpertea/gffread)
- [centrifuger](https://github.com/mourisl/centrifuger) (if used)

## License

MIT License

## Author

Evan Ernst (eernst@cshl.edu)
