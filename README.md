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

**Required:**
- Python 3.8+
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
conda create -n final_finalizer python=3.10 miniprot gffread blast minimap2 centrifuger -c bioconda -c conda-forge
conda activate final_finalizer
```

## Quick Start

```bash
./final_finalizer.py \
    -r reference.fasta \
    -q assembly.fasta \
    -o output_prefix \
    --ref-gff3 reference.gff3 \
    --plot
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
| `-C, --chr-like-minlen` | Min contig length for chromosome classification | 80% of smallest ref chromosome |
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
| `classification` | Category: chrom, organelle_complete, rDNA, contaminant, debris, unclassified |
| `reversed` | Whether contig was reverse-complemented |
| `contaminant_taxid` | NCBI taxonomy ID (for contaminants) |
| `contaminant_sci` | Scientific name (for contaminants) |
| `assigned_ref_id` | Best-matching reference chromosome |
| `genes_per_Mbp` | Gene density for chromosome contigs |

### Visualization

| File | Description |
|------|-------------|
| `*.chromosome_overview.pdf` | Multi-panel plot showing contig composition, subgenome support, and alignment identity |

## Classification Pipeline

The tool runs these phases in order:

1. **Reference protein extraction** - Extract proteins from GFF3 using gffread
2. **Protein-anchor synteny** - Align proteins to query assembly (miniprot)
3. **Synteny block building** - Chain alignments into synteny blocks
4. **Chromosome assignment** - Assign contigs to reference chromosomes
5. **Organelle detection** - BLAST against chrC/chrM references
6. **rDNA detection** - BLAST against rDNA reference
7. **Chromosome debris detection** - High-identity matches to assembled chromosomes
8. **Contaminant detection** - Centrifuger taxonomic classification
9. **Debris classification** - Reference/protein-based debris detection
10. **Orientation determination** - Determine strand for chromosome contigs
11. **Output generation** - Write classified FASTAs and summary tables

## Contig Classification Categories

| Category | Description |
|----------|-------------|
| `chrom` | Chromosome-like contig assigned to a reference chromosome |
| `organelle_complete` | Complete organelle genome (chrC or chrM) |
| `organelle_debris` | Partial organelle sequence |
| `rDNA` | Ribosomal DNA repeat unit |
| `contaminant` | Sequence from contaminating organism |
| `chrom_debris` | Fragment of assembled chromosome (high-identity duplicate) |
| `debris` | Assembly debris with protein/reference support |
| `unclassified` | Could not be classified |

## Thresholds

### Chromosome assignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assign-min-frac` | 0.10 | Min synteny coverage of contig |
| `--assign-min-ratio` | 1.25 | Min best/second score ratio |
| `--miniprot-min-genes` | 3 | Min unique genes for assignment |
| `--miniprot-min-segments` | 5 | Min synteny segments |
| `--miniprot-min-span-frac` | 0.20 | Min span fraction of contig |

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
    --plot
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
    --plot
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
    --plot
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

1. **Chromosome debris** - Contigs with high-identity matches (>90%) to already-assembled chromosome contigs (detected via minimap2 asm5 alignment)

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
