# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`final_finalizer` is a bioinformatics genome assembly finalization tool that classifies contigs from *de novo* genome assemblies into biological categories (chromosomes, organelles, rDNA, contaminants, debris) using protein-anchored synteny evidence, organelle/rDNA alignments, and taxonomic classification.

**Key concept**: The tool uses protein homology (via miniprot) as the primary evidence source because proteins are more conserved than nucleotide sequences, work across distantly related species, and are robust to repetitive sequences.

## Development Commands

### Running tests
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_final_finalizer.py

# Run specific test function
pytest tests/test_final_finalizer.py::test_parse_paf_chain_evidence_for_subgenome_assignment

# Run with coverage
pytest --cov=final_finalizer --cov-report=html

# Run only integration tests (may download large reference files)
pytest -m integration

# Skip integration tests
pytest -m "not integration"
```

### Running the tool
```bash
# Basic run
./final_finalizer.py -r reference.fasta -q assembly.fasta -o output --ref-gff3 reference.gff3

# With all features enabled
./final_finalizer.py -r ref.fasta -q assembly.fasta -o output \
    --ref-gff3 ref.gff3 \
    --reads hifi.fastq.gz \
    --centrifuger-idx /path/to/index \
    --plot --plot-html \
    -v --log-file analysis.log

# Dump config template
./final_finalizer.py --dump-config > config.toml

# Run with config file
./final_finalizer.py --config config.toml
```

### Code style
```bash
# No formal linters configured - follow existing code style
# Use Python 3.9+ type hints with `from __future__ import annotations`
```

## Architecture Overview

### Package Structure

```
final_finalizer/
├── models.py              # Data models (Block, Chain, ContigClassification, etc.)
├── cli.py                 # Main entry point and argument parsing
├── alignment/             # Alignment and synteny block building
│   ├── external_tools.py  # Wrappers for miniprot, minimap2, gffread, BLAST
│   ├── chain_parsing.py   # PAF/miniprot parsing → synteny blocks (O(n log n) interval tree)
│   └── stats.py           # Alignment statistics
├── classification/        # Contig classification logic
│   └── classifier.py      # Gate-based assignment, orientation, confidence levels
├── detection/             # Feature detection modules
│   ├── organelle.py       # chrC/chrM detection via BLAST
│   ├── rdna.py            # rDNA detection via BLAST
│   ├── debris.py          # Chromosome debris detection via minimap2
│   └── contaminant.py     # Centrifuger taxonomic classification
├── analysis/              # Read depth analysis
│   └── read_depth.py      # Depth calculation via mosdepth with caching
├── output/                # Output generation
│   ├── fasta_output.py    # Write classified FASTA files
│   ├── tsv_output.py      # TSV summary tables
│   └── plotting.py        # R/ggplot2 visualization
└── utils/                 # Shared utilities
    ├── io_utils.py        # File I/O, gzip handling, interval merging
    ├── sequence_utils.py  # FASTA reading/writing, reverse complement
    ├── reference_utils.py # Reference chromosome ID normalization
    ├── logging_config.py  # Custom logging with [info]/[warn]/[error] format
    └── config.py          # TOML configuration file support
```

### Data Flow

1. **Reference Preparation** (`cli.py:main`)
   - Extract proteins from GFF3 using gffread
   - Build reference ID maps (handles polyploid naming like chr1A/chr1B)
   - Prepare organelle references (chrC/chrM)

2. **Protein-Anchor Synteny** (`alignment/`)
   - Run miniprot to align proteins to query assembly
   - Parse PAF output in `chain_parsing.py`
   - Filter overlapping hits using **interval tree** (O(n log n), not O(n²))
   - Chain alignments into synteny blocks
   - Aggregate evidence per contig × reference chromosome

3. **Detection Pipeline** (`detection/`)
   - Organelle detection via BLAST (chrC, chrM)
   - rDNA detection via BLAST against reference rDNA
   - Chromosome debris via minimap2 (high coverage/identity duplicates)
   - Contaminant screening via centrifuger

4. **Classification** (`classification/classifier.py`)
   - **Gate-based assignment**: contigs must pass multiple criteria (min genes, segments, span) to be assigned
   - Orientation determination via synteny strand votes
   - Confidence level assignment (high/medium/low) based on multiple evidence factors
   - Debris/unclassified handling for remaining contigs

5. **Read Depth Analysis** (`analysis/read_depth.py`) - OPTIONAL
   - Auto-detect read format (FASTQ, BAM, CRAM)
   - Align reads with minimap2 (if not pre-aligned)
   - Calculate depth with mosdepth
   - **Caching**: depth results cached with metadata validation
   - Optional BAM cleanup (default: delete BAM after depth calculation)

6. **Output Generation** (`output/`)
   - Classified FASTA files (chromosomes, organelles, rDNA, etc.)
   - TSV summary tables (contig_summary, evidence_summary, segments, macro_blocks)
   - PDF/HTML visualizations (if R/ggplot2 available)

### Key Design Patterns

**Gate-based assignment**: Contigs must satisfy ALL criteria for chromosome assignment (AND logic). This prevents spurious assignments from single conserved genes or repetitive sequences. See `--assign-min-frac`, `--miniprot-min-genes`, `--miniprot-min-segments`, `--miniprot-min-span-frac`.

**Interval tree optimization**: `chain_parsing.py:_filter_overlapping_hits_by_identity()` uses the `intervaltree` package for O(n log n) overlap detection instead of O(n²). This is critical for large assemblies with many protein alignments.

**Reference ID normalization**: `reference_utils.py` handles polyploid naming conventions (chr1A/chr1B) and organelle aliases (chrC/Pt, chrM/Mt). Configurable via `set_ref_id_patterns()`.

**Depth caching**: `read_depth.py` caches alignment and depth results with metadata validation (query MD5, window size, etc.) to avoid re-running expensive alignment steps. Cache is automatically invalidated if inputs change.

**Logging infrastructure**: Custom formatter in `logging_config.py` provides `[info]`, `[warn]`, `[error]`, `[done]` prefixes. Use `logger.phase()` for major pipeline phases and `logger.done()` for completion messages.

**TOML configuration**: `config.py` provides schema-based config file support. CLI arguments override config file values.

### Critical Security Notes

**Command injection prevention**: `alignment/external_tools.py` validates extra arguments against shell metacharacters using regex pattern `[;&|` `$(){}[\]<>!\\'"#]`. This prevents command injection via user-provided arguments like `--miniprot-extra-args`.

**TOCTOU race condition fix**: `utils/io_utils.py:file_exists_and_valid()` uses a single `stat()` call instead of separate `exists()` + `getsize()` to avoid time-of-check-time-of-use vulnerabilities.

**Process cleanup**: After killing subprocesses, code calls `kill()` then `wait()` to reap zombie processes and prevent resource leaks.

## Common Gotchas

**GC deviation returns None**: `classifier.py:compute_gc_deviation()` returns `None` (not 0.0) when reference GC std < 0.001. This prevents division-by-zero and ensures low-confidence classification for contigs when GC comparison is meaningless.

**Subsample fraction validation**: `read_depth.py` validates `--depth-target-coverage` is non-negative. Invalid fractions raise `ValueError` with clear error messages.

**Organelle reference extraction**: If `--chrC-ref` or `--chrM-ref` not provided, code attempts to extract from main reference FASTA. If organelle sequences not found, organelle detection is skipped (not an error).

**PAF vs miniprot parsing**: The tool supports both nucleotide synteny (minimap2 PAF) and protein-anchor synteny (miniprot PAF). `chain_parsing.py` has separate functions: `parse_paf_chain_evidence_and_segments()` for minimap2 and `parse_miniprot_synteny_evidence_and_segments()` for miniprot.

**Reference length normalization**: `reference_utils.py:normalize_ref_lengths()` filters out organelles (chrC/chrM) and unlocalized scaffolds (chr*_random) from nuclear chromosome length calculations. This ensures `--chr-like-minlen` thresholds are computed correctly.

## Testing Philosophy

Tests are in `tests/` directory. 21 tests total:
- Unit tests: Test individual functions with minimal dependencies
- Integration tests: Marked with `@pytest.mark.integration`, may download large reference files (Arabidopsis TAIR10)

**Integration test caching**: Downloaded reference files are cached in `data/` directory to speed up subsequent test runs.

**Test data**: Tests use real biological data (Arabidopsis TAIR10 reference) to ensure the tool works with real-world inputs, not just synthetic test cases.

## External Dependencies

**Required command-line tools**:
- miniprot - protein-to-genome alignment
- gffread - GFF3/FASTA processing
- blastn - organelle/rDNA detection

**Optional command-line tools**:
- minimap2/mm2plus - nucleotide synteny and read alignment
- samtools - BAM/CRAM handling
- mosdepth - depth calculation
- rasusa - FASTQ downsampling
- centrifuger - contaminant detection
- Rscript (with ggplot2, dplyr, etc.) - visualization

**Python packages**:
- tomli - TOML config file parsing (Python <3.11 only)
- intervaltree - efficient overlap detection (recommended for performance)

All external tools are called via subprocess with proper error handling. Use `utils/io_utils.py:have_exe()` to check availability before calling.

## Output Files Reference

**FASTA outputs**: `*.chrs.fasta`, `*.organelles.fasta`, `*.rdna.fasta`, `*.contaminants.fasta`, `*.debris.fasta`, `*.unclassified.fasta`

**TSV outputs**:
- `*.contig_summary.tsv` - Per-contig classification with all evidence fields
- `*.evidence_summary.tsv` - Synteny evidence per contig × reference
- `*.segments.tsv` - Individual synteny segments from chain parsing
- `*.macro_blocks.tsv` - Aggregated synteny macro-blocks
- `*.ref_lengths.tsv` - Reference chromosome lengths

**Visualization**: `*.chromosome_overview.pdf`, `*.depth_overview.pdf`, `*.depth_overview.html` (if `--plot-html`)

See `docs/output_formats.md` for complete column documentation.

## Modifying Classification Logic

Classification happens in `classification/classifier.py:classify_all_contigs()`. The function applies a decision tree:

1. Organelle complete (if BLAST coverage ≥80%)
2. rDNA (if BLAST coverage ≥50%)
3. Chromosome assigned (if passes ALL synteny gates)
4. Chromosome unassigned (if length ≥ `--chr-like-minlen` but failed synteny gates)
5. Contaminant (if centrifuger score ≥ threshold)
6. Chromosome debris (if high coverage/identity vs assembled chromosomes)
7. Organelle debris (if partial organelle match)
8. Debris (if reference coverage >50% or protein hits ≥2)
9. Unclassified (no evidence)

**Confidence levels** are assigned by `classifier.py:assign_classification_confidence()` based on category-specific criteria (gene proportion, GC deviation, coverage, identity, protein hits).

## When Modifying Synteny Block Building

The synteny block building algorithm in `chain_parsing.py` is critical for performance. Key optimizations:

- **Interval tree filtering**: Uses `intervaltree` package for O(n log n) overlap detection
- **Fallback to O(n²)**: If `intervaltree` not installed, falls back to slower nested loop
- **Chain scoring**: Uses `score_topk=10` to consider top-K chains per contig × reference, not just the best chain

If modifying chain parsing, ensure you handle both PAF (minimap2) and miniprot PAF formats correctly. Miniprot PAF includes gene IDs in alignment records.

## Configuration File Schema

`utils/config.py:CONFIG_SCHEMA` maps TOML sections to argparse argument names. To add new configurable parameters:

1. Add CLI argument in `cli.py:parse_args()`
2. Add mapping in `config.py:CONFIG_SCHEMA`
3. Update `--dump-config` template in `config.py:dump_config_template()`

CLI arguments always override config file values.
