# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`final_finalizer` is a bioinformatics genome assembly finalization tool that classifies contigs from *de novo* genome assemblies into biological categories (chromosomes, organelles, rDNA, contaminants, debris) using synteny evidence (protein-anchored or nucleotide whole-genome alignment), organelle/rDNA alignments, and taxonomic classification.

**Key concepts**:
- **Nucleotide mode** (default): Uses whole-genome nucleotide alignment (via minimap2) for structural composition analysis. Detects actual sequence-level identity, making it ideal for identifying structural features like chromosomal fusions, homeologous recombination, or introgression events.
- **Protein mode**: Uses protein homology (via miniprot) as the primary evidence source. Proteins are more conserved than nucleotide sequences, work across distantly related species, and are robust to repetitive sequences. Ideal for detecting conserved gene content.

## Development Commands

### Running tests

**Important:** Use the `final_finalizer` conda environment for testing. Integration tests require external tools (gffread, miniprot, etc.) that are installed in this environment:

```bash
conda activate final_finalizer
```

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
# Basic run (nucleotide mode - default)
# Output goes to output_dir/assembly/assembly.* (name derived from query filename)
./final_finalizer.py -r reference.fasta -q assembly.fasta -o results/

# Protein mode for gene-level classification (requires --ref-gff3)
./final_finalizer.py -r reference.fasta -q assembly.fasta -o results/ --synteny-mode protein --ref-gff3 reference.gff3

# With all features enabled (protein mode)
./final_finalizer.py -r ref.fasta -q assembly.fasta -o results/ \
    --synteny-mode protein --ref-gff3 ref.gff3 \
    --reads hifi.fastq.gz \
    --centrifuger-idx /path/to/index \
    --compleasm-lineage embryophyta \
    -v --log-file analysis.log

# Multi-assembly mode: analyze multiple assemblies against a common reference
# Using a file-of-filenames (TSV with columns: path, name, reads)
./final_finalizer.py -r ref.fasta --ref-gff3 ref.gff3 \
    --fofn assemblies.tsv -o multi_output/

# Multi-assembly mode: scan a directory of FASTA files
./final_finalizer.py -r ref.fasta --ref-gff3 ref.gff3 \
    --assembly-dir /path/to/assemblies/ -o multi_output/

# Distributed mode: submit compute phases as SLURM jobs via executorlib
./final_finalizer.py -r ref.fasta -q assembly.fasta -o results/ --ref-gff3 ref.gff3 \
    --cluster --partition cpuq --max-threads-dist 32 --max-mem-dist 64

# Multi-assembly + distributed: assemblies run concurrently, each phase is a SLURM job
./final_finalizer.py -r ref.fasta --ref-gff3 ref.gff3 \
    --fofn assemblies.tsv -o multi_output/ --cluster

# Dump config template
./final_finalizer.py --dump-config > config.toml

# Run with config file
./final_finalizer.py --config config.toml
```

### Code style
```bash
# No formal linters configured - follow existing code style
# Use Python 3.11+ type hints with `from __future__ import annotations`
```

## Architecture Overview

### Package Structure

```
final_finalizer/
├── models.py              # Data models (ReferenceContext, Block, Chain, ContigClassification, etc.)
├── cli.py                 # Main entry point: prepare_reference(), run_assembly(), main()
├── data/                  # Bundled reference data
│   └── rfam/              # Rfam 15.0 subset: eukaryotic rRNA covariance models (auto-pressed on first use)
├── alignment/             # Alignment and synteny block building
│   ├── external_tools.py  # Wrappers for miniprot, minimap2, gffread
│   ├── chain_parsing.py   # PAF/miniprot parsing → synteny blocks (O(n log n) interval tree)
│   ├── stats.py           # Alignment statistics
│   └── pairwise.py        # Pairwise assembly-vs-assembly synteny (nucleotide mode)
├── classification/        # Contig classification logic
│   └── classifier.py      # Gate-based assignment, orientation, confidence levels
├── detection/             # Feature detection modules
│   ├── blast.py           # BLAST wrappers (makeblastdb, blastn megablast)
│   ├── organelle.py       # chrC/chrM detection via BLAST
│   ├── rdna.py            # rDNA detection via BLAST
│   ├── rdna_consensus.py  # Consensus 45S rDNA building and sub-feature annotation
│   ├── telomere.py        # Telomere repeat detection at contig ends
│   ├── debris.py          # Chromosome debris detection via minimap2
│   ├── contaminant.py     # Centrifuger taxonomic classification
│   └── compleasm.py       # Compleasm (BUSCO completeness) evaluation
├── analysis/              # Read depth analysis
│   └── read_depth.py      # Depth calculation via mosdepth with caching
├── output/                # Output generation
│   ├── fasta_output.py    # Write classified FASTA files
│   ├── tsv_output.py      # TSV summary tables
│   ├── plotting.py        # R/ggplot2 visualization
│   ├── scaffolding.py     # Reference-guided scaffolded chromosome sequences
│   └── comparison.py      # Cross-assembly comparison TSVs
└── utils/                 # Shared utilities
    ├── io_utils.py        # File I/O, gzip handling, interval merging
    ├── sequence_utils.py  # FASTA reading/writing, reverse complement
    ├── reference_utils.py # Reference chromosome ID normalization
    ├── logging_config.py  # Custom logging with [info]/[warn]/[error] format
    ├── config.py          # TOML configuration file support
    ├── multi_assembly.py  # FOFN/directory parsing for multi-assembly mode
    ├── distributed.py     # Executor abstraction (LocalExecutor / executorlib SLURM)
    └── resource_estimation.py  # Per-phase resource estimation for SLURM jobs
```

### Data Flow

1. **Reference Preparation** (`cli.py:prepare_reference()` → `ReferenceContext`)
   - Build reference ID maps (handles polyploid naming like chr1A/chr1B)
   - Compute reference GC baseline
   - Prepare organelle references (chrC/chrM)
   - Prepare rDNA reference
   - Compute `chr_like_minlen`
   - **Protein mode only**: Extract proteins from GFF3 using gffread
   - Returns `ReferenceContext` dataclass shared across all assemblies

2. **Synteny Analysis** (`alignment/`) - **Mode-dependent**
   - **Nucleotide mode** (`--synteny-mode nucleotide`, default):
     - Run minimap2 whole-genome alignment (query vs reference)
     - Parse PAF output in `chain_parsing.py:parse_paf_chain_evidence_and_segments()`
     - Chain alignments into synteny blocks based on nucleotide identity
     - Gate-based filtering requires minimum segments and span (no gene count requirement)
   - **Protein mode** (`--synteny-mode protein`):
     - Run miniprot to align proteins to query assembly
     - Parse PAF output in `chain_parsing.py:parse_miniprot_synteny_evidence_and_segments()`
     - Filter overlapping hits using **interval tree** (O(n log n), not O(n²))
     - Chain alignments into synteny blocks based on protein hits
     - Gate-based filtering requires minimum genes, segments, and span

3. **Detection Pipeline** (`detection/`)
   - Organelle detection via BLAST (chrC, chrM)
   - rDNA detection via BLAST against reference rDNA
   - Chromosome debris via minimap2 (high coverage/identity duplicates)
   - Contaminant screening via centrifuger with two-gate filtering:
     - Score threshold (default 1000, ~1kb matching sequence with k=31)
     - Coverage threshold (default 0.50, filters low-coverage hits that may represent conserved genes rather than contamination)

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
   - TSV summary tables (contig_summary, evidence_summary, segments, macro_blocks, contaminants with taxonomic lineage)
   - Unified HTML report with embedded plots (enabled by default; skip with `--skip-plot`; requires rmarkdown + pandoc):
     - Chromosome overview, classification bar, read depth overview, contaminant table

### Pipeline Phases

`run_assembly()` logs each phase with a numbered prefix. Phases are sequential integers; optional phases are skipped with an info message.

| Phase | Description | Optional? |
|-------|-------------|-----------|
| 1 | Reading query assembly | |
| 2 | Synteny analysis (protein or nucleotide, mode-dependent) | |
| 3 | Organelle detection (BLAST) | `--skip-organelles` |
| 4 | rDNA detection (BLAST) | `--skip-rdna` |
| 5 | Chromosome debris detection (minimap2) | |
| 6 | Contaminant detection (centrifuger) | requires `--centrifuger-idx` |
| 7 | Debris/unclassified classification | |
| 8 | Gene count statistics | requires `--ref-gff3` |
| 9 | Orientation determination | |
| 10 | Telomere detection | `--skip-telomeres` |
| 11 | Classification (assign all contigs) | |
| 12 | Read depth analysis | requires `--reads` |
| 13 | rDNA consensus building | `--skip-rdna-consensus` |
| 14 | Reference-guided scaffolding | `--scaffold` |
| 15 | Writing FASTA outputs | |
| 16 | Writing summary TSV | |
| 17 | Compleasm (BUSCO completeness) evaluation | requires `--compleasm-lineage` / `--skip-compleasm` |

Reference preparation (`prepare_reference()`) runs before phase 1 and is not numbered — it produces the shared `ReferenceContext`. Multi-assembly orchestration (assembly headers, pairwise synteny, comparison report) also runs outside the per-assembly phase sequence.

The rDNA consensus module (`rdna_consensus.py`) uses internal "Step 1-4" numbering within phase 13.

### Key Design Patterns

**Unified assembly architecture**: `cli.py` is structured as three functions: `prepare_reference()` computes shared reference state once (returns `ReferenceContext`), `run_assembly()` runs the full per-assembly pipeline, and `main()` is a thin dispatcher. Single-assembly mode (`-q`) and multi-assembly mode (`--fofn`/`--assembly-dir`) flow through the same code path — single-assembly simply produces a 1-element assembly list. Both use `--output-dir` (`-o`) with a uniform layout: reference files go under `{output_dir}/reference/` and each assembly gets its own `{output_dir}/{name}/` subfolder. In single-assembly mode, the assembly name is derived from `--assembly-name` or the query filename stem. If one assembly fails, others continue; a summary is logged at the end. Cross-assembly comparison outputs are produced when there are ≥2 assemblies.

**Gate-based assignment**: Contigs must satisfy ALL criteria for chromosome assignment (AND logic). This prevents spurious assignments from single conserved genes or repetitive sequences. Gate requirements differ by mode:
- **Protein mode**: Requires ≥5 segments, ≥3 genes, ≥50kb span, ≥20% span fraction (configurable via `--miniprot-min-*` flags)
- **Nucleotide mode**: Requires ≥1 segment (hardcoded), ≥50kb span, ≥20% span fraction. Lower segment threshold because perfect full-length alignments produce fewer segments than fragmented protein hits.

**Interval tree optimization**: `chain_parsing.py:_filter_overlapping_hits_by_identity()` uses the `intervaltree` package for O(n log n) overlap detection instead of O(n²). This is critical for large assemblies with many protein alignments.

**Reference ID normalization**: `reference_utils.py` handles polyploid naming conventions (chr1A/chr1B) and organelle aliases (chrC/Pt, chrM/Mt). Configurable via `set_ref_id_patterns()`.

**Depth caching**: `read_depth.py` caches alignment and depth results with metadata validation (query MD5, window size, etc.) to avoid re-running expensive alignment steps. Cache is automatically invalidated if inputs change.

**Logging infrastructure**: Custom formatter in `logging_config.py` provides `[info]`, `[warn]`, `[error]`, `[done]` prefixes. Use `logger.phase()` for major pipeline phases and `logger.done()` for completion messages.

**TOML configuration**: `config.py` provides schema-based config file support. CLI arguments override config file values.

**Distributed computing** (`--cluster`): Optional SLURM job submission via [executorlib](https://github.com/pyiron/executorlib). When `--cluster` is set, compute-intensive phases (synteny alignment, BLAST detection, debris detection, contaminant screening, read depth, compleasm) are submitted as individual SLURM jobs with per-job resource control. The coordinator process orchestrates submission and waits on futures. Two levels of parallelism: intra-assembly (independent phases like organelle + rDNA BLAST run as parallel SLURM jobs) and inter-assembly (multiple `run_assembly()` calls run concurrently via ThreadPoolExecutor, each submitting its own SLURM jobs). When `--cluster` is not set, `LocalExecutor` provides synchronous execution with zero overhead — behavior is identical to the non-distributed code path. If `--cluster` is set but executorlib (or its optional dependencies pysqa/h5py) is not installed, the tool exits with a clear error message and install instructions. Resource estimation (`resource_estimation.py`) sizes each job based on input file sizes and caps against `--max-threads-dist`, `--max-mem-dist`, `--max-time-dist`.

### Critical Security Notes

**Command injection prevention**: `alignment/external_tools.py` validates extra arguments against shell metacharacters using regex pattern `[;&|` `$(){}[\]<>!\\'"#]`. This prevents command injection via user-provided arguments like `--miniprot-extra-args`.

**TOCTOU race condition fix**: `utils/io_utils.py:file_exists_and_valid()` uses a single `stat()` call instead of separate `exists()` + `getsize()` to avoid time-of-check-time-of-use vulnerabilities.

**Process cleanup**: After killing subprocesses, code calls `kill()` then `wait()` to reap zombie processes and prevent resource leaks.

## Common Gotchas

**GC deviation returns None**: The GC deviation helpers `_gc_deviation_vs_ref()` and `_gc_deviation_vs_asm()` (inner functions of `classifier.py:classify_all_contigs()`) return `None` (not 0.0) when GC std < 0.001. This prevents division-by-zero and ensures low-confidence classification for contigs when GC comparison is meaningless.

**Subsample fraction validation**: `read_depth.py` validates `--depth-target-coverage` is non-negative. Invalid fractions raise `ValueError` with clear error messages.

**Organelle reference extraction**: If `--chrC-ref` or `--chrM-ref` not provided, code attempts to extract from main reference FASTA. If organelle sequences not found, organelle detection is skipped (not an error).

**Synteny mode selection**: The tool supports two synteny modes controlled by `--synteny-mode`:
- **nucleotide** (default): Uses minimap2 whole-genome nucleotide alignment with permissive chaining for structural composition analysis. Calls `parse_paf_chain_evidence_and_segments()`.
- **protein**: Uses miniprot protein-anchored synteny. Requires `--ref-gff3`. Calls `parse_miniprot_synteny_evidence_and_segments()`. Both parsers in `chain_parsing.py` produce compatible `ChainEvidenceResult` objects.

**Permissive synteny chaining in nucleotide mode**: Nucleotide mode always uses permissive minimap2 chaining parameters (`--max-chain-skip=300 -z 10000,1000 -r 50000`) to create megabase-scale synteny blocks by chaining through repetitive regions and small gaps (kb-scale). This approach:
- Works well for both within-species and cross-species comparisons
- Chains through homopolymer runs, tandem repeats, and ambiguous regions that would otherwise fragment alignments
- Produces clean chromosome-scale blocks suitable for compositional analysis and visualization
- Is appropriate for final_finalizer's goal (chromosome classification and architecture) rather than fine-scale SV detection

The permissive parameters are balanced by downstream filtering (identity thresholds, minimum alignment length, gate filtering) to prevent spurious assignments.

**Segment counts differ by mode**: In nucleotide mode, perfect full-length alignments often produce very few segments (even just 1!) because they're continuous matches. In protein mode, the same chromosome produces many segments (fragmented protein hits). The gate filtering accounts for this: nucleotide mode requires ≥1 segment (hardcoded), while protein mode requires ≥5 segments (configurable).

**Why permissive chaining**: Even for byte-for-byte identical sequences, minimap2's conservative default parameters split alignments at complex repetitive regions (long homopolymer runs, tandem repeats). For chromosome classification and structural composition visualization, these kb-scale gaps are noise that obscures megabase-scale patterns. Permissive chaining creates continuous blocks that better represent biological chromosome architecture.

**Reference length normalization**: `reference_utils.py:get_min_nuclear_chrom_length()` filters out organelles (chrC/chrM) and non-chromosome sequences via `is_nuclear_chromosome()` when computing the smallest nuclear chromosome length. This ensures `--chr-like-minlen` thresholds are computed correctly.

**rDNA consensus building**: The tool builds a species-specific 45S rDNA consensus from the query assembly by default and annotates rRNA sub-features (18S, 5.8S, 25S/28S, ITS1, ITS2). Use `--skip-rdna-consensus` to disable. Sub-feature annotation uses Infernal/cmscan with bundled Rfam 15.0 covariance models for structure-based boundary detection. The bundled Rfam database (`data/rfam/euk-rrna.cm`) contains 4 eukaryotic rRNA models (5S, 5.8S, 18S, 28S) and is automatically pressed (indexed) on first use. Requires Infernal (`conda install -c bioconda infernal`). Output includes a GFF3 file (`*.rdna_annotations.gff3`) with hierarchical features using proper Sequence Ontology terms (SO:0001637 for rRNA_gene, SO:0000252 for rRNA, SO:0000635 for ITS).

**Rfam database auto-pressing**: The bundled Rfam covariance models are stored in text format and automatically pressed to binary indices (`.i1f`, `.i1i`, `.i1m`, `.i1p`) by Infernal on first use. The tool checks for existing indices and only presses if they're missing or outdated. This eliminates the need for manual database preparation.

**Compleasm (BUSCO) evaluation**: Phase 17 runs compleasm on two FASTA subsets: `chrs.fasta` (chromosome-assigned contigs) and `non_chrs.fasta` (debris + unclassified + contaminants). Requires `--compleasm-lineage` (e.g., `eukaryota`, `viridiplantae`, `embryophyta`). Optionally specify `--compleasm-library` for pre-downloaded lineage files to avoid runtime downloads. Both compleasm runs are submitted in parallel via the executor. Results are stored as `CompleasmResult` dataclass fields on `AssemblyResult` (`compleasm_chrs`, `compleasm_non_chrs`) and included in `comparison_summary.tsv` columns (`compleasm_lineage`, `compleasm_chrs_S/D/F/I/M`, `compleasm_non_chrs_S/D/F/I/M`). The detection module is in `detection/compleasm.py`, which parses compleasm's `summary.txt` output format. Cached results are reused if the output directory already contains a `summary.txt`.

## Testing Philosophy

Tests are in `tests/` directory (7 test files, ~165 tests):
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
- infernal (cmscan) - structure-based rRNA annotation with Rfam models (rDNA consensus building; enabled by default, skip with --skip-rdna-consensus)
- compleasm - BUSCO completeness evaluation (requires `--compleasm-lineage`; skip with `--skip-compleasm`; **install in a separate conda environment** due to dependency conflicts — see README.md)
- Rscript (with ggplot2, dplyr, etc.) - visualization

**Python packages**:
- intervaltree - efficient overlap detection
- executorlib + pysqa + h5py (optional) - SLURM job submission for `--cluster` mode (`conda install -c conda-forge executorlib pysqa h5py`)

All external tools are called via subprocess with proper error handling. Use `utils/io_utils.py:have_exe()` to check availability before calling.

## Output Files Reference

**FASTA outputs**: `*.chrs.fasta`, `*.organelles.fasta`, `*.rdna.fasta`, `*.contaminants.fasta`, `*.debris.fasta`, `*.unclassified.fasta`, `*.non_chrs.fasta` (combined non-chromosome contigs, produced when compleasm is enabled)

**TSV outputs**:
- `*.contig_summary.tsv` - Per-contig classification with all evidence fields
- `*.evidence_summary.tsv` - Synteny evidence per contig × reference
- `*.segments.tsv` - Individual synteny segments from chain parsing
- `*.macro_blocks.tsv` - Aggregated synteny macro-blocks
- `*.ref_lengths.tsv` - Reference chromosome lengths
- `*.contaminants.tsv` - Detailed contaminant summary with taxonomic lineage (if contaminants detected)
- `*.rdna_annotations.tsv` - rRNA sub-feature annotations in TSV format (produced by default; skip with `--skip-rdna-consensus`)
- `*.rdna_arrays.tsv` - rDNA array locations per contig (produced by default when arrays are detected; skip with `--skip-rdna-consensus`)

**GFF3 outputs**:
- `*.rdna_annotations.gff3` - Hierarchical rRNA gene annotations with 18S, 5.8S, 25S, ITS1, ITS2 sub-features (produced by default; skip with `--skip-rdna-consensus`)

**Visualization** (enabled by default; skip with `--skip-plot`; requires rmarkdown + pandoc):
- `*.unified_report.html` - Self-contained HTML report with all plots (chromosome overview, classification bar, depth overview, contaminant table)
- Individual PDFs (`*.chromosome_overview.pdf`, etc.) are also exported from within the report

**Output directory structure** (uniform for single- and multi-assembly):
```
{output_dir}/
  reference/
    reference.ref_lengths.tsv
    reference.ref_proteins.fa           # protein mode
    reference.ref_gff3.filtered.gff3    # if GFF3 provided
  {assembly_name}/
    {assembly_name}.contig_summary.tsv
    {assembly_name}.chrs.fasta
    ...
  {assembly_name_2}/                    # multi-assembly only
    ...
  comparison_summary.tsv                # multi-assembly only (≥2 assemblies)
  chromosome_completeness.tsv           # multi-assembly only (≥2 assemblies)
```

See `docs/output_formats.md` for complete column documentation.

## Plot Design Guidelines

When creating or modifying visualizations:

- **Maximum width**: 7.2 inches for all PDF outputs (fits standard page margins)
- **DPI**: 300 for print-quality output
- **Device**: Use `cairo_pdf` for better font rendering
- **HTML versions**: Create interactive HTML versions using ggiraph where appropriate
- **Font**: Use a consistent base font family (e.g., "Helvetica" or system default)
- **Colors**: Use colorblind-friendly palettes where possible

## Modifying Classification Logic

Classification happens in `classification/classifier.py:classify_all_contigs()`. The function applies a decision tree:

Before the decision tree, reference assignment runs for contigs that have synteny evidence. In nucleotide mode, each contig is assigned to the reference chromosome with the highest combined span fraction: `max(ref_span_bp / ref_length, query_union_bp / query_length)`. Reference span fraction (`ref_span_bp / ref_length`) is size-normalized and correctly handles simple translocations between chromosomes of unequal size. Query span fraction (`query_union_bp / query_length`) handles the complementary case where a translocation occupies most of the contig but covers only a small fraction of a large reference chromosome. Taking the maximum of both metrics ensures correct assignment in either scenario. This combined metric is computed directly in `chain_parsing.py` from reference lengths extracted from the PAF file. In protein mode, the combined span fraction falls back to raw synteny score because miniprot PAF does not carry reference chromosome lengths.

Each contig is scored independently. Multiple contigs can still be assigned to the same reference (valid for polyploids).

1. Organelle complete (if BLAST coverage ≥80%)
2. rDNA (if BLAST coverage ≥50%)
3. Chromosome assigned (if passes ALL synteny gates)
4. Chromosome unassigned (if length ≥ `--chr-like-minlen` but failed synteny gates)
5. Contaminant (if centrifuger score ≥1000 AND coverage ≥0.50; two-gate filtering prevents false positives from conserved genes)
6. Chromosome debris (if high coverage/identity vs assembled chromosomes)
7. Organelle debris (if partial organelle match)
8. Debris (if reference coverage >50% or protein hits ≥2)
9. Unclassified (no evidence)

**Confidence levels** are assigned by `classifier.py:assign_classification_confidence()` based on category-specific criteria (gene proportion, GC deviation, coverage, identity, protein hits).

## When Modifying Synteny Block Building

The synteny block building algorithm in `chain_parsing.py` is critical for performance. Key optimizations:

- **Interval tree filtering**: Uses `intervaltree` package for O(n log n) overlap detection
- **Chain scoring**: Uses `assign_chain_topk=3` (configurable via `--assign-chain-topk`) to consider top-K chains per contig × reference, not just the best chain

If modifying chain parsing, ensure you handle both PAF (minimap2) and miniprot PAF formats correctly. Miniprot PAF includes gene IDs in alignment records.

## Configuration File Schema

`utils/config.py:CONFIG_SCHEMA` maps TOML sections to argparse argument names. To add new configurable parameters:

1. Add CLI argument in `cli.py:parse_args()`
2. Add mapping in `config.py:CONFIG_SCHEMA`
3. Update `--dump-config` template in `config.py:dump_config_template()`

CLI arguments always override config file values.
