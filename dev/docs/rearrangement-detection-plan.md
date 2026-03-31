# Rearrangement Detection and Reporting — Plan

## Overview

Add a rearrangement detection module that classifies structural variants
per contig × reference chromosome pair, producing a typed annotation table
with confidence levels. Results are written as a TSV per assembly and
summarized in both the assembly report and comparison report.

## Rearrangement Types

### High confidence (clear macro_block signatures)

| Type | Signature | Evidence source |
|------|-----------|-----------------|
| Reciprocal translocation | Contig has substantial blocks to two refs; partner contig shows the reciprocal pattern | macro_blocks: off-target blocks to a consistent partner ref from both contigs |
| Whole-arm translocation | One arm of a contig maps to a different ref than the rest | macro_blocks: contiguous off-target region at one end of the contig |
| Chromosome fusion | Single contig covers two reference chromosomes end-to-end | macro_blocks: two refs with complementary, collinear coverage; one contig assigned, other ref partially/fully covered |
| Large inversion | Contiguous block(s) on the same ref but opposite strand | macro_blocks/segments: strand flip within otherwise collinear on-target blocks |
| Fission | One reference chromosome represented by two or more query contigs with complementary coordinate ranges | Assignment: multiple contigs assigned to same ref; macro_blocks: non-overlapping ref coordinate ranges |

### Moderate confidence

| Type | Signature | Limitations |
|------|-----------|-------------|
| Unbalanced translocation | Off-target blocks from a single ref, no reciprocal pattern in the partner | Hard to distinguish from insertions of duplicated/transposed sequence |
| Robertsonian translocation | Fusion signature specifically at centromeric regions | Requires centromere annotation to distinguish from generic fusion |

### Low confidence / currently limited

| Type | Limitation |
|------|-----------|
| Small inversions (<50kb) | Below macro_block min_span threshold; would need segment-level analysis |
| Insertions/deletions | Not represented as synteny blocks; detectable only as coverage gaps or unmapped query regions |
| Tandem duplications | Multiple blocks to same ref region; confused with assembly artifacts |

## Detection Algorithm

### Input
- `macro_block_rows` from chain parsing (per-assembly, ref→query)
- `best_ref` assignments
- `ref_lengths` for coordinate normalization
- `query_lengths` for contig size context

### Per-contig analysis

For each `chrom_assigned` contig:

1. **Collect blocks by reference**: Group macro_blocks by target ref_id.
   Identify on-target (matches assigned ref) vs off-target blocks.

2. **Translocation detection**:
   - For each off-target ref with total span ≥ threshold (e.g., 500kb or 5% of contig):
     - Record as translocation candidate with ref origin, query coordinate range, and span
     - Check if the off-target region is contiguous (whole-arm) or interstitial
     - Check for reciprocal pattern: does the partner contig (assigned to the off-target ref) have reciprocal off-target blocks?
   - Classify: reciprocal translocation, whole-arm translocation, or unbalanced translocation

3. **Inversion detection**:
   - Within on-target blocks, identify strand changes
   - Group consecutive same-strand blocks into segments
   - Strand-flip boundaries between segments indicate inversions
   - Record inversion coordinates, size, and whether it's terminal or interstitial

4. **Fusion detection**:
   - Contig maps substantially (>30% span) to two or more reference chromosomes
   - The coverage is collinear and complementary (ref A on one end, ref B on the other)
   - Distinguished from translocation by the scale (fusion = most of both refs represented)

5. **Fission detection**:
   - Two or more contigs assigned to the same ref with non-overlapping ref coordinate ranges
   - Combined coverage approaches full ref length
   - Neither contig alone covers the full ref

### Confidence scoring

Each call gets a confidence level based on:
- **Size**: Larger rearrangements are higher confidence (less likely to be noise)
- **Reciprocity**: Reciprocal translocations (evidence from both partners) are higher confidence than one-sided
- **Collinearity**: Rearrangements with clear collinear blocks on both sides of the breakpoint are higher confidence
- **Coverage**: Higher alignment coverage in the rearranged region increases confidence
- **Consistency**: Same rearrangement observed in multiple assemblies (comparison report) increases confidence

Proposed levels:
- **high**: Large (>1 Mb), reciprocal or collinear, good coverage
- **medium**: Moderate size (50kb-1 Mb), or large but one-sided
- **low**: Near the macro_block threshold, or ambiguous pattern

No additional threshold flags.  The detection inherits the pipeline's
existing macro_block minimum span (controlled by `--min-span-bp`,
default 50kb).  Any macro_block that survived chain parsing is already
credible alignment evidence.  Confidence tiers describe trust in the
biological interpretation, not the alignment quality.

## Output

### Per-assembly TSV: `*.rearrangements.tsv`

| Column | Type | Description |
|--------|------|-------------|
| contig | string | Query contig name (renamed) |
| original_name | string | Original contig name |
| assigned_ref_id | string | Contig's assigned reference chromosome |
| rearrangement_type | string | `reciprocal_translocation`, `translocation`, `whole_arm_translocation`, `inversion`, `fusion`, `fission` |
| partner_ref_id | string | Partner reference chromosome (for translocations/fusions) |
| query_start | int | Start position on query contig |
| query_end | int | End position on query contig |
| ref_start | int | Start position on reference (for the rearranged region) |
| ref_end | int | End position on reference |
| span_bp | int | Size of the rearranged region |
| strand | string | `+` or `-` (for inversions) |
| confidence | string | `high`, `medium`, `low` |
| evidence | string | Brief description of the evidence supporting the call |

### Assembly report integration

Add a "Rearrangements" tab or section showing:
- Per-chromosome rearrangement summary table
- Count of each type with confidence breakdown
- Visual indicators in the chromosome overview plot (already partially done with red off-target blocks), need a legend for this plot

### Comparison report integration

Add a rearrangement summary to the comparison report:
- Per-chromosome × per-assembly matrix showing rearrangement types detected
- Highlight rearrangements shared across multiple assemblies vs assembly-specific
- Could use icons or color-coded cells (e.g., arrows for translocations, U-turn for inversions)

## Implementation

### New module: `detection/rearrangements.py`

Functions:
- `detect_rearrangements(macro_block_rows, best_ref, ref_lengths, query_lengths)` → list of RearrangementCall dataclass instances (inherits macro_block min_span from pipeline)
- `classify_translocation(contig_blocks, assigned_ref, partner_ref)` → type + confidence
- `detect_inversions(on_target_blocks)` → list of inversion calls
- `detect_fusions(contig_blocks, ref_lengths)` → fusion call or None
- `detect_fissions(ref_contigs, macro_blocks, ref_lengths)` → list of fission calls

### New dataclass: `RearrangementCall` in `models.py`

### TSV writer: `write_rearrangements_tsv()` in `tsv_output.py`

### Pipeline integration

- Run after Phase 11 (classification) since it needs `best_ref` assignments
- Could be Phase 11b or a new phase number - use a new phase number, we shouldn't have sub-phases (a, b, etc.)
- Results stored on `AssemblyResult` for comparison report access

## Dependencies

No new external tools required. All detection is from existing macro_block
and segment data. The module is pure Python analysis of alignment evidence.

## Testing

- Synthetic test cases with known rearrangements (from existing TODO)
- Compare detected rearrangements against known karyotype differences in the development test dataset
- Edge cases: very small rearrangements near threshold, overlapping rearrangements, assembly artifacts that mimic rearrangements
