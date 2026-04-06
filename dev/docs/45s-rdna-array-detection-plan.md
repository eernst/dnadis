# 45S rDNA Array Detection and Annotation

## Context

The pipeline detects individual 45S rDNA repeat copies (`RdnaLocus` objects) by BLASTing a consensus 45S sequence against all assembly contigs. A basic `_mark_nor_candidates()` function sets a boolean `is_nor_candidate` flag on loci that cluster (>=3 copies within 50kb gaps), but there's no array-level data structure, no array-level GFF3 features, and no array summary in the report.

**Goal**: Detect tandem 45S rDNA repeat arrays as first-class objects, annotate them in GFF3 as `repeat_region` parent features, produce an array summary TSV, and display array-level summaries in the HTML report. Rename "NOR candidate" terminology to "45S rDNA array" throughout (NOR implies functional data we don't have).

## Changes

### 1. Add `RdnaArray` dataclass, modify `RdnaLocus`

**File**: `dnadis/models.py`

Add `RdnaArray` dataclass after `RdnaLocus` (~line 309):
```python
@dataclass
class RdnaArray:
    array_id: str                # "array_1", "array_2", ...
    contig: str
    start: int                   # 0-based, first locus start
    end: int                     # 0-based exclusive, last locus end
    strand: str                  # majority vote of constituent loci
    loci: List[RdnaLocus]        # constituent loci, sorted by start
    n_total: int
    n_full: int                  # copy_type == "full"
    n_partial: int
    n_fragment: int
    identity_median: float
    identity_min: float
    identity_max: float

    @property
    def span(self) -> int:
        return self.end - self.start
```

Modify `RdnaLocus`:
- Replace `is_nor_candidate: bool` -> `array_id: Optional[str] = None`
- Add backward-compat property:
```python
@property
def is_nor_candidate(self) -> bool:
    return self.array_id is not None
```

### 2. Replace `_mark_nor_candidates` with `_detect_arrays`

**File**: `dnadis/detection/rdna_consensus.py`

Replace `_mark_nor_candidates()` with:

- `_detect_arrays(loci, min_tandem_copies, max_tandem_gap) -> List[RdnaArray]`: Same clustering logic (group by contig, sort by start, cluster by gap <= max_tandem_gap, filter by >= min_tandem_copies), but builds `RdnaArray` objects and sets `locus.array_id` on each constituent locus.
- `_build_array(array_id, contig, cluster) -> RdnaArray`: Helper that constructs an `RdnaArray` from a cluster of loci, computing copy type counts, identity stats (via `statistics.median`), and strand majority vote.

Array IDs: sequential `array_1`, `array_2`, ... numbered by contig sort order then start position.

Update `annotate_contigs_with_consensus()`:
- Return type: `Tuple[List[RdnaLocus], List[RdnaArray]]`
- Call `_detect_arrays()` instead of `_mark_nor_candidates()`
- Remove the `is_nor_candidate=False` default in locus construction (field no longer exists; `array_id` defaults to `None`)

Update `build_rdna_consensus()`:
- Return type: `Tuple[Optional[RdnaConsensus], List[RdnaLocus], List[RdnaArray]]`
- All early returns become 3-tuples: `return None, [], []`

### 3. Update TSV and GFF3 output

**File**: `dnadis/output/tsv_output.py`

**Per-locus TSV** (`write_rdna_annotations_tsv`):
- Replace `is_nor_candidate` column with `array_id` (empty string if None)

**New array summary TSV** (`write_rdna_arrays_tsv`):
- Columns: `array_id, contig, start, end, span_kb, strand, n_total, n_full, n_partial, n_fragment, identity_median, identity_min, identity_max, contig_classification`
- Output file: `*.rdna_arrays.tsv`

**GFF3** (`write_rdna_annotations_gff3`):
- Add `arrays: Optional[List[RdnaArray]] = None` parameter
- For each array, emit a `repeat_region` (SO:0000657) parent feature with attributes: `ID=rdna_{array_id}`, `Name=45S_rDNA_array`, `repeat_type=rDNA_45S`, `n_copies`, `n_full`, `identity_median`
- For loci in arrays, add `Parent=rdna_{array_id}` to the rRNA_gene attributes
- Replace `is_nor_candidate` attribute with `array_id` (only when set)

### 4. Wire up CLI and plotting

**File**: `dnadis/cli.py`
- Unpack 3-tuple from `build_rdna_consensus()`
- Write `*.rdna_arrays.tsv` via `write_rdna_arrays_tsv()`
- Pass `arrays` to `write_rdna_annotations_gff3()`
- Pass `rdna_arrays_tsv` path to `run_unified_report()` and `run_plot()`

**File**: `dnadis/output/plotting.py`
- Add `rdna_arrays_tsv: Optional[Path] = None` parameter to `run_unified_report()` and `run_plot()`
- Pass through as `__RDNA_ARRAYS__` placeholder

### 5. Update report template

**File**: `dnadis/output/unified_report.tmpl.Rmd`

- Add `rdna_arrays_tsv: "__RDNA_ARRAYS__"` to YAML params
- Load arrays TSV in data setup chunk
- Replace "NOR Candidates" tab with "45S rDNA Arrays" tab showing array-level summary:
  - Columns: array_id, contig, span (kb), copy count breakdown, strand, identity median, identity range
- Update "Other Annotations" tab to filter by `array_id` instead of `is_nor_candidate`
- Update `rdna_annot` col_types to use `array_id` instead of `is_nor_candidate`

**File**: `dnadis/output/chromosome_overview.tmpl.R` (line 156)
- Update col_types: `is_nor_candidate` -> `array_id`

### 6. Add tests

**File**: `tests/test_dnadis.py`

- `test_detect_arrays_basic`: 5 loci on same contig within gap -> 1 array, correct metrics
- `test_detect_arrays_multiple_on_contig`: 2 clusters separated by large gap -> 2 arrays
- `test_detect_arrays_below_min_copies`: 2 loci -> no arrays, `array_id` stays None
- `test_detect_arrays_sets_array_id`: Verify constituent loci get `array_id` set
- `test_rdna_locus_is_nor_candidate_compat`: Property returns True when `array_id` set
- `test_detect_arrays_strand_majority`: Mixed strands -> array strand is majority

## Implementation order

1. `models.py` -- Add `RdnaArray`, modify `RdnaLocus`
2. `rdna_consensus.py` -- Replace `_mark_nor_candidates` -> `_detect_arrays`
3. `tsv_output.py` -- Update writers, add `write_rdna_arrays_tsv`
4. `plotting.py` -- Add `rdna_arrays_tsv` parameter
5. `cli.py` -- Wire everything together
6. `unified_report.tmpl.Rmd` + `chromosome_overview.tmpl.R` -- Update templates
7. `tests/test_dnadis.py` -- Add tests

## Verification

1. `conda activate dnadis && pytest tests/test_dnadis.py -v` -- all tests pass
2. Render on a dataset with `--build-rdna-consensus --plot --plot-html`:
   - Verify "45S rDNA Arrays" tab appears with array summary
   - Verify `*.rdna_arrays.tsv` has correct columns and array data
   - Verify GFF3 has `repeat_region` parent features wrapping constituent `rRNA_gene` features
   - Verify per-locus TSV has `array_id` column instead of `is_nor_candidate`
3. Render without `--build-rdna-consensus`: verify rDNA section shows stub gracefully
