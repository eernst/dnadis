# Implementation plan: read-depth & alignment evidence for fusion / misjoin discrimination

Status: **planned** (not yet implemented).
Governing TODO: `dev/docs/TODO.md` → "Read depth evidence for rearrangement confidence" (Detection section).

## 1. Motivation

In polyploid (subgenome) assemblies a chromosome copy can appear as one full-length
contig plus separate fragments. Sometimes a subgenomic copy is **fused** to a chunk of
a *different* chromosome in the primary assembly (a chimeric contig). The fusion leaves
the copy's other half stranded as a fragment that reference-guided scaffolding cannot
join, because the fused contig has no free terminus on the relevant side.

The scientific question at each such junction: **real biological fusion** vs **assembly
misjoin**. Two independent lines of evidence discriminate them:

- **Read-depth uniformity across the junction.** Uniform long-read depth spanning the junction
  supports a real fused chromosome; a depth drop or spike at the junction suggests a
  misassembly.
- **Alignment continuity across the junction.** Reads (or the reference alignment)
  spanning the junction continuously supports a real fusion; a pileup of clipped/split
  reads exactly at the breakpoint suggests a misjoin.

This feature detects contig-internal fusion junctions, quantifies depth uniformity and
alignment continuity across each, and reports a `biological_fusion` /
`assembly_artifact` / `inconclusive` verdict with confidence — kept **separate** from the
existing combined rearrangement confidence score, per the TODO.

Telomere geometry is a useful complementary prior (already surfaced in the synteny plot):
a fusion contig capped by telomeres on both ends favors "leave it," but cannot by itself
distinguish a genuine fused chromosome from a misjoin that welded two telomere-bearing
ends — which is exactly what this feature adds.

## 2. Existing machinery (verified against current code)

### 2.1 `dnadis/detection/rearrangements.py` — fusion detection & the junction gap
- `detect_rearrangements()` groups `macro_block_rows` by `contig → ref_id → [block dict]`,
  then per contig runs `_detect_translocations`, `_detect_inversions`, `_detect_fusions`,
  and cross-contig `_detect_fissions`.
- `_detect_fusions()` (line ~352): a contig is a fusion candidate when ≥2 references each
  cover ≥ `_MIN_FUSION_FRAC` (0.30) of the contig. It builds
  `significant_refs = [(ref_id, span, q_min, q_max), …]`, sorts by query position, and for
  each pair emits a `RearrangementCall(rearrangement_type="fusion", …)`.
  **Gap:** the call sets `ref_start=0, ref_end=0` and `query_start=min(q_min)`,
  `query_end=max(q_max)` — spanning the whole contig. **It never computes the junction
  coordinate.** The junction lives between the sorted pair's query intervals: around
  `qa_max` (5′ ref) ↔ `qb_min` (3′ ref). The per-ref `q_min/q_max` needed to locate it are
  already in scope inside the pair loop.
- `_unpack_row()` (line ~102) parses the 19-field macro-block tuple (17-field legacy also
  handled), exposing per-chain `qstart, qend, ref_start, ref_end, strand, identity,
  n_segments`. `_block_bounds()`/`_union_query_span()` give the merged spans.
- `RearrangementCall` (models.py, ~line 192) is the carrier; `_score_confidence()` scores
  from span + collinear/reciprocal flags only. We extend the call and adjust confidence.

### 2.2 `dnadis/analysis/read_depth.py` — depth computation
- `calculate_depth_metrics()` (~line 442+) is the phase-13 entry: format detection,
  optional downsample, minimap2 → sorted BAM, `run_mosdepth()`, `parse_mosdepth_regions()`.
  Returns `Dict[contig → DepthStats]`.
- `run_mosdepth()` (~line 799) runs `mosdepth --by <window_size> --no-per-base …`, producing
  `<work_dir>/read_depth/depth.regions.bed.gz` (windowed depth; **no per-base track**).
  Default `window_size=1000` (CLI `--depth-window-size`).
- `parse_mosdepth_regions()` (~line 858) reads the gzipped BED into `contig_windows`
  (per-contig lists of window depths) and **aggregates to per-contig `DepthStats`,
  discarding the per-window values.** There is currently **no** accessor to retrieve the
  depth profile across a coordinate range — the main new primitive (§3.2).
- BAM lifecycle: the sorted BAM is deleted after mosdepth unless a keep flag is set, so
  split-read/clipping analysis (needs the BAM) is opt-in; depth-based evidence works from
  the persisted `depth.regions.bed.gz` with **no BAM** (hard TODO requirement).
- Caching validates `window_size` (among other keys); reducing the default to 100 will
  correctly invalidate stale caches.

### 2.3 Pipeline phases in `cli.py`
- Phase 12 = rearrangement detection (`detect_rearrangements(...)` → `write_rearrangements_tsv`).
- Phase 13 = read depth (only when `reads and not args.skip_depth`), submitted on the
  executor (may be a SLURM future), then `DepthStats` attached to each `ContigClassification`.
- **Depth runs strictly after rearrangement detection today**, so the TODO's phase swap is
  required to make depth available when scoring fusions.

### 2.4 `dnadis/alignment/chain_parsing.py`
- Macro-block rows carry per-chain contig coords (`qstart/qend`) and reference coords
  (`ref_start/ref_end`); this data already flows into `detect_rearrangements`. No change
  needed to locate junctions — only to consume it.

### 2.5 `dnadis/output/scaffolding.py` — fragment withholding
- `scaffold_chromosomes()` groups contigs by `(assigned_ref_id, query_subgenome_grp)`. A
  fused contig keeps `assigned_ref_id` but may sit in a different subgenome group than its
  partner fragment, so they never co-scaffold; there is **no explicit "withheld" list** —
  such fragments simply emit as trivial single-contig scaffolds. We will **annotate/flag**
  the withheld fragment, not change grouping (reassignment is a separate TODO).

### 2.6 Models, outputs, config, report
- `RearrangementCall` (models.py), `DepthStats` (models.py), `ContigClassification` depth
  fields.
- `write_rearrangements_tsv()` (tsv_output.py) — append evidence columns.
- Report `assembly_report.tmpl.Rmd` §"Rearrangements / Misassemblies" — reads the
  rearrangements TSV; `has_depth` already gates depth-dependent content.
- Config: `utils/config.py:CONFIG_SCHEMA` maps toml sections → argparse dests;
  `dump_config_template()` auto-generates from it. Note pre-existing gap:
  `--rearrangement-threshold` / `--rearrangement-density-frac` are CLI-only and absent from
  `CONFIG_SCHEMA` — backfill them while adding the new section.

## 3. Design

### 3.1 Junction detection (new module `dnadis/analysis/junction_evidence.py`)
`find_fusion_junctions(macro_block_rows, best_ref, query_lengths, min_fusion_frac,
min_junction_flank_bp) -> list[FusionJunction]`:
- Reuse the fusion-candidate logic from `_detect_fusions` (promote `_unpack_row` →
  public `unpack_macro_block_row`, keep a private alias; import `_union_query_span`/
  `_block_bounds` or promote likewise — avoid duplicating the parser).
- For each contig, collect refs whose merged query span ≥ `min_fusion_frac × contig_len`.
  Sort by query midpoint; for each adjacent pair define the **junction interval**
  `[left.q_max, right.q_min]`; `junction_point` = its midpoint (handle overlap where
  `right.q_min < left.q_max` by taking the overlap midpoint).
- Emit `FusionJunction(contig, left_ref, right_ref, junction_start, junction_end,
  junction_point, contig_len, left_span, right_span)`.

### 3.2 Depth statistics across the junction
New accessor in `read_depth.py` (the current parser throws windows away):
`read_depth_windows(regions_bed_gz, contig, start, end) -> list[tuple[int,int,float]]`
— stream the gzipped BED (same format `parse_mosdepth_regions` handles) and return
`(win_start, win_end, depth)` for windows overlapping `[start,end)` on `contig`.
Allocation-light scan (no full-file dict) so it is cheap per junction. Returns `[]` when
the file is absent (no-`--reads` path) — never raises.

Then `score_depth_uniformity(...) -> DepthEvidence` in `junction_evidence.py`:
- Three zones around `junction_point`: **left flank** `[jp − flank − pad, jp − pad)`,
  **junction** `[jp − pad, jp + pad)`, **right flank** `[jp + pad, jp + flank + pad)`.
  Defaults: `flank = 5 kb` (TODO), `pad = 1 window`.
- Metrics: `depth_left_median`, `depth_right_median`, `depth_junction_median`;
  **flank ratio** `min(L,R)/max(L,R)` (≈1 = symmetric); **junction ratio**
  `depth_junction_median / mean(L, R)` (≈1 supports fusion; <drop or >spike supports
  artifact); **CV** of junction-zone windows vs flanks (elevated junction CV = discontinuity).
- Guard rails: require min flank median depth (`--fusion-min-evidence-depth`, default ~10×)
  and a min number of informative windows; else `verdict="insufficient_depth"` →
  inconclusive.
- Return `DepthEvidence(depth_flank_ratio, depth_junction_ratio, depth_junction_cv,
  depth_flank_cv, n_windows, verdict)`.

Graceful degradation: no `depth.regions.bed.gz` → depth fields `None`, verdict falls back
to alignment-only.

### 3.3 Alignment continuity across the junction
- **Tier A — reference-PAF continuity (always available):** from the macro-blocks, assess
  whether the left/right ref segments abut cleanly (junction interval ≈ 0, collinear ref
  coords each side) or leave a large unaligned contig gap at the switch. Needs no reads;
  reuses `ref_start/ref_end`.
- **Tier B — split-read / clipping pileup (opt-in, needs BAM):** only when the sorted BAM
  is retained (`--keep-depth-bam` or a new `--fusion-split-reads`). Add
  `read_depth.py:clipping_pileup(bam_path, contig, junction_start, junction_end, threads)`
  using `samtools view` over the junction; count reads spanning unclipped vs clipped/split
  (`SA`-tagged) at the breakpoint. Guard with `have_exe("samtools")`; validate any extra
  args per the `external_tools.py` shell-metacharacter convention. Absent BAM/samtools →
  `None`, degrade to Tier A + depth.

### 3.4 Scoring & classification
`classify_fusion_evidence(depth_ev, aln_ev, clip_ev, call, thresholds) -> FusionVerdict`:
- `fusion_call ∈ {biological_fusion, assembly_artifact, inconclusive}`
- `fusion_confidence ∈ {high, medium, low}`
- `evidence_codes`: subset of `{Aln, Depth, Split}` (exact TODO vocab) for which evidence
  was available and positive.

Rules (all thresholds configurable):
- **biological_fusion**: depth uniform (`depth_junction_ratio ∈ [drop, spike]`,
  `depth_flank_ratio ≥ ~0.7`, junction CV not elevated) AND (Tier A clean abutment OR
  Tier B spanning reads). Confidence `high` if depth+split agree; `medium` if depth-only +
  clean Tier A; `low` if a single weak signal.
- **assembly_artifact**: depth drop/spike at junction OR clip pileup (Tier B) OR large
  unaligned contig gap with asymmetric flank depth. Confidence scales with agreeing signals.
- **inconclusive**: otherwise, or insufficient depth and no BAM.

Wire back: a fusion `RearrangementCall`'s existing `confidence` is upgraded/downgraded by
the verdict (TODO: uniform → upgrade, anomaly → downgrade), while the *new* fields carry
the discriminative detail. Combined confidence and evidence codes stay separate.

### 3.5 CLI flags (+ config schema + dump-config)
New argument group "Fusion / misjoin evidence" in `cli.py:parse_args()`:
- `--fusion-junction-flank` (int, 5000) — flank bp each side.
- `--fusion-depth-drop-ratio` (float, 0.65) — junction/flank ratio below → artifact.
- `--fusion-depth-spike-ratio` (float, 1.5) — above → artifact.
- `--fusion-min-evidence-depth` (float, 10.0) — min flank median depth to trust depth.
- `--fusion-min-frac` (float, 0.30) — expose `_MIN_FUSION_FRAC`.
- `--fusion-split-reads` (flag) — enable Tier B (requires/forces retaining the BAM).
- `--skip-fusion-evidence` (flag) — opt out.
- Consider `--depth-window-size` default 1000 → 100 (TODO); flag as a behavior change.

Per CLAUDE.md config workflow: add to `parse_args`; add a `fusion_evidence` section to
`CONFIG_SCHEMA` (and backfill `rearrangement_threshold`/`rearrangement_density_frac`);
`dump_config_template()` auto-updates.

### 3.6 Pipeline placement & scaffolding interaction
- **Swap phases 12 and 13** (12 → read depth, 13 → rearrangement detection) so depth is
  available when scoring fusions. Renumber log strings and update the phase table in
  CLAUDE.md. Depth stays optional; without `--reads`, phase 12 logs "skipping" and fusion
  scoring degrades to alignment-only.
- After `detect_rearrangements(...)`, run
  `annotate_fusion_evidence(rearrangement_calls, macro_block_rows, best_ref,
  query_lengths, depth_regions_path, bam_path_or_None, thresholds)` from
  `junction_evidence.py`. It mutates each fusion `RearrangementCall` (populates new fields,
  adjusts confidence) and returns per-junction records.
- **Scaffolding annotation (not re-grouping):** flag any `chrom_fragment` whose
  `assigned_ref_id` matches a `biological_fusion` junction's `left_ref`/`right_ref` but
  which landed in a different subgenome group as a rescue candidate
  (`fusion_partner_withheld=True`), and `logger.info` naming the withheld fragment and the
  fusion contig. Actual rescue = the separate "subgenome-aware fragment reassignment" TODO.

### 3.7 Outputs
- **Extend `<prefix>.rearrangements.tsv`**: append `junction_start`, `junction_end`,
  `depth_flank_ratio`, `depth_junction_ratio`, `depth_junction_cv`, `fusion_call`,
  `fusion_confidence`, `evidence_codes` (populated for fusion rows; empty otherwise). Mirror
  on `RearrangementCall`.
- **New `<prefix>.fusion_junctions.tsv`** (only when junctions exist): one row per junction
  with full metrics (left/right ref, junction coords, per-zone medians/CVs, clip counts if
  Tier B, verdict, evidence codes). New `write_fusion_junctions_tsv()` in `tsv_output.py`.
  Document both in `docs/output_formats.md`.
- **HTML report**: add `evidence_codes` + `fusion_call`/`fusion_confidence` to the
  Rearrangements table, and an optional per-junction depth-profile panel
  (left-flank / junction / right-flank medians) gated on `has_rearrangements && has_depth`.
  Follow plot guidelines (≤7.2 in, cairo_pdf, colorblind palette).

## 4. Build sequence (dependency order)

1. **models.py** — extend `RearrangementCall` with the new optional (defaulted) fields; add
   `FusionJunction` and evidence dataclasses (`DepthEvidence`/`ClipEvidence`/`FusionVerdict`
   may stay internal to `junction_evidence.py`). `from __future__ import annotations`.
2. **analysis/read_depth.py** — add `read_depth_windows()` and, for Tier B,
   `clipping_pileup()`. Optionally change `--depth-window-size` default to 100.
3. **detection/rearrangements.py** — promote `_unpack_row` → public
   `unpack_macro_block_row` (private alias retained) for sharing.
4. **analysis/junction_evidence.py** (new) — `find_fusion_junctions`,
   `score_depth_uniformity`, `assess_alignment_continuity` (Tier A),
   `classify_fusion_evidence`, orchestrator `annotate_fusion_evidence`. Pure functions where
   possible.
5. **output/tsv_output.py** — extend `write_rearrangements_tsv`; add
   `write_fusion_junctions_tsv`.
6. **cli.py** — new argument group; swap phases 12/13 and renumber; call
   `annotate_fusion_evidence`; write `fusion_junctions.tsv` (add `fusion_junctions_tsv` to
   `AssemblyResult`); scaffolding rescue-candidate annotation/log.
7. **utils/config.py** — `fusion_evidence` section (+ backfill rearrangement flags).
8. **output/reports/assembly_report.tmpl.Rmd** — table columns/badges + optional depth panel,
   gated on `has_depth`.
9. **CLAUDE.md / docs/output_formats.md** — phase-table swap, new flags, new/extended TSV
   columns; mark the TODO item done/partial.

## 5. Testing (`tests/`, pytest)

New `tests/test_junction_evidence.py`:
- **Junction location** — synthetic `macro_block_rows` for a fused contig (ref A `[0,12Mb]`,
  ref B `[12Mb,18Mb]`); assert `junction_point ≈ 12Mb`, correct left/right ref; cover
  abutting, gapped, overlapping, and reverse-strand right-arm layouts.
- **Depth scoring** — synthetic window lists: flat → biological (ratio ≈ 1); 50% drop →
  artifact; 2× spike → artifact; too-few/too-shallow → inconclusive. Assert ratio/CV math
  and threshold edges.
- **`read_depth_windows`** — tiny gzipped BED via `tmp_path`; assert slice/boundary
  handling and empty-list on missing file.
- **Classification matrix** — `classify_fusion_evidence` over present/absent depth, Tier A,
  Tier B; assert `evidence_codes ⊆ {Aln,Depth,Split}` and up/downgrade logic.
- **Tier B parsing** — hand-built SAM string for `clipping_pileup` (avoid a hard `samtools`
  unit-test dep; mark the tool-dependent path `@pytest.mark.integration`).
- **TSV round-trip** — write/re-read both TSVs; header order + empty-field handling for
  non-fusion rows.
- **Config** — `fusion_evidence` appears in `dump_config_template`; CLI overrides config.

Run `pytest -m "not integration"` for the fast suite; Tier B under `-m integration`.
