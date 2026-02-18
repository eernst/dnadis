# Distributed Computing Support via executorlib

## Context

final_finalizer processes genome assemblies sequentially, both across assemblies (multi-assembly mode) and within each assembly (tool phases). On SLURM clusters, this leaves significant capacity unused. This plan adds optional distributed execution via [executorlib](https://github.com/pyiron/executorlib), which provides a `concurrent.futures`-compatible interface that transparently submits Python functions as SLURM jobs with per-job resource control. When `--cluster` is not set, behavior is identical to current code.

## Design Overview

**Coordinator pattern**: The main process orchestrates job submission and result collection. Compute-intensive phases (synteny alignment, BLAST detection, read depth) are submitted as individual SLURM jobs. The coordinator is lightweight — mostly waiting on futures.

**Key insight**: The existing detection functions (`detect_organelles()`, `detect_rdna_contigs()`, etc.) already have fully serializable signatures (Path, dict, set, primitives) and are self-contained top-level module functions. They can be submitted directly to executorlib with **no wrapper task modules needed**.

**Two levels of parallelism**:
1. **Inter-assembly**: In multi-assembly mode, multiple `run_assembly()` calls run concurrently (coordinator threads), each submitting SLURM jobs
2. **Intra-assembly**: Within each assembly, independent phases run in parallel (e.g., organelle + rDNA BLAST)

## Design Decisions

- **Time units**: Minutes throughout (CLI, ResourceSpec, ClusterConfig). Convert to seconds only at the executorlib boundary in `_resource_spec_to_dict()`.
- **Synteny task scope**: SLURM job runs only the alignment tool (run_miniprot/run_minimap2_synteny). PAF parsing is done locally by the coordinator (takes seconds vs minutes/hours for alignment).
- **Phase 8 (debris/unclassified)**: Stays local. It's typically fast and depends on many prior results.

## New Files

### 1. `final_finalizer/utils/distributed.py` (~150 lines)

Core executor abstraction:

- **`ResourceSpec`** dataclass: `cores`, `memory_gb`, `time_minutes`, `job_name`, `partition`, `qos`
- **`ClusterConfig`** dataclass: `enabled`, `max_threads`, `max_mem_gb`, `max_time_minutes`, `partition`, `qos`
- **`clamp_resources(spec, config)`**: Cap a ResourceSpec to cluster-wide maximums, apply default partition/qos
- **`LocalFuture`**: Synchronous Future-like object — executes function immediately on `submit()`, provides `.result()`/`.done()` interface
- **`LocalExecutor`**: Synchronous executor using `LocalFuture`. Context manager. Used when `--cluster` is not set — zero overhead, identical behavior to current code
- **`create_executor(config)`**: Factory function. Returns `LocalExecutor` if cluster disabled. Returns `SlurmClusterExecutor` if executorlib available. Falls back to `LocalExecutor` with warning if executorlib not installed

### 2. `final_finalizer/utils/resource_estimation.py` (~200 lines)

Functions that estimate ResourceSpec per phase based on input sizes:

| Function | Phase | Memory Formula | Time Formula | Cores |
|----------|-------|---------------|-------------|-------|
| `estimate_synteny_resources()` | 2 | minimap2: 8 bytes/ref_base + 2GB; miniprot: 6 bytes/qry_base + 1GB | Scales with ref_size × qry_size | min(16, max_threads) |
| `estimate_blast_resources()` | 4, 5 | 2GB (constant for small DBs) | Scales with query size | min(8, max_threads) |
| `estimate_debris_resources()` | 6 | 8 bytes/qry_base + 1GB | Scales with query size | min(8, max_threads) |
| `estimate_contaminant_resources()` | 7 | Based on index file size × 1.2 + 2GB | 30 min (typical) | min(16, max_threads) |
| `estimate_depth_resources()` | 11.5 | 8 bytes/qry_base + 4GB sort buffer | Scales with reads file size | min(16, max_threads) |

Helper `_estimate_genome_bp_from_filesize(path)` estimates genome size from FASTA file size (accounts for gzip).

All estimates are capped by `clamp_resources()` against `--max-threads-dist`, `--max-mem-dist`, `--max-time-dist`.

## Modified Files

### 3. `final_finalizer/cli.py` — Pipeline restructuring

**`parse_args()`** — Add new argument group:

```
Distributed computing (SLURM cluster):
  --cluster                    Enable distributed execution via executorlib
  --max-threads-dist INT       Max threads per distributed job [64]
  --max-mem-dist FLOAT         Max memory (GB) per distributed job [128]
  --max-time-dist INT          Max wall time (minutes) per distributed job [720]
  --partition STR              SLURM partition for distributed jobs [cpuq]
  --qos STR                    SLURM QOS for distributed jobs [default]
```

**`run_assembly()`** — Add optional `executor` and `cluster_config` parameters:

```python
def run_assembly(args, ref_ctx, qry, assembly_name, reads, outprefix,
                 executor=None, cluster_config=None) -> AssemblyResult:
```

Restructured phase dispatch (cluster mode):

```
Phase 1 (local): Read query FASTA, GC content
    │
Phase 2 (SUBMIT): Synteny alignment (run_miniprot or run_minimap2_synteny)
    │ ← future.result() blocks
Phase 2b (local): Parse PAF, gate filtering, identify chromosome contigs
    │
    ├──→ Phase 4 (SUBMIT): detect_organelles()     ─┐
    └──→ Phase 5 (SUBMIT): detect_rdna_contigs()   ─┤ parallel
         (exclude_contigs=set() for Phase 5)         │
         ↓                                           │
    Phase 4-5 results collected ←────────────────────┘
    Reconcile: apply organelle exclusion to rDNA results locally
    │
Phase 6 (SUBMIT): detect_chromosome_debris()
    │ ← future.result()
Phase 7 (SUBMIT): detect_contaminants()
    │ ← future.result()
Phase 8-11 (local): Classification, orientation, naming
    │
Phase 11.5 (SUBMIT): calculate_depth_metrics()
    │ ← future.result()
Phase 12-13 (local): Output generation, plots
```

**Thread semantics**: When `--cluster` is active, each submitted function receives `spec.cores` as its thread count (from resource estimation). The `--threads` parameter only applies to local phases. When `--cluster` is not active, `--threads` applies everywhere (unchanged).

**Phases 4-5 parallel trick**: Currently Phase 5 (rDNA) depends on Phase 4 (organelle) via `exclude_contigs`. In cluster mode, submit both with `exclude_contigs=set()`, then reconcile locally: remove any organelle-classified contigs from rDNA results. The classifier already prioritizes organelle over rDNA, so this is safe. Avoids sequential SLURM job dependency.

**`main()`** — Multi-assembly dispatch with cluster mode:

```python
with create_executor(cluster_config) as executor:
    if cluster_config.enabled:
        # Run assemblies concurrently via ThreadPoolExecutor
        # Each run_assembly() is mostly I/O-bound (waiting on SLURM futures)
        with ThreadPoolExecutor(max_workers=len(assemblies)) as pool:
            futures = [pool.submit(run_assembly, ..., executor, cluster_config)
                       for asm in assemblies]
            for name, f in zip(names, futures):
                try: results.append(f.result())
                except Exception as e: failures.append(...)
    else:
        # Sequential loop (unchanged)
        for asm in assemblies:
            result = run_assembly(...)
```

### 4. `final_finalizer/utils/config.py`

Add to `CONFIG_SCHEMA`:
```python
"distributed": {
    "cluster": "cluster",
    "max_threads_dist": "max_threads_dist",
    "max_mem_dist": "max_mem_dist",
    "max_time_dist": "max_time_dist",
    "partition": "partition",
    "qos": "qos",
},
```

Update `dump_config_template()` to include the new section.

### 5. `CLAUDE.md` and `README.md`

- Document `--cluster` flag, new parameters, executorlib as optional dependency
- Add executorlib to the optional conda dependencies section

## Implementation Order

1. **`distributed.py`**: ResourceSpec, ClusterConfig, LocalFuture, LocalExecutor, create_executor()
2. **`resource_estimation.py`**: All estimation functions
3. **CLI arguments**: Add `--cluster` group to parse_args(), update config.py schema
4. **Phase 2** (synteny): Modify run_assembly() to submit alignment via executor. Test with LocalExecutor (zero behavioral change)
5. **Phases 4-5** (organelle + rDNA): Submit in parallel, reconcile locally
6. **Phase 6** (debris): Submit via executor
7. **Phase 7** (contaminant): Submit via executor
8. **Phase 11.5** (depth): Submit via executor
9. **Multi-assembly parallelism**: ThreadPoolExecutor in main()
10. **Tests**: Unit tests for resource estimation, integration test with LocalExecutor verifying identical output, test graceful fallback when executorlib not installed
11. **Documentation**: CLAUDE.md, README.md

## Verification

1. **Baseline**: Run unit tests (`pytest -m "not integration"`) — all must pass unchanged
2. **LocalExecutor parity**: Run with `--cluster` but without executorlib installed → falls back to LocalExecutor → identical output to non-cluster run
3. **Resource estimation**: Unit tests for each estimation function with known file sizes
4. **SLURM integration** (manual): Run on actual SLURM cluster with executorlib installed, verify jobs appear in `squeue`, correct resource requests, results match local run
