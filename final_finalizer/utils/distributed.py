"""Distributed execution abstraction for final_finalizer.

Provides a ``concurrent.futures``-compatible interface that transparently
dispatches work to either:

* **LocalExecutor** — synchronous, zero-overhead execution (default)
* **executorlib SlurmClusterExecutor** — submits each call as a SLURM job
  with per-job resource control (requires ``--cluster`` flag)

When ``--cluster`` is set but executorlib is not installed the factory
:func:`create_executor` exits with a clear error message.
"""
from __future__ import annotations

import logging
import math
import random
import subprocess
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class ResourceSpec:
    """Resource requirements for a single distributed job."""

    cores: int = 1
    memory_gb: float = 2.0
    time_minutes: int = 60
    job_name: str = ""
    partition: str = ""
    qos: str = ""


@dataclass
class ClusterConfig:
    """Cluster-wide settings derived from CLI flags."""

    enabled: bool = False
    max_threads: int = 64
    max_mem_gb: float = 128.0
    max_time_minutes: int = 720
    partition: str = "cpuq"
    qos: str = ""


def clamp_resources(spec: ResourceSpec, config: ClusterConfig) -> ResourceSpec:
    """Cap *spec* values to cluster-wide maximums and fill defaults."""

    spec.cores = min(spec.cores, config.max_threads)
    spec.memory_gb = min(spec.memory_gb, config.max_mem_gb)
    spec.time_minutes = min(spec.time_minutes, config.max_time_minutes)
    if not spec.partition:
        spec.partition = config.partition
    if not spec.qos:
        spec.qos = config.qos
    return spec


# ---------------------------------------------------------------------------
# Local (synchronous) executor — zero overhead
# ---------------------------------------------------------------------------
class LocalFuture:
    """Synchronous future that executes immediately on construction."""

    def __init__(self, fn: Callable, *args: Any, **kwargs: Any) -> None:
        self._exception: Optional[BaseException] = None
        self._result: Any = None
        try:
            self._result = fn(*args, **kwargs)
        except BaseException as exc:
            self._exception = exc

    def result(self, timeout: float | None = None) -> Any:
        if self._exception is not None:
            raise self._exception
        return self._result

    def done(self) -> bool:
        return True


class LocalExecutor:
    """Synchronous executor using :class:`LocalFuture`.

    Implements the ``concurrent.futures.Executor`` interface (submit only).
    Used when ``--cluster`` is not set — identical behaviour to direct calls.
    """

    def submit(
        self,
        fn: Callable,
        /,
        *args: Any,
        resource_spec: ResourceSpec | None = None,
        **kwargs: Any,
    ) -> LocalFuture:
        # resource_spec is accepted but ignored in local mode
        return LocalFuture(fn, *args, **kwargs)

    def __enter__(self) -> LocalExecutor:
        return self

    def __exit__(self, *exc: Any) -> None:
        pass


# ---------------------------------------------------------------------------
# executorlib helpers
# ---------------------------------------------------------------------------
def _resource_spec_to_dict(spec: ResourceSpec) -> dict:
    """Convert a :class:`ResourceSpec` to the dict format executorlib/pysqa expects.

    executorlib uses ``cores`` to decide between MPI (``cache_parallel.py``)
    and single-process (``cache_serial.py``) execution.  Since our submitted
    functions are single-process wrappers around multi-threaded tools
    (minimap2, blastn, etc.), we always set ``cores=1`` to select the serial
    backend and use ``threads_per_core`` to request the actual CPU count.
    pysqa multiplies ``cores × threads_per_core`` for ``--cpus-per-task``.

    Keys passed through to the pysqa SLURM template:
    - ``cores × threads_per_core`` → ``#SBATCH --cpus-per-task``
    - ``memory_max`` → ``#SBATCH --mem`` (GB)
    - ``run_time_max`` → ``#SBATCH --time`` (seconds, converted to minutes)
    - ``partition`` → ``#SBATCH --partition``
    - ``qos`` → ``#SBATCH --qos`` (requires :func:`_patch_pysqa_template`)
    - ``job_name`` → ``#SBATCH --job-name``
    """
    d: dict = {
        "cores": 1,                              # 1 process → cache_serial.py (no MPI)
        "threads_per_core": spec.cores,           # pysqa: 1 × N → --cpus-per-task=N
        "memory_max": math.ceil(spec.memory_gb),  # pysqa renders --mem={{memory_max}}G; SLURM needs integer
        "run_time_max": spec.time_minutes * 60,   # pysqa expects seconds
    }
    if spec.partition:
        d["partition"] = spec.partition
    if spec.qos:
        d["qos"] = spec.qos
    if spec.job_name:
        d["job_name"] = spec.job_name
    return d


class _ExecutorlibWrapper:
    """Thin wrapper around executorlib's ``Executor`` so that callers can pass
    :class:`ResourceSpec` via a ``resource_spec`` keyword argument on
    ``submit()``.
    """

    def __init__(self, executor: Any) -> None:
        self._executor = executor

    def submit(
        self,
        fn: Callable,
        /,
        *args: Any,
        resource_spec: ResourceSpec | None = None,
        **kwargs: Any,
    ) -> Any:
        if resource_spec is not None:
            kwargs["resource_dict"] = _resource_spec_to_dict(resource_spec)
        return self._executor.submit(fn, *args, **kwargs)

    def __enter__(self) -> _ExecutorlibWrapper:
        self._executor.__enter__()
        return self

    def __exit__(self, *exc: Any) -> None:
        self._executor.__exit__(*exc)


# ---------------------------------------------------------------------------
# pysqa template patching
# ---------------------------------------------------------------------------
_QOS_DIRECTIVE = """\
{%- if qos %}
#SBATCH --qos={{qos}}
{%- endif %}"""


def _patch_pysqa_template() -> None:
    """Add ``--qos`` support to pysqa's built-in SLURM template.

    The default pysqa SLURM template does not include a ``--qos`` directive.
    This function injects a conditional ``#SBATCH --qos={{qos}}`` line so that
    the QoS value from :class:`ResourceSpec` flows through to ``sbatch``.

    Safe to call multiple times — patches only once.
    """
    try:
        from pysqa.wrapper import slurm as _slurm_mod  # type: ignore[import-untyped]
    except ImportError:
        return  # pysqa not installed; executor init will fail with its own error
    if "qos" not in _slurm_mod.template:
        _slurm_mod.template = _slurm_mod.template.replace(
            "#SBATCH --cpus-per-task={{cores}}",
            _QOS_DIRECTIVE + "\n#SBATCH --cpus-per-task={{cores}}",
        )


# ---------------------------------------------------------------------------
# sbatch retry patching
# ---------------------------------------------------------------------------
def _patch_sbatch_retry(max_retries: int = 8, base_delay: float = 15.0) -> None:
    """Add retry-with-backoff to executorlib's sbatch submission.

    When many SLURM jobs are submitted concurrently (e.g. 16 assemblies each
    submitting multiple phases), sbatch can fail transiently with exit code 1
    due to ``QOSMaxJobsPerUserLimit`` or similar rate limits.  executorlib's
    internal thread does not retry, so one failure kills the entire thread
    and hangs all pending futures.

    pysqa stores a direct reference to ``pysqa_execute_command`` at
    construction time (``self._execute_command_function``), so replacing the
    module attribute does not affect calls made through the stored reference.
    Instead, we replace the ``subprocess`` module in the scheduler's namespace
    with a copy whose ``check_output`` retries on ``CalledProcessError``.
    The original function always looks up ``subprocess`` from its own module
    globals, so the retry wrapper is invoked regardless of who holds the
    function reference.

    Safe to call multiple times — patches only once.
    """
    try:
        from executorlib.standalone import scheduler as _sched_mod  # type: ignore[import-untyped]
    except ImportError:
        return

    if getattr(_sched_mod, "_sbatch_retry_patched", False):
        return  # already patched

    import types as _types

    _real_subprocess = _sched_mod.subprocess
    _orig_check_output = _real_subprocess.check_output

    def _check_output_with_retry(*args: Any, **kwargs: Any) -> Any:
        last_exc: subprocess.CalledProcessError | None = None
        for attempt in range(max_retries):
            try:
                return _orig_check_output(*args, **kwargs)
            except subprocess.CalledProcessError as exc:
                last_exc = exc
                if attempt < max_retries - 1:
                    delay = base_delay * (2 ** attempt) + random.uniform(0, 5)
                    logger.warning(
                        f"sbatch failed (exit {exc.returncode}), "
                        f"retrying in {delay:.0f}s "
                        f"(attempt {attempt + 1}/{max_retries})"
                    )
                    time.sleep(delay)
        raise last_exc  # type: ignore[misc]

    # Build a drop-in replacement for the subprocess module with retry logic
    _retry_mod = _types.ModuleType("subprocess")
    _retry_mod.__dict__.update(_real_subprocess.__dict__)
    _retry_mod.check_output = _check_output_with_retry  # type: ignore[attr-defined]
    _sched_mod.subprocess = _retry_mod  # type: ignore[attr-defined]
    _sched_mod._sbatch_retry_patched = True  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
def create_executor(config: ClusterConfig) -> LocalExecutor | _ExecutorlibWrapper:
    """Return an executor matching *config*.

    * ``config.enabled == False`` → :class:`LocalExecutor`
    * ``config.enabled == True`` and executorlib installed →
      :class:`_ExecutorlibWrapper` around ``SlurmClusterExecutor``
    * ``config.enabled == True`` but executorlib missing →
      raises :class:`SystemExit` with install instructions
    """
    if not config.enabled:
        return LocalExecutor()

    missing: list[str] = []
    try:
        from executorlib import SlurmClusterExecutor  # type: ignore[import-untyped]
    except ImportError:
        missing.append("executorlib")
    try:
        import mpi4py  # noqa: F401  # type: ignore[import-untyped]
    except ImportError:
        missing.append("mpi4py")

    if missing:
        logger.error(
            f"--cluster requires packages that are not installed: {', '.join(missing)}\n"
            "Install with:  conda install -c conda-forge executorlib mpi4py\n"
            "Or run without --cluster for local execution."
        )
        raise SystemExit(1)

    _patch_pysqa_template()
    _patch_sbatch_retry()

    try:
        inner = SlurmClusterExecutor()
        logger.info(
            "Distributed mode: executorlib SlurmClusterExecutor "
            f"(partition={config.partition}, qos={config.qos})"
        )
        return _ExecutorlibWrapper(inner)
    except Exception as exc:
        logger.error(f"Failed to initialise SlurmClusterExecutor: {exc}")
        raise SystemExit(1) from exc
