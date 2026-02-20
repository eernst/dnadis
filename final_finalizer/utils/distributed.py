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
    qos: str = "default"


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

    Keys are passed through to the pysqa SLURM template
    (``pysqa.wrapper.slurm.template``):
    - ``cores`` → ``#SBATCH --cpus-per-task``
    - ``memory_max`` → ``#SBATCH --mem`` (GB)
    - ``run_time_max`` → ``#SBATCH --time`` (seconds, converted to minutes)
    - ``partition`` → ``#SBATCH --partition``
    - ``qos`` → ``#SBATCH --qos`` (requires :func:`_patch_pysqa_template`)
    - ``job_name`` → ``#SBATCH --job-name``
    """
    d: dict = {
        "cores": spec.cores,
        "memory_max": spec.memory_gb,            # pysqa expects GB
        "run_time_max": spec.time_minutes * 60,  # pysqa expects seconds
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
