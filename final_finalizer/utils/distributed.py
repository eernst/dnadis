"""Distributed execution abstraction for final_finalizer.

Provides a ``concurrent.futures``-compatible interface that transparently
dispatches work to either:

* **LocalExecutor** — synchronous, zero-overhead execution (default)
* **executorlib SlurmClusterExecutor** — submits each call as a SLURM job
  with per-job resource control (requires ``--cluster`` flag)

When ``--cluster`` is set but executorlib is not installed the factory
:func:`create_executor` falls back to *LocalExecutor* with a warning so
that the pipeline always runs.
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
    """Convert a :class:`ResourceSpec` to the dict format executorlib expects."""
    return {
        "cores": spec.cores,
        "memory": spec.memory_gb,           # executorlib expects GB
        "run_time": spec.time_minutes * 60, # executorlib expects seconds
        "partition": spec.partition,
        "qos": spec.qos,
    }


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

    try:
        from executorlib import Executor  # type: ignore[import-untyped]

        inner = Executor(backend="slurm_submission")
        logger.info(
            "Distributed mode: executorlib SlurmClusterExecutor "
            f"(partition={config.partition}, qos={config.qos})"
        )
        return _ExecutorlibWrapper(inner)
    except ImportError:
        logger.error(
            "--cluster requires the executorlib package, which is not installed.\n"
            "Install with:  conda install -c conda-forge executorlib\n"
            "Or run without --cluster for local execution."
        )
        raise SystemExit(1)
