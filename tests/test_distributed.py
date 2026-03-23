"""Unit tests for distributed computing support."""
from __future__ import annotations

import importlib.util
import sys
import textwrap
import types
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from final_finalizer.utils.distributed import (
    ClusterConfig,
    LocalExecutor,
    LocalFuture,
    ResourceSpec,
    _ExecutorlibWrapper,
    _patch_pysqa_template,
    _patch_sbatch_retry,
    _resource_spec_to_dict,
    clamp_resources,
    create_executor,
)
from final_finalizer.utils.resource_estimation import (
    _estimate_genome_bp_from_filesize,
    _file_size_bytes,
    _scale_time,
    estimate_blast_resources,
    estimate_contaminant_resources,
    estimate_debris_resources,
    estimate_depth_resources,
    estimate_pairwise_resources,
    estimate_synteny_resources,
)


# ===================================================================
# distributed.py tests
# ===================================================================
class TestResourceSpec:
    def test_defaults(self):
        spec = ResourceSpec()
        assert spec.cores == 1
        assert spec.memory_gb == 2.0
        assert spec.time_minutes == 60
        assert spec.job_name == ""
        assert spec.partition == ""
        assert spec.qos == ""

    def test_custom(self):
        spec = ResourceSpec(cores=16, memory_gb=64.0, time_minutes=360, job_name="test")
        assert spec.cores == 16
        assert spec.memory_gb == 64.0
        assert spec.time_minutes == 360
        assert spec.job_name == "test"


class TestClusterConfig:
    def test_defaults(self):
        cfg = ClusterConfig()
        assert cfg.enabled is False
        assert cfg.max_threads == 64
        assert cfg.max_mem_gb == 128.0
        assert cfg.max_time_minutes == 720
        assert cfg.partition == "cpuq"
        assert cfg.qos == ""


class TestClampResources:
    def test_caps_to_maximums(self):
        spec = ResourceSpec(cores=128, memory_gb=500.0, time_minutes=2000)
        cfg = ClusterConfig(max_threads=32, max_mem_gb=64.0, max_time_minutes=720)
        result = clamp_resources(spec, cfg)
        assert result.cores == 32
        assert result.memory_gb == 64.0
        assert result.time_minutes == 720

    def test_fills_default_partition_qos(self):
        spec = ResourceSpec(partition="", qos="")
        cfg = ClusterConfig(partition="gpuq", qos="highpri")
        result = clamp_resources(spec, cfg)
        assert result.partition == "gpuq"
        assert result.qos == "highpri"

    def test_preserves_explicit_partition_qos(self):
        spec = ResourceSpec(partition="custom", qos="custom")
        cfg = ClusterConfig(partition="gpuq", qos="highpri")
        result = clamp_resources(spec, cfg)
        assert result.partition == "custom"
        assert result.qos == "custom"

    def test_no_clamping_when_under_limits(self):
        spec = ResourceSpec(cores=4, memory_gb=8.0, time_minutes=30)
        cfg = ClusterConfig(max_threads=64, max_mem_gb=128.0, max_time_minutes=720)
        result = clamp_resources(spec, cfg)
        assert result.cores == 4
        assert result.memory_gb == 8.0
        assert result.time_minutes == 30


class TestLocalFuture:
    def test_success(self):
        f = LocalFuture(lambda x: x * 2, 21)
        assert f.done()
        assert f.result() == 42

    def test_exception(self):
        def bad():
            raise ValueError("boom")
        f = LocalFuture(bad)
        assert f.done()
        with pytest.raises(ValueError, match="boom"):
            f.result()

    def test_kwargs(self):
        def add(a, b=0):
            return a + b
        f = LocalFuture(add, 3, b=7)
        assert f.result() == 10


class TestLocalExecutor:
    def test_submit_and_result(self):
        with LocalExecutor() as ex:
            f = ex.submit(sum, [1, 2, 3])
            assert f.result() == 6

    def test_submit_ignores_resource_spec(self):
        spec = ResourceSpec(cores=99)
        with LocalExecutor() as ex:
            f = ex.submit(lambda: 42, resource_spec=spec)
            assert f.result() == 42

    def test_context_manager(self):
        ex = LocalExecutor()
        with ex as ctx:
            assert ctx is ex


class TestResourceSpecToDict:
    def test_conversion(self):
        spec = ResourceSpec(
            cores=8, memory_gb=16.5, time_minutes=30,
            partition="cpuq", qos="default",
        )
        d = _resource_spec_to_dict(spec)
        assert d["cores"] == 1  # always 1 → serial backend (no MPI)
        assert d["threads_per_core"] == 8  # pysqa: 1×8 → --cpus-per-task=8
        assert d["memory_max"] == 17  # ceil(16.5) → integer GB for SLURM
        assert d["run_time_max"] == 1800  # 30 * 60 seconds
        assert d["partition"] == "cpuq"
        assert d["qos"] == "default"

    def test_empty_partition_and_qos_omitted(self):
        spec = ResourceSpec(cores=4, partition="", qos="")
        d = _resource_spec_to_dict(spec)
        assert "partition" not in d
        assert "qos" not in d

    def test_memory_ceiled_to_integer(self):
        spec = ResourceSpec(memory_gb=7.281214928)
        d = _resource_spec_to_dict(spec)
        assert d["memory_max"] == 8  # ceil → integer GB for SLURM


@pytest.mark.skipif(
    not importlib.util.find_spec("pysqa"),
    reason="pysqa not installed",
)
class TestPatchPysqaTemplate:
    def test_adds_qos_directive(self):
        """Patching inserts --qos into the pysqa SLURM template."""
        from pysqa.wrapper import slurm as _slurm_mod

        original = _slurm_mod.template
        try:
            # Reset to a template without qos
            _slurm_mod.template = original.replace(
                "{%- if qos %}\n#SBATCH --qos={{qos}}\n{%- endif %}\n", ""
            )
            assert "qos" not in _slurm_mod.template
            _patch_pysqa_template()
            assert "#SBATCH --qos={{qos}}" in _slurm_mod.template
        finally:
            _slurm_mod.template = original

    def test_idempotent(self):
        """Calling _patch_pysqa_template twice doesn't duplicate the directive."""
        from pysqa.wrapper import slurm as _slurm_mod

        original = _slurm_mod.template
        try:
            _patch_pysqa_template()
            after_first = _slurm_mod.template
            _patch_pysqa_template()
            assert _slurm_mod.template == after_first
        finally:
            _slurm_mod.template = original


class TestCreateExecutor:
    def test_disabled_returns_local(self):
        cfg = ClusterConfig(enabled=False)
        ex = create_executor(cfg)
        assert isinstance(ex, LocalExecutor)

    def test_enabled_without_executorlib_raises(self):
        """When executorlib is not installed, exits with error."""
        cfg = ClusterConfig(enabled=True)
        with patch.dict("sys.modules", {"executorlib": None}):
            with pytest.raises(SystemExit):
                create_executor(cfg)

@pytest.mark.skipif(
    not importlib.util.find_spec("executorlib"),
    reason="executorlib not installed",
)
class TestSbatchRetry:
    """Tests for the subprocess.check_output retry wrapper.

    The retry patch replaces the ``subprocess`` module in the scheduler's
    namespace so that ``check_output`` retries on ``CalledProcessError``.
    This works even when pysqa stores a direct reference to
    ``pysqa_execute_command`` because the function still looks up
    ``subprocess`` from its own module globals on every call.
    """

    def test_retries_on_failure_then_succeeds(self):
        """Retry logic retries on CalledProcessError and returns on success."""
        import subprocess as _subprocess
        from executorlib.standalone import scheduler as _sched_mod

        # Clean slate: remove any prior patch
        _orig_subprocess = _sched_mod.subprocess
        if hasattr(_sched_mod, "_sbatch_retry_patched"):
            delattr(_sched_mod, "_sbatch_retry_patched")
        _sched_mod.subprocess = _subprocess

        call_count = 0
        _real_check_output = _subprocess.check_output

        def _flaky(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            if call_count < 3:
                raise _subprocess.CalledProcessError(1, "sbatch")
            return b"12345\n"

        _subprocess.check_output = _flaky
        try:
            _patch_sbatch_retry(max_retries=5, base_delay=0.01)
            # The patched check_output retries _flaky via the module replacement
            result = _sched_mod.subprocess.check_output(["sbatch", "test.sh"])
            assert result == b"12345\n"
            assert call_count == 3
        finally:
            _subprocess.check_output = _real_check_output
            if hasattr(_sched_mod, "_sbatch_retry_patched"):
                delattr(_sched_mod, "_sbatch_retry_patched")
            _sched_mod.subprocess = _orig_subprocess

    def test_raises_after_max_retries(self):
        """Raises CalledProcessError after exhausting retries."""
        import subprocess as _subprocess
        from executorlib.standalone import scheduler as _sched_mod

        _orig_subprocess = _sched_mod.subprocess
        if hasattr(_sched_mod, "_sbatch_retry_patched"):
            delattr(_sched_mod, "_sbatch_retry_patched")
        _sched_mod.subprocess = _subprocess

        _real_check_output = _subprocess.check_output

        def _always_fail(*args, **kwargs):
            raise _subprocess.CalledProcessError(1, "sbatch")

        _subprocess.check_output = _always_fail
        try:
            _patch_sbatch_retry(max_retries=2, base_delay=0.01)
            with pytest.raises(_subprocess.CalledProcessError):
                _sched_mod.subprocess.check_output(["sbatch", "test.sh"])
        finally:
            _subprocess.check_output = _real_check_output
            if hasattr(_sched_mod, "_sbatch_retry_patched"):
                delattr(_sched_mod, "_sbatch_retry_patched")
            _sched_mod.subprocess = _orig_subprocess


# ===================================================================
# resource_estimation.py tests
# ===================================================================
class TestFileHelpers:
    def test_file_size_nonexistent(self):
        assert _file_size_bytes(Path("/nonexistent/path.fa")) == 0

    def test_file_size_real(self, tmp_path):
        p = tmp_path / "test.fa"
        p.write_text(">seq\nACGT\n")
        assert _file_size_bytes(p) == 10

    def test_estimate_bp_plain_fasta(self, tmp_path):
        p = tmp_path / "genome.fa"
        p.write_bytes(b"x" * 1000)  # 1KB
        bp = _estimate_genome_bp_from_filesize(p)
        assert bp == 500  # 50% of 1000

    def test_estimate_bp_gzip(self, tmp_path):
        p = tmp_path / "genome.fa.gz"
        p.write_bytes(b"x" * 1000)
        bp = _estimate_genome_bp_from_filesize(p)
        assert bp == 1500  # 3 * 0.5 * 1000

    def test_estimate_bp_nonexistent(self):
        assert _estimate_genome_bp_from_filesize(Path("/no/file")) == 0


class TestScaleTime:
    def test_base_case(self):
        # genome == scale_bp → factor == 1 → returns base × safety
        from final_finalizer.utils.resource_estimation import _TIME_SAFETY_FACTOR as SF
        assert _scale_time(60, 500_000_000) == 60 * SF

    def test_scales_up(self):
        # 2x genome → at least 2x base × safety
        from final_finalizer.utils.resource_estimation import _TIME_SAFETY_FACTOR as SF
        assert _scale_time(60, 1_000_000_000) == 120 * SF

    def test_minimum(self):
        # Small genome → still returns base_minutes × safety
        from final_finalizer.utils.resource_estimation import _TIME_SAFETY_FACTOR as SF
        assert _scale_time(60, 1000) == 60 * SF

    def test_zero_bp(self):
        from final_finalizer.utils.resource_estimation import _TIME_SAFETY_FACTOR as SF
        assert _scale_time(60, 0) == 60 * SF


class TestEstimateSyntenyResources:
    def test_protein_mode(self, tmp_path):
        ref = tmp_path / "ref.fa"
        qry = tmp_path / "qry.fa"
        ref.write_bytes(b"x" * 100_000)
        qry.write_bytes(b"x" * 200_000)
        cfg = ClusterConfig(max_threads=32)
        spec = estimate_synteny_resources(ref, qry, "protein", cfg)
        assert spec.cores <= 32
        assert spec.memory_gb >= 1.0
        assert spec.time_minutes >= 30
        assert spec.job_name == "synteny"

    def test_nucleotide_mode(self, tmp_path):
        ref = tmp_path / "ref.fa"
        qry = tmp_path / "qry.fa"
        ref.write_bytes(b"x" * 100_000)
        qry.write_bytes(b"x" * 100_000)
        cfg = ClusterConfig(max_threads=16)
        spec = estimate_synteny_resources(ref, qry, "nucleotide", cfg)
        assert spec.cores <= 16
        assert spec.memory_gb >= 4.0
        assert spec.job_name == "synteny"

    def test_clamped_to_config(self, tmp_path):
        ref = tmp_path / "ref.fa"
        qry = tmp_path / "qry.fa"
        ref.write_bytes(b"x" * 100)
        qry.write_bytes(b"x" * 100)
        cfg = ClusterConfig(max_threads=4, max_mem_gb=3.0, max_time_minutes=10)
        spec = estimate_synteny_resources(ref, qry, "protein", cfg)
        assert spec.cores <= 4
        assert spec.memory_gb <= 3.0
        assert spec.time_minutes <= 10


class TestEstimateBlastResources:
    def test_basic(self, tmp_path):
        qry = tmp_path / "qry.fa"
        qry.write_bytes(b"x" * 50_000)
        cfg = ClusterConfig()
        spec = estimate_blast_resources(qry, cfg, job_name="organelle")
        assert spec.cores <= 8
        assert spec.memory_gb == 2.0
        assert spec.job_name == "organelle"


class TestEstimateDebrisResources:
    def test_basic(self, tmp_path):
        qry = tmp_path / "qry.fa"
        qry.write_bytes(b"x" * 50_000)
        cfg = ClusterConfig()
        spec = estimate_debris_resources(qry, cfg)
        assert spec.job_name == "debris"
        assert spec.memory_gb >= 1.0


class TestEstimateContaminantResources:
    def test_with_index(self, tmp_path):
        idx_prefix = str(tmp_path / "cfr_idx")
        idx_file = tmp_path / "cfr_idx.1.cfr"
        idx_file.write_bytes(b"x" * 4_000_000_000)  # 4GB index
        cfg = ClusterConfig()
        spec = estimate_contaminant_resources(idx_prefix, cfg)
        assert spec.memory_gb > 4.0  # Should be index_size * 1.2 + 2
        assert spec.job_name == "contaminant"

    def test_without_index(self, tmp_path):
        idx_prefix = str(tmp_path / "missing_idx")
        cfg = ClusterConfig()
        spec = estimate_contaminant_resources(idx_prefix, cfg)
        assert spec.memory_gb >= 4.0  # Minimum


class TestEstimatePairwiseResources:
    def test_basic(self, tmp_path):
        left = tmp_path / "left.fa"
        right = tmp_path / "right.fa"
        left.write_bytes(b"x" * 100_000)
        right.write_bytes(b"x" * 200_000)
        cfg = ClusterConfig(max_threads=32)
        spec = estimate_pairwise_resources(left, right, cfg)
        assert spec.cores <= 32
        assert spec.memory_gb >= 4.0
        assert spec.time_minutes >= 30
        assert spec.job_name == "pairwise"

    def test_clamped_to_config(self, tmp_path):
        left = tmp_path / "left.fa"
        right = tmp_path / "right.fa"
        left.write_bytes(b"x" * 100)
        right.write_bytes(b"x" * 100)
        cfg = ClusterConfig(max_threads=4, max_mem_gb=3.0, max_time_minutes=10)
        spec = estimate_pairwise_resources(left, right, cfg)
        assert spec.cores <= 4
        assert spec.memory_gb <= 3.0
        assert spec.time_minutes <= 10

    def test_nonexistent_files(self):
        cfg = ClusterConfig()
        spec = estimate_pairwise_resources(Path("/no/left.fa"), Path("/no/right.fa"), cfg)
        assert spec.memory_gb >= 4.0  # Minimum
        assert spec.job_name == "pairwise"


class TestEstimateDepthResources:
    def test_basic(self, tmp_path):
        reads = tmp_path / "reads.fq.gz"
        assembly = tmp_path / "asm.fa"
        reads.write_bytes(b"x" * 1_000_000)
        assembly.write_bytes(b"x" * 100_000)
        cfg = ClusterConfig()
        spec = estimate_depth_resources(reads, assembly, cfg)
        assert spec.cores <= 32
        assert spec.memory_gb >= 32.0
        assert spec.job_name == "depth"
