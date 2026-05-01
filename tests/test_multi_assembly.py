"""Tests for multi-assembly input parsing utilities."""
import os
import tempfile
from pathlib import Path

import pytest

from dnadis.utils.multi_assembly import (
    _strip_fasta_extension,
    parse_fofn,
    scan_assembly_dir,
    resolve_assemblies,
)


# ---------------------------------------------------------------------------
# Helper to create temp FASTA files
# ---------------------------------------------------------------------------
def _make_fasta(path: Path, seqname: str = "chr1"):
    """Write a minimal FASTA file at *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{seqname}\nACGTACGT\n")


# ---------------------------------------------------------------------------
# _strip_fasta_extension
# ---------------------------------------------------------------------------
class TestStripFastaExtension:
    def test_fasta(self):
        assert _strip_fasta_extension("sample.fasta") == "sample"

    def test_fa(self):
        assert _strip_fasta_extension("sample.fa") == "sample"

    def test_fna(self):
        assert _strip_fasta_extension("sample.fna") == "sample"

    def test_fasta_gz(self):
        assert _strip_fasta_extension("sample.fasta.gz") == "sample"

    def test_fa_gz(self):
        assert _strip_fasta_extension("sample.fa.gz") == "sample"

    def test_no_extension(self):
        assert _strip_fasta_extension("sample") == "sample"

    def test_complex_name(self):
        assert _strip_fasta_extension("my.assembly.v2.fasta") == "my.assembly.v2"


# ---------------------------------------------------------------------------
# parse_fofn
# ---------------------------------------------------------------------------
class TestParseFofn:
    def test_basic_path_only(self, tmp_path):
        """FOFN with only 'path' column."""
        asm1 = tmp_path / "asm1.fasta"
        asm2 = tmp_path / "asm2.fa"
        _make_fasta(asm1)
        _make_fasta(asm2)

        fofn = tmp_path / "assemblies.tsv"
        fofn.write_text(f"path\n{asm1}\n{asm2}\n")

        result = parse_fofn(fofn)
        assert len(result) == 2
        assert result[0] == (asm1, "asm1", None)
        assert result[1] == (asm2, "asm2", None)

    def test_with_name_and_reads(self, tmp_path):
        """FOFN with path, name, and reads columns."""
        asm = tmp_path / "raw.fasta"
        reads = tmp_path / "reads.fastq"
        _make_fasta(asm)
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        fofn = tmp_path / "assemblies.tsv"
        fofn.write_text(f"path\tname\treads\n{asm}\tcool_asm\t{reads}\n")

        result = parse_fofn(fofn)
        assert len(result) == 1
        assert result[0] == (asm, "cool_asm", reads)

    def test_relative_paths(self, tmp_path):
        """Relative paths are resolved relative to the FOFN directory."""
        subdir = tmp_path / "data"
        subdir.mkdir()
        asm = subdir / "asm.fasta"
        _make_fasta(asm)

        fofn = subdir / "list.tsv"
        fofn.write_text("path\nasm.fasta\n")

        result = parse_fofn(fofn)
        assert result[0][0] == asm.resolve()

    def test_headerless_path_only(self, tmp_path):
        """Headerless FOFN with just paths (one per line)."""
        asm1 = tmp_path / "asm1.fasta"
        asm2 = tmp_path / "asm2.fa"
        _make_fasta(asm1)
        _make_fasta(asm2)

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"{asm1}\n{asm2}\n")

        result = parse_fofn(fofn)
        assert len(result) == 2
        assert result[0] == (asm1, "asm1", None)
        assert result[1] == (asm2, "asm2", None)

    def test_headerless_with_name_and_reads(self, tmp_path):
        """Headerless FOFN where columns 1 and 2 are interpreted as name/reads."""
        asm = tmp_path / "raw.fasta"
        reads = tmp_path / "reads.fastq"
        _make_fasta(asm)
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"{asm}\tcool_asm\t{reads}\n")

        result = parse_fofn(fofn)
        assert result[0] == (asm, "cool_asm", reads)

    def test_missing_file(self, tmp_path):
        """FOFN referencing a non-existent file raises ValueError."""
        fofn = tmp_path / "bad.tsv"
        fofn.write_text(f"path\n{tmp_path / 'nonexistent.fa'}\n")

        with pytest.raises(ValueError, match="not found"):
            parse_fofn(fofn)

    def test_duplicate_names(self, tmp_path):
        """Duplicate assembly names raise ValueError."""
        asm1 = tmp_path / "a" / "sample.fasta"
        asm2 = tmp_path / "b" / "sample.fasta"
        _make_fasta(asm1)
        _make_fasta(asm2)

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"path\n{asm1}\n{asm2}\n")

        with pytest.raises(ValueError, match="duplicate"):
            parse_fofn(fofn)

    def test_empty_fofn(self, tmp_path):
        """FOFN with header but no entries raises ValueError."""
        fofn = tmp_path / "empty.tsv"
        fofn.write_text("path\n")

        with pytest.raises(ValueError, match="no assembly entries"):
            parse_fofn(fofn)

    def test_comment_lines_skipped(self, tmp_path):
        """Lines starting with # are skipped."""
        asm = tmp_path / "asm.fasta"
        _make_fasta(asm)

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"path\n# this is a comment\n{asm}\n")

        result = parse_fofn(fofn)
        assert len(result) == 1

    def test_empty_reads_field(self, tmp_path):
        """Empty reads field results in None."""
        asm = tmp_path / "asm.fasta"
        _make_fasta(asm)

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"path\tname\treads\n{asm}\ttest\t\n")

        result = parse_fofn(fofn)
        assert result[0][2] is None


# ---------------------------------------------------------------------------
# scan_assembly_dir
# ---------------------------------------------------------------------------
class TestScanAssemblyDir:
    def test_finds_fasta_files(self, tmp_path):
        """Finds .fasta, .fa, .fna files."""
        _make_fasta(tmp_path / "one.fasta")
        _make_fasta(tmp_path / "two.fa")
        _make_fasta(tmp_path / "three.fna")
        # Non-FASTA file should be ignored
        (tmp_path / "readme.txt").write_text("hello")

        result = scan_assembly_dir(tmp_path)
        assert len(result) == 3
        names = [n for _, n in result]
        assert sorted(names) == ["one", "three", "two"]

    def test_strips_extensions(self, tmp_path):
        """Names are derived by stripping FASTA extensions."""
        _make_fasta(tmp_path / "genomeA.fasta.gz")
        _make_fasta(tmp_path / "genomeB.fa")

        result = scan_assembly_dir(tmp_path)
        names = sorted(n for _, n in result)
        assert names == ["genomeA", "genomeB"]

    def test_empty_directory(self, tmp_path):
        """Empty directory raises ValueError."""
        with pytest.raises(ValueError, match="No FASTA files"):
            scan_assembly_dir(tmp_path)

    def test_not_a_directory(self, tmp_path):
        """Non-directory path raises ValueError."""
        f = tmp_path / "file.txt"
        f.write_text("hi")
        with pytest.raises(ValueError, match="Not a directory"):
            scan_assembly_dir(f)

    def test_name_collision_uses_stem(self, tmp_path):
        """When stripped names collide, falls back to full stem."""
        # Both would strip to "sample" — collision
        _make_fasta(tmp_path / "sample.fasta")
        _make_fasta(tmp_path / "sample.fa")

        result = scan_assembly_dir(tmp_path)
        names = sorted(n for _, n in result)
        # Falls back to stem: "sample" for both .fasta and .fa
        # Then numeric suffix to break tie
        assert len(set(names)) == 2  # Names must be unique


# ---------------------------------------------------------------------------
# resolve_assemblies
# ---------------------------------------------------------------------------
class TestResolveAssemblies:
    def test_fofn_mode(self, tmp_path):
        """resolve_assemblies dispatches to parse_fofn."""
        asm = tmp_path / "asm.fasta"
        _make_fasta(asm)

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"path\n{asm}\n")

        class Args:
            pass

        args = Args()
        args.fofn = fofn
        args.assembly_dir = None
        args.reads = None

        result = resolve_assemblies(args)
        assert len(result) == 1
        assert result[0][1] == "asm"

    def test_dir_mode(self, tmp_path):
        """resolve_assemblies dispatches to scan_assembly_dir."""
        _make_fasta(tmp_path / "a.fasta")
        _make_fasta(tmp_path / "b.fasta")

        class Args:
            pass

        args = Args()
        args.fofn = None
        args.assembly_dir = tmp_path
        args.reads = None

        result = resolve_assemblies(args)
        assert len(result) == 2

    def test_global_reads_fallback(self, tmp_path):
        """Global --reads is used as fallback for FOFN entries without reads."""
        asm = tmp_path / "asm.fasta"
        reads = tmp_path / "reads.fq"
        _make_fasta(asm)
        reads.write_text("@r\nACGT\n+\nIIII\n")

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"path\n{asm}\n")

        class Args:
            pass

        args = Args()
        args.fofn = fofn
        args.assembly_dir = None
        args.reads = reads

        result = resolve_assemblies(args)
        assert result[0][2] == reads

    def test_fofn_reads_override_global(self, tmp_path):
        """Per-entry reads in FOFN take priority over global --reads."""
        asm = tmp_path / "asm.fasta"
        global_reads = tmp_path / "global.fq"
        local_reads = tmp_path / "local.fq"
        _make_fasta(asm)
        global_reads.write_text("@r\nACGT\n+\nIIII\n")
        local_reads.write_text("@r\nACGT\n+\nIIII\n")

        fofn = tmp_path / "list.tsv"
        fofn.write_text(f"path\treads\n{asm}\t{local_reads}\n")

        class Args:
            pass

        args = Args()
        args.fofn = fofn
        args.assembly_dir = None
        args.reads = global_reads

        result = resolve_assemblies(args)
        assert result[0][2] == local_reads.resolve()

    def test_dir_mode_with_global_reads(self, tmp_path):
        """Directory mode applies global reads to all assemblies."""
        _make_fasta(tmp_path / "a.fasta")
        reads = tmp_path / "reads.fq"
        reads.write_text("@r\nACGT\n+\nIIII\n")

        class Args:
            pass

        args = Args()
        args.fofn = None
        args.assembly_dir = tmp_path
        args.reads = reads

        result = resolve_assemblies(args)
        assert result[0][2] == reads
