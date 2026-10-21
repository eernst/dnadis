"""Tests for ``dnadis cleanup``."""
from __future__ import annotations

from pathlib import Path

import pytest

from dnadis.cleanup import build_plan, main, run_cleanup


def _touch(path: Path, size: int = 10) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)
    return path


@pytest.fixture
def results(tmp_path):
    root = tmp_path / "out"
    asm = root / "asm1"
    work = asm / "asm1_classification"
    lin = work / "compleasm_full" / "liliopsida_odb12"
    remove = [
        work / "non_chromosome_contigs.fa",
        work / "chr_debris" / "chromosome_contigs.fa",
        work / "chr_debris" / "debris_candidates.fa",
        work / "scaffolding" / "chr1A" / "contigs.fa",
        work / "rdna_consensus" / "reannotation" / "query_db.nsq",
        work / "cobionts" / "contigs.fq",
        work / "read_depth" / "depth.mosdepth.global.dist.txt",
        work / "read_depth" / "aligned_reads.sorted.bam",
        lin / "miniprot_output.gff",
        lin / "hmmer_output" / "x.tbl",
        lin / "miniprot.done",
        root / "pairwise" / "a_vs_b.paf.gz",
    ]
    keep = [
        asm / "asm1.chrs.fasta.gz",
        asm / "asm1.contig_summary.tsv",
        work / "chr_debris" / "debris_vs_chrs.paf",
        work / "cobionts" / "centrifuger_results.tsv",
        work / "read_depth" / "depth.regions.bed.gz",
        work / "read_depth" / "depth.regions.bed.gz.csi",
        work / "read_depth" / "alignment_metadata.json",
        work / "rdna_consensus" / "cache_key.txt",
        work / "compleasm_full" / "summary.txt",
        lin / "full_table_busco_format.tsv",
        lin / "translated_protein.fasta",
        root / "pairwise" / "a_vs_b.macro_blocks.tsv",
        root / "sister_anchoring" / "anchors" / "asm1.fa",
        root / "sister_anchoring" / "candidates.fa",
        root / "sister_anchoring" / "sister_anchoring.tsv",
        root / "scaffold_join_check" / "notes.fa",
    ]
    for p in remove + keep:
        _touch(p)
    return root, remove, keep


def test_dry_run_removes_nothing(results):
    root, remove, keep = results
    freed = run_cleanup(root, dry_run=True)
    assert freed == 10 * len(remove)
    assert all(p.exists() for p in remove + keep)


def test_cleanup_removes_intermediates_only(results):
    root, remove, keep = results
    run_cleanup(root)
    assert not [p for p in remove if p.exists()]
    assert all(p.exists() for p in keep)
    assert not (root / "asm1" / "asm1_classification" / "scaffolding").exists()


def test_sister_anchors_opt_in(results):
    root, _, _ = results
    run_cleanup(root, sister_anchors=True)
    assert not (root / "sister_anchoring" / "anchors").exists()
    assert not (root / "sister_anchoring" / "candidates.fa").exists()
    assert (root / "sister_anchoring" / "sister_anchoring.tsv").exists()


def test_legacy_outputs_need_replacement(results):
    root, _, _ = results
    asm = root / "asm1"
    replaced = _touch(asm / "asm1.chrs.fasta")              # chrs.fasta.gz exists
    orphan = _touch(asm / "asm1.debris.fasta")              # no .gz replacement
    old_agp = _touch(asm / "asm1.scaffolded.agp")
    consensus = _touch(asm / "asm1.rdna_consensus.fasta")
    plan = build_plan(root, legacy_outputs=True)
    assert plan.groups["legacy-name outputs"] == [replaced]
    _touch(asm / "asm1.chrs.scaffolds.agp")
    plan = build_plan(root, legacy_outputs=True)
    assert set(plan.groups["legacy-name outputs"]) == {replaced, old_agp}
    assert orphan.exists() and consensus.exists()


def test_symlinked_dir_is_unlinked_not_followed(results, tmp_path):
    root, _, _ = results
    target = tmp_path / "elsewhere"
    _touch(target / "precious.fa")
    work = root / "asm1" / "asm1_classification"
    (work / "organelles").symlink_to(target, target_is_directory=True)
    run_cleanup(root)
    assert not (work / "organelles").is_symlink()
    assert (target / "precious.fa").exists()


def test_keep_bams(results):
    root, _, _ = results
    bam = root / "asm1" / "asm1_classification" / "read_depth" / "aligned_reads.sorted.bam"
    bai = _touch(Path(f"{bam}.bai"))
    run_cleanup(root, keep_bams=True)
    assert bam.exists() and bai.exists()
    run_cleanup(root)
    assert not bam.exists() and not bai.exists()


def test_rejects_non_results_dir(tmp_path):
    with pytest.raises(SystemExit):
        main([str(tmp_path)])
