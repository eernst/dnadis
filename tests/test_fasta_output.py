"""Tests for classified FASTA outputs: BGZF writing, file layout, symlinks."""
from __future__ import annotations

import gzip
import shutil
import subprocess

import pytest

from dnadis.models import ContigClassification
from dnadis.output.fasta_output import (
    find_stale_outputs,
    write_classified_fastas,
    write_per_subgenome_chrs_fastas,
)
from dnadis.utils.io_utils import BgzfWriter
from dnadis.utils.sequence_utils import (
    fasta_has_records,
    read_fasta_sequences,
    reverse_complement,
    write_fasta_if_changed,
)


def _clf(name, new_name, cls, ref=None, **kw):
    return ContigClassification(
        original_name=name, new_name=new_name, classification=cls,
        reversed=False, cobiont_taxid=None, cobiont_sci=None,
        assigned_ref_id=ref, ref_gene_proportion=None, contig_len=0, **kw,
    )


@pytest.fixture
def query(tmp_path):
    seqs = {
        "ctg1": "ACGTTGCA" * 50,
        "ctg2": "GGGCCCAT" * 40,
        "ctg3": "TTTTAAAAC" * 30,
        "circ1": "ATATGCGC" * 20,
        "deb1": "CCCCGGGGA" * 10,
    }
    path = tmp_path / "query.fa"
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return path, seqs


def test_bgzf_roundtrip_and_determinism(tmp_path):
    seqs = {f"s{i}": "ACGTN" * (20000 * (i + 1)) for i in range(4)}
    one = tmp_path / "one.fasta.gz"
    four = tmp_path / "four.fasta.gz"
    assert write_fasta_if_changed(seqs, one, threads=1)
    assert write_fasta_if_changed(seqs, four, threads=4)
    assert one.read_bytes() == four.read_bytes()
    assert read_fasta_sequences(one) == seqs
    # Unchanged content is not rewritten
    assert not write_fasta_if_changed(seqs, one)
    assert write_fasta_if_changed({"s0": "A"}, one)
    if shutil.which("bgzip"):
        assert subprocess.run(["bgzip", "-t", str(four)]).returncode == 0


def test_empty_bgzf_has_no_records(tmp_path):
    path = tmp_path / "empty.fasta.gz"
    with BgzfWriter(path):
        pass
    assert gzip.open(path).read() == b""
    assert path.stat().st_size > 0
    assert not fasta_has_records(path)
    assert not fasta_has_records(tmp_path / "missing.fasta.gz")


def test_classified_fasta_layout(tmp_path, query):
    qpath, seqs = query
    clfs = [
        _clf("ctg1", "chr1A", "chrom_assigned", "chr1A"),
        _clf("ctg2", "chr2A", "chrom_assigned", "chr2A"),
        _clf("ctg3", "chr2A_f1", "chrom_fragment", "chr2A"),
        _clf("circ1", "contig_1", "circular_element"),
        _clf("deb1", "contig_2", "chrom_debris"),
    ]
    prefix = tmp_path / "asm"
    paths = write_classified_fastas(
        qpath, clfs, {"ctg2": True}, prefix, scaffolded_seqs=None,
    )
    assert paths["chrs"].name == "asm.chrs.fasta.gz"
    assert paths["circles"].name == "asm.circles.fasta.gz"
    assert "chrs_scaffolds" not in paths
    chrs = read_fasta_sequences(paths["chrs"])
    assert chrs == {
        "chr1A": seqs["ctg1"],
        "chr2A": reverse_complement(seqs["ctg2"]),
        "chr2A_f1": seqs["ctg3"],
    }
    assert read_fasta_sequences(paths["circles"]) == {"contig_1": seqs["circ1"]}
    assert read_fasta_sequences(paths["debris"]) == {"contig_2": seqs["deb1"]}
    assert not fasta_has_records(paths["cobionts"])

    plain = write_classified_fastas(qpath, clfs, {}, tmp_path / "p", compress=False)
    assert plain["chrs"].name == "p.chrs.fasta"


def test_scaffolds_link_to_chrs_when_identical(tmp_path, query):
    qpath, seqs = query
    clfs = [
        _clf("ctg1", "chr1A", "chrom_assigned", "chr1A"),
        _clf("ctg2", "chr2A", "chrom_assigned", "chr2A"),
    ]
    prefix = tmp_path / "asm"
    trivial = {"chr2A": seqs["ctg2"], "chr1A": seqs["ctg1"]}
    paths = write_classified_fastas(qpath, clfs, {}, prefix, scaffolded_seqs=trivial)
    scaf = paths["chrs_scaffolds"]
    assert scaf.name == "asm.chrs.scaffolds.fasta.gz"
    assert scaf.is_symlink() and scaf.resolve() == paths["chrs"].resolve()

    # A real join replaces the link with its own file; chrs keeps the contigs.
    joined = {"chr1A": seqs["ctg1"] + "N" * 100 + seqs["ctg2"]}
    paths = write_classified_fastas(qpath, clfs, {}, prefix, scaffolded_seqs=joined)
    scaf = paths["chrs_scaffolds"]
    assert not scaf.is_symlink()
    assert read_fasta_sequences(scaf) == joined
    assert set(read_fasta_sequences(paths["chrs"])) == {"chr1A", "chr2A"}


def test_per_subgenome_single_set_is_symlink(tmp_path):
    chrs = tmp_path / "asm.chrs.fasta.gz"
    write_fasta_if_changed({"chr1A": "ACGT", "chr2A": "GGCC"}, chrs)
    out = write_per_subgenome_chrs_fastas(chrs, tmp_path / "asm")
    assert list(out) == ["A"]
    assert out["A"].name == "asm.chrs.A.fasta.gz"
    assert out["A"].is_symlink() and out["A"].resolve() == chrs.resolve()


def test_per_subgenome_multiple_sets_are_files(tmp_path):
    chrs = tmp_path / "asm.chrs.scaffolds.fasta.gz"
    write_fasta_if_changed({"chr1A": "ACGT", "chr1B": "GGCC", "chr2A": "TTAA"}, chrs)
    out = write_per_subgenome_chrs_fastas(chrs, tmp_path / "asm", kind="chrs.scaffolds")
    assert sorted(out) == ["A", "B"]
    assert out["A"].name == "asm.chrs.scaffolds.A.fasta.gz"
    assert not out["A"].is_symlink()
    assert read_fasta_sequences(out["A"]) == {"chr1A": "ACGT", "chr2A": "TTAA"}
    assert read_fasta_sequences(out["B"]) == {"chr1B": "GGCC"}


def test_find_stale_outputs(tmp_path):
    prefix = tmp_path / "asm"
    names = [
        "asm.chrs.fasta.gz", "asm.chrs.fasta", "asm.chrs.A.fasta",
        "asm.scaffolded.fasta", "asm.scaffolded.agp", "asm.circular_elements.fasta",
        "asm.rdna_consensus.fasta", "asm.contig_summary.tsv", "asm2.chrs.fasta",
    ]
    for n in names:
        (tmp_path / n).write_text("")
    stale = find_stale_outputs(prefix, [tmp_path / "asm.chrs.fasta.gz"])
    assert [p.name for p in stale] == [
        "asm.chrs.A.fasta", "asm.chrs.fasta", "asm.circular_elements.fasta",
        "asm.scaffolded.agp", "asm.scaffolded.fasta",
    ]


def test_scaffolds_keep_fragments_beside_t2t(tmp_path):
    """A fragment that is not joined to its chromosome's T2T contig is still
    emitted in the scaffold set, under its contig name."""
    from dnadis.output.scaffolding import scaffold_chromosomes

    seqs = {"t2t": "ACGT" * 500, "frag": "TTGCA" * 100}
    qfa = tmp_path / "q.fa"
    qfa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    rfa = tmp_path / "r.fa"
    rfa.write_text(">chr1A\n" + "A" * 3000 + "\n")
    clfs = [
        _clf("t2t", "chr1A", "chrom_assigned", "chr1A", is_full_length=True,
             has_5p_telomere=True, has_3p_telomere=True,
             query_subgenome_grp=1, seq_identity_vs_ref=0.95),
        _clf("frag", "chr1A_f1", "chrom_fragment", "chr1A",
             query_subgenome_grp=1, seq_identity_vs_ref=0.95),
    ]
    best_ref = {"t2t": "chr1A", "frag": "chr1A"}
    scaffolded, agp, _ = scaffold_chromosomes(
        query_fasta=qfa, classifications=clfs,
        contig_orientations={"t2t": False, "frag": True},
        qr_ref_ranges={(k, v): (1, 2000) for k, v in best_ref.items()},
        ref_fasta=rfa, ref_lengths={"chr1A": 3000},
        best_ref=best_ref, contig_refs={k: {v} for k, v in best_ref.items()},
        work_dir=tmp_path / "scaf",
        qr_best_chain_ident={(k, v): 0.95 for k, v in best_ref.items()},
    )
    assert scaffolded == {"chr1A": seqs["t2t"], "chr1A_f1": reverse_complement(seqs["frag"])}
    assert {line.split("\t")[5] for line in agp} == {"t2t", "frag"}


def test_no_join_contigs_are_scaffolded_alone(tmp_path, monkeypatch):
    """A no_join contig (sister-placed) is never joined to its chromosome's
    other contigs; it is emitted on its own, oriented as given."""
    from dnadis.output import scaffolding

    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs = {"anchor": "ACGT" * 500, "frag": "TTGCA" * 100, "sister": "GGATC" * 80}
    qfa = tmp_path / "q.fa"
    qfa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    rfa = tmp_path / "r.fa"
    rfa.write_text(">chr1A\n" + "A" * 3000 + "\n")
    clfs = [
        _clf("anchor", "chr1A", "chrom_assigned", "chr1A", is_full_length=True,
             query_subgenome_grp=1, seq_identity_vs_ref=0.95),
        _clf("frag", "chr1A_f1", "chrom_fragment", "chr1A",
             query_subgenome_grp=1, seq_identity_vs_ref=0.95),
        _clf("sister", "chr1A_f2", "chrom_fragment", "chr1A",
             query_subgenome_grp=1, seq_identity_vs_ref=0.95),
    ]
    best_ref = {c.original_name: "chr1A" for c in clfs}
    kwargs = dict(
        query_fasta=qfa, classifications=clfs,
        contig_orientations={"anchor": False, "frag": False, "sister": True},
        qr_ref_ranges={("anchor", "chr1A"): (1, 2000), ("frag", "chr1A"): (2100, 2500),
                       ("sister", "chr1A"): (1, 10)},
        ref_fasta=rfa, ref_lengths={"chr1A": 3000},
        best_ref=best_ref, contig_refs={k: {v} for k, v in best_ref.items()},
        qr_best_chain_ident={(k, v): 0.95 for k, v in best_ref.items()},
        gap_size=10,
    )

    # Without no_join, all three are joined into one scaffold.
    joined, _, _ = scaffolding.scaffold_chromosomes(work_dir=tmp_path / "a", **kwargs)
    assert set(joined) == {"chr1A"}

    scaffolded, agp, _ = scaffolding.scaffold_chromosomes(
        work_dir=tmp_path / "b", no_join={"sister"}, **kwargs,
    )
    assert scaffolded == {
        "chr1A": seqs["anchor"] + "N" * 10 + seqs["frag"],
        "chr1A_f2": reverse_complement(seqs["sister"]),
    }
    sister_lines = [line for line in agp if line.split("\t")[5] == "sister"]
    assert len(sister_lines) == 1 and sister_lines[0].split("\t")[0] == "chr1A_f2"


def test_fragments_copying_own_chromosome_are_not_joined(tmp_path, monkeypatch, caplog):
    """A chrom_fragment that is a copy of sequence in a chromosome-length
    contig of the same assembly (haplotig: both copies at about half depth,
    or a near-identical duplicate) is emitted on its own, not joined; a
    full-depth fragment matching at the same identity (homeologous region)
    is still joined."""
    import logging

    from dnadis.output import scaffolding

    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs = {
        "anchor": "ACGT" * 500, "hap": "TTGCA" * 100,
        "dup": "GGATC" * 80, "homeo": "CATG" * 90,
    }
    qfa = tmp_path / "q.fa"
    qfa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    rfa = tmp_path / "r.fa"
    rfa.write_text(">chr1A\n" + "A" * 3000 + "\n")
    frag = dict(query_subgenome_grp=1, seq_identity_vs_ref=0.95,
                redundancy_target="chr1A", redundancy_region="1-500")
    clfs = [
        _clf("anchor", "chr1A", "chrom_assigned", "chr1A", is_full_length=True,
             has_5p_telomere=True, query_subgenome_grp=1, seq_identity_vs_ref=0.95),
        _clf("hap", "chr1A_f1", "chrom_fragment", "chr1A", redundancy_cov=0.71,
             redundancy_identity=0.963, depth_ratio=0.46, target_depth_ratio=0.44, **frag),
        _clf("dup", "chr1A_f2", "chrom_fragment", "chr1A", redundancy_cov=0.96,
             redundancy_identity=0.9995, **frag),
        _clf("homeo", "chr1A_f3", "chrom_fragment", "chr1A", redundancy_cov=0.71,
             redundancy_identity=0.963, depth_ratio=0.98, target_depth_ratio=1.02, **frag),
    ]
    best_ref = {c.original_name: "chr1A" for c in clfs}
    with caplog.at_level(logging.INFO):
        scaffolded, agp, _ = scaffolding.scaffold_chromosomes(
            query_fasta=qfa, classifications=clfs,
            contig_orientations={k: False for k in seqs},
            qr_ref_ranges={("anchor", "chr1A"): (1, 2000), ("hap", "chr1A"): (2100, 2400),
                           ("dup", "chr1A"): (2500, 2600), ("homeo", "chr1A"): (2700, 2900)},
            ref_fasta=rfa, ref_lengths={"chr1A": 3000},
            best_ref=best_ref, contig_refs={k: {v} for k, v in best_ref.items()},
            work_dir=tmp_path / "scaf",
            qr_best_chain_ident={(k, v): 0.95 for k, v in best_ref.items()},
            gap_size=10,
        )
    assert scaffolded == {
        "chr1A": seqs["anchor"] + "N" * 10 + seqs["homeo"],
        "chr1A_f1": seqs["hap"],
        "chr1A_f2": seqs["dup"],
    }
    assert {line.split("\t")[5] for line in agp if line.split("\t")[4] == "W"} == set(seqs)
    assert [c.scaffold_exclusion for c in clfs] == [
        None, "own_assembly_copy", "own_assembly_copy", None,
    ]
    assert caplog.text.count("chr1A_f1: copy of chr1A:1-500 (coverage=0.71") == 1
    assert "chr1A_f2: copy of chr1A:1-500" in caplog.text
    assert "chr1A_f3: copy of" not in caplog.text


def _copy_frag(name, ref, target, **kw):
    ev = dict(redundancy_target=target, redundancy_region="1001-101000",
              redundancy_cov=0.8, redundancy_identity=0.96,
              depth_ratio=0.5, target_depth_ratio=0.5, contig_len=100_000)
    ev.update(kw)
    return ContigClassification(
        original_name=name, new_name=f"{ref}_f_{name}", classification="chrom_fragment",
        reversed=False, cobiont_taxid=None, cobiont_sci=None,
        assigned_ref_id=ref, ref_gene_proportion=None, **ev,
    )


def test_own_assembly_copy_rules():
    """Haplotig copies count only against a contig on the fragment's own
    reference chromosome; near-identical copies count against any; a copy
    overlapping a chromosome end and extending past it is excluded too."""
    from dnadis.output.scaffolding import _own_assembly_copies

    def chrom(name, ref, tel5=True, tel3=True):
        return ContigClassification(
            original_name=name, new_name=name, classification="chrom_assigned",
            reversed=False, cobiont_taxid=None, cobiont_sci=None,
            assigned_ref_id=ref, ref_gene_proportion=None, contig_len=5_000_000,
            has_5p_telomere=tel5, has_3p_telomere=tel3,
        )

    clfs = [
        chrom("chr1A", "chr1A"), chrom("chr2A", "chr2A"),
        chrom("chr3A", "chr3A", tel3=False),
        # haplotig of its own chromosome: excluded
        _copy_frag("own", "chr1A", "chr1A"),
        # half-depth match to another chromosome (gametolog/paralog): joinable
        _copy_frag("other", "chr1A", "chr2A"),
        # near-identical match to another chromosome: excluded
        _copy_frag("dup", "chr1A", "chr2A", redundancy_identity=0.9995),
        # half-depth fragment, full-depth region: joinable
        _copy_frag("homeo", "chr1A", "chr1A", target_depth_ratio=1.0),
        # overlaps the telomere-less 3' end of chr3A, 20 kb past it: excluded
        _copy_frag("ext", "chr3A", "chr3A", redundancy_region="4920001-5000000"),
        # fragment with no reference assignment uses best_ref
        _copy_frag("bestref", None, "chr1A"),
    ]
    out = _own_assembly_copies(clfs, best_ref={"bestref": "chr1A"})
    assert set(out) == {"own", "dup", "ext", "bestref"}
    assert out["own"].startswith("chr1A:1001-101000 (coverage=0.80, identity=0.9600")


def test_own_assembly_copies_in_t2t_and_fragment_only_groups(tmp_path, monkeypatch):
    """An own-assembly copy is kept out of a fragment-only group (the other
    fragments are still joined) and is emitted on its own beside a T2T copy."""
    from dnadis.output import scaffolding

    monkeypatch.setattr(scaffolding, "have_exe", lambda exe: False)
    seqs = {"t2t": "ACGT" * 500, "hapT": "TTGCA" * 100,
            "f1": "GGATC" * 80, "f2": "CATG" * 90, "hapF": "AACCG" * 70}
    qfa = tmp_path / "q.fa"
    qfa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    rfa = tmp_path / "r.fa"
    rfa.write_text(">chr1A\n" + "A" * 3000 + "\n>chr2A\n" + "C" * 3000 + "\n")
    hap = dict(redundancy_cov=0.9, redundancy_identity=0.96, depth_ratio=0.5,
               target_depth_ratio=0.5, redundancy_region="1-400", query_subgenome_grp=1)
    clfs = [
        _clf("t2t", "chr1A", "chrom_assigned", "chr1A", is_full_length=True,
             has_5p_telomere=True, has_3p_telomere=True, query_subgenome_grp=1),
        _clf("hapT", "chr1A_f1", "chrom_fragment", "chr1A", redundancy_target="chr1A", **hap),
        # chr2A survives only as fragments; hapF copies chr1A at near identity
        _clf("f1", "chr2A_f1", "chrom_fragment", "chr2A", query_subgenome_grp=1),
        _clf("f2", "chr2A_f2", "chrom_fragment", "chr2A", query_subgenome_grp=1),
        _clf("hapF", "chr2A_f3", "chrom_fragment", "chr2A", redundancy_target="chr1A",
             **{**hap, "redundancy_identity": 0.9995}),
    ]
    best_ref = {c.original_name: c.assigned_ref_id for c in clfs}
    scaffolded, agp, _ = scaffolding.scaffold_chromosomes(
        query_fasta=qfa, classifications=clfs,
        contig_orientations={k: False for k in seqs},
        qr_ref_ranges={("t2t", "chr1A"): (1, 2000), ("hapT", "chr1A"): (1, 400),
                       ("f1", "chr2A"): (1, 400), ("f2", "chr2A"): (1000, 1400),
                       ("hapF", "chr2A"): (2000, 2400)},
        ref_fasta=rfa, ref_lengths={"chr1A": 3000, "chr2A": 3000},
        best_ref=best_ref, contig_refs={k: {v} for k, v in best_ref.items()},
        work_dir=tmp_path / "scaf",
        qr_best_chain_ident={(k, v): 0.95 for k, v in best_ref.items()},
        gap_size=10,
    )
    assert scaffolded == {
        "chr1A": seqs["t2t"],
        "chr1A_f1": seqs["hapT"],
        "chr2A": seqs["f1"] + "N" * 10 + seqs["f2"],
        "chr2A_f3": seqs["hapF"],
    }
    assert [c.original_name for c in clfs if c.scaffold_exclusion] == ["hapT", "hapF"]
