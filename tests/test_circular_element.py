"""Tests for circular-element detection and classification."""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnadis.models import ContigClassification
from dnadis.classification.classifier import apply_circular_reclassification
from dnadis.utils.sequence_utils import resolve_circular_contigs


def _clf(name: str, classification: str, ref: str | None = None) -> ContigClassification:
    return ContigClassification(
        original_name=name,
        new_name="",
        classification=classification,
        reversed=False,
        cobiont_taxid=None,
        cobiont_sci=None,
        assigned_ref_id=ref,
        ref_gene_proportion=None,
        contig_len=1000,
    )


# ---------------------------------------------------------------------------
# apply_circular_reclassification
# ---------------------------------------------------------------------------
class TestCircularReclassification:
    def test_reroutes_eligible_categories(self):
        clfs = [
            _clf("f", "chrom_fragment", "chr1"),
            _clf("cd", "chrom_debris", "chr1"),
            _clf("cu", "chrom_unassigned"),
            _clf("d", "debris"),
            _clf("u", "unclassified"),
        ]
        n = apply_circular_reclassification(clfs, {"f", "cd", "cu", "d", "u"})
        assert n == 5
        assert all(c.classification == "circular_element" for c in clfs)
        assert all(c.is_circular for c in clfs)
        # assigned_ref_id preserved as a "homologous to" annotation
        assert clfs[0].assigned_ref_id == "chr1"
        assert all(c.classification_confidence == "high" for c in clfs)

    def test_preserves_organelle_cobiont_rdna_chromosome(self):
        # A circular chloroplast / mito / cobiont / rDNA / real chromosome must
        # keep its classification even though it is flagged circular.
        preserved = ["organelle_complete", "organelle_debris", "rDNA",
                     "cobiont", "chrom_assigned"]
        clfs = [_clf(name, name) for name in preserved]
        n = apply_circular_reclassification(clfs, {c.original_name for c in clfs})
        assert n == 0
        assert [c.classification for c in clfs] == preserved
        assert all(c.is_circular for c in clfs)  # still flagged circular

    def test_noncircular_contigs_flagged_false_not_rerouted(self):
        clfs = [_clf("f", "chrom_fragment", "chr1")]
        n = apply_circular_reclassification(clfs, set())  # nothing circular
        assert n == 0
        assert clfs[0].classification == "chrom_fragment"
        assert clfs[0].is_circular is False


# ---------------------------------------------------------------------------
# resolve_circular_contigs
# ---------------------------------------------------------------------------
class TestResolveCircularContigs:
    def test_circular_list_source(self, tmp_path):
        q = tmp_path / "asm.fasta"
        q.write_text(">a\nAC\n>b\nGT\n")
        lst = tmp_path / "circ.txt"
        lst.write_text("# names\nb\nmissing\n")
        names, known = resolve_circular_contigs(
            q, {"a", "b"}, circular_list=lst)
        assert known is True
        assert names == {"b"}  # intersected with real names

    def test_circular_fasta_source(self, tmp_path):
        q = tmp_path / "asm.fasta"
        q.write_text(">a\nAC\n>b\nGT\n")
        cf = tmp_path / "circ.fasta"
        cf.write_text(">b desc\nGT\n")
        names, known = resolve_circular_contigs(
            q, {"a", "b"}, circular_fasta=cf)
        assert known is True and names == {"b"}

    def test_sibling_autodetect_gz(self, tmp_path):
        q = tmp_path / "x.bp.p_ctg.fasta.gz"
        with gzip.open(q, "wt") as fh:
            fh.write(">a\nAC\n>b\nGT\n")
        sib = tmp_path / "x.bp.p_ctg.circ.fasta.gz"
        with gzip.open(sib, "wt") as fh:
            fh.write(">b\nGT\n")
        names, known = resolve_circular_contigs(q, {"a", "b"})
        assert known is True and names == {"b"}

    def test_hifiasm_name_heuristic(self, tmp_path):
        q = tmp_path / "asm.fasta"
        q.write_text(">ptg000001l\nAC\n>ptg000002c\nGT\n")
        names, known = resolve_circular_contigs(q, {"ptg000001l", "ptg000002c"})
        assert known is True and names == {"ptg000002c"}

    def test_unknown_when_no_source(self, tmp_path):
        q = tmp_path / "asm.fasta"
        q.write_text(">chr1\nAC\n>chr2\nGT\n")
        names, known = resolve_circular_contigs(q, {"chr1", "chr2"})
        assert known is False and names == set()
