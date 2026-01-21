import re
import textwrap
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import final_finalizer as ff


def test_normalize_ref_id_case_variations():
    """Test that normalize_ref_id handles various case patterns."""
    # Lowercase already - unchanged
    assert ff.normalize_ref_id("chr1") == "chr1"
    assert ff.normalize_ref_id("chrC") == "chrC"

    # Uppercase Chr -> chr
    assert ff.normalize_ref_id("Chr1") == "chr1"
    assert ff.normalize_ref_id("ChrC") == "chrC"
    assert ff.normalize_ref_id("Chr10") == "chr10"

    # All uppercase CHR -> chr
    assert ff.normalize_ref_id("CHR1") == "chr1"
    assert ff.normalize_ref_id("CHRM") == "chrM"

    # Non-chr prefixes unchanged
    assert ff.normalize_ref_id("scaffold_1") == "scaffold_1"
    assert ff.normalize_ref_id("contig123") == "contig123"
    assert ff.normalize_ref_id("Mt") == "Mt"


def test_normalize_organelle_id_aliases():
    """Test that normalize_organelle_id recognizes various naming conventions."""
    # Standard names
    assert ff.normalize_organelle_id("chrC") == "chrC"
    assert ff.normalize_organelle_id("chrM") == "chrM"
    assert ff.normalize_organelle_id("ChrC") == "chrC"
    assert ff.normalize_organelle_id("ChrM") == "chrM"

    # Common aliases
    assert ff.normalize_organelle_id("Mt") == "chrM"
    assert ff.normalize_organelle_id("Pt") == "chrC"
    assert ff.normalize_organelle_id("CpDNA") == "chrC"
    assert ff.normalize_organelle_id("Chloroplast") == "chrC"
    assert ff.normalize_organelle_id("Mitochondrion") == "chrM"
    assert ff.normalize_organelle_id("mitochondria") == "chrM"

    # Non-organelles return None
    assert ff.normalize_organelle_id("chr1") is None
    assert ff.normalize_organelle_id("Chr5") is None
    assert ff.normalize_organelle_id("scaffold_1") is None


def test_set_ref_id_patterns_and_reset():
    """Test that custom patterns can be set and reset."""
    # Save original state
    original_patterns = ff.get_ref_id_patterns()

    try:
        # Set custom patterns
        custom = [re.compile(r"^(?P<chrom>scaffold\d+)(?P<sg>[XY])$")]
        ff.set_ref_id_patterns(custom)
        assert ff.get_ref_id_patterns() == custom

        # Test split_chrom_subgenome uses custom patterns
        assert ff.split_chrom_subgenome("scaffold1X") == ("scaffold1", "X")
        assert ff.split_chrom_subgenome("scaffold2Y") == ("scaffold2", "Y")

        # Non-matching IDs fall through
        assert ff.split_chrom_subgenome("chr1A") == ("chr1A", "NA")

        # Reset to defaults
        ff.set_ref_id_patterns(None)
        assert ff.get_ref_id_patterns() == ff._REF_ID_PATTERNS

        # Verify default patterns work again
        assert ff.split_chrom_subgenome("chr1A") == ("chr1", "A")
    finally:
        # Restore original state
        ff.set_ref_id_patterns(None)


def test_split_chrom_subgenome_with_explicit_patterns():
    """Test that split_chrom_subgenome accepts explicit patterns parameter."""
    # Custom pattern for cotton-style naming
    cotton_patterns = [
        re.compile(r"^(?P<chrom>chr\d+)(?P<sg>At|Dt)$", re.IGNORECASE),
    ]

    assert ff.split_chrom_subgenome("chr1At", patterns=cotton_patterns) == ("chr1", "At")
    assert ff.split_chrom_subgenome("chr5Dt", patterns=cotton_patterns) == ("chr5", "Dt")

    # Without explicit patterns, uses global/default
    assert ff.split_chrom_subgenome("chr1A") == ("chr1", "A")


def test_split_chrom_subgenome_arabidopsis_conventions():
    assert ff.split_chrom_subgenome("chr1") == ("chr1", "NA")
    assert ff.split_chrom_subgenome("chr1A") == ("chr1", "A")
    assert ff.split_chrom_subgenome("chr5B") == ("chr5", "B")
    assert ff.split_chrom_subgenome("chrC") == ("chrC", "NA")
    assert ff.split_chrom_subgenome("chrM") == ("chrM", "NA")
    assert ff.split_chrom_subgenome("scaffold_12") == ("scaffold_12", "NA")


def test_read_fasta_lengths_and_ref_lengths_tsv(tmp_path):
    fasta_path = tmp_path / "arabidopsis_mock.fa"
    fasta_path.write_text(
        textwrap.dedent(
            """
            >chr1
            ACGTACGTAC
            >chr2
            ACGTAC
            >chrC
            ACGT
            >chrM
            ACGTA
            """
        ).strip()
        + "\n"
    )

    lengths = ff.read_fasta_lengths(fasta_path)
    assert lengths["chr1"] == 10
    assert lengths["chr2"] == 6
    assert lengths["chrC"] == 4
    assert lengths["chrM"] == 5

    assert ff.get_min_nuclear_chrom_length(lengths) == 6

    tsv_path = tmp_path / "ref_lengths.tsv"
    ff.write_ref_lengths_tsv(tsv_path, fasta_path)
    tsv = tsv_path.read_text().splitlines()
    assert tsv[0].split("\t") == ["ref_id", "chrom_id", "subgenome", "ref_len"]
    assert "chr1\tchr1\tNA\t10" in tsv
    assert "chrC\tchrC\tNA\t4" in tsv


def test_parse_gff3_transcript_coords_araport11(tmp_path):
    gff3_path = tmp_path / "araport11_mock.gff3"
    gff3_path.write_text(
        textwrap.dedent(
            """
            ##gff-version 3
            chr1\tAraport11\tmRNA\t10\t20\t.\t+\t.\tID=AT1G01010.1;Parent=AT1G01010
            chr1\tAraport11\ttranscript\t30\t40\t.\t-\t.\tID=AT1G01020.1
            """
        ).strip()
        + "\n"
    )

    tx2loc, tx2gene = ff.parse_gff3_transcript_coords(gff3_path)
    assert tx2loc["AT1G01010.1"] == ("chr1", 9, 20, "+")
    assert tx2gene["AT1G01010.1"] == "AT1G01010"

    assert tx2loc["AT1G01020.1"] == ("chr1", 29, 40, "-")
    assert tx2gene["AT1G01020.1"] == "AT1G01020.1"


def test_parse_blast_coverage(tmp_path):
    blast_path = tmp_path / "blast.tsv"
    blast_path.write_text(
        textwrap.dedent(
            """
            contigA\tsubject1\t99.0\t21\t0\t0\t10\t30\t1\t21\t1e-20\t100\t3702
            contigA\tsubject2\t98.0\t11\t0\t0\t40\t50\t5\t15\t1e-10\t80\t3702
            contigB\tsubject3\t97.0\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\t3702
            """
        ).strip()
        + "\n"
    )
    query_lengths = {"contigA": 100, "contigB": 50}

    results = ff.parse_blast_coverage(blast_path, query_lengths)

    contigA = results["contigA"]
    assert contigA.total_aligned_bp == 30
    assert contigA.total_coverage == 0.3
    assert contigA.best_hit_subject == "subject1"
    assert contigA.best_hit_taxid == 3702

    contigB = results["contigB"]
    assert contigB.total_aligned_bp == 9
    assert contigB.total_coverage == 0.18


def test_parse_paf_chain_evidence_for_subgenome_assignment(tmp_path):
    paf_path = tmp_path / "arabidopsis_mock.paf"
    paf_path.write_text(
        textwrap.dedent(
            """
            contigA\t1000\t0\t200\t+\tchr1A\t1000\t0\t200\t180\t200\t60\ttp:A:P
            contigA\t1000\t300\t500\t+\tchr1A\t1000\t300\t500\t180\t200\t60\ttp:A:P
            contigA\t1000\t0\t150\t+\tchr1B\t1000\t0\t150\t120\t150\t60\ttp:A:P
            """
        ).strip()
        + "\n"
    )
    contig_lengths = {"contigA": 1000}

    ev = ff.parse_paf_chain_evidence_and_segments(
        paf_gz_path=paf_path,
        contig_lengths=contig_lengths,
        assign_minlen=100,
        assign_minmapq=20,
        assign_tp="P",
        chain_q_gap=200,
        chain_r_gap=200,
        chain_diag_slop=10,
        assign_min_ident=0.8,
        assign_chain_topk=1,
        assign_chain_score="matches",
        assign_chain_min_bp=100,
        assign_ref_score="topk",
    )

    assert ev.best_ref["contigA"] == "chr1A"
    assert ev.best_bp["contigA"] == 400
    assert ev.contig_total["contigA"] == 550


def test_filter_overlapping_hits_by_identity():
    """Test that overlapping hits are filtered, keeping highest identity."""
    # Create overlapping blocks from two different ref_ids (simulating homeologs)
    blocks = {
        ("contigA", "chr1A", "+"): [
            ff.Block(qs=100, qe=200, rs=1000, re=1100, matches=95, aln_len=100, mapq=60, strand="+", gene_id="geneA1"),
            ff.Block(qs=300, qe=400, rs=2000, re=2100, matches=90, aln_len=100, mapq=60, strand="+", gene_id="geneA2"),
        ],
        ("contigA", "chr1P", "+"): [
            # Overlaps with geneA1 at same position but lower identity
            ff.Block(qs=100, qe=200, rs=1000, re=1100, matches=85, aln_len=100, mapq=60, strand="+", gene_id="geneP1"),
            # Non-overlapping hit should be kept
            ff.Block(qs=500, qe=600, rs=3000, re=3100, matches=80, aln_len=100, mapq=60, strand="+", gene_id="geneP2"),
        ],
    }

    filtered, n_before, n_removed, removed_per_contig = ff._filter_overlapping_hits_by_identity(blocks)

    assert n_before == 4
    assert n_removed == 1  # geneP1 should be removed (overlaps geneA1 with lower identity)
    assert removed_per_contig.get("contigA", 0) == 1

    # chr1A should keep both its hits
    assert len(filtered[("contigA", "chr1A", "+")]) == 2

    # chr1P should keep only geneP2 (non-overlapping)
    assert len(filtered[("contigA", "chr1P", "+")]) == 1
    assert filtered[("contigA", "chr1P", "+")][0].gene_id == "geneP2"


def test_filter_overlapping_hits_higher_identity_wins():
    """Test that when P subgenome has higher identity, it wins."""
    blocks = {
        ("contigA", "chr1A", "+"): [
            ff.Block(qs=100, qe=200, rs=1000, re=1100, matches=80, aln_len=100, mapq=60, strand="+", gene_id="geneA1"),
        ],
        ("contigA", "chr1P", "+"): [
            # Overlaps with geneA1 but has HIGHER identity - should win
            ff.Block(qs=100, qe=200, rs=1000, re=1100, matches=95, aln_len=100, mapq=60, strand="+", gene_id="geneP1"),
        ],
    }

    filtered, n_before, n_removed, removed_per_contig = ff._filter_overlapping_hits_by_identity(blocks)

    assert n_before == 2
    assert n_removed == 1  # geneA1 should be removed
    assert removed_per_contig.get("contigA", 0) == 1

    # chr1P should keep its hit (higher identity)
    assert len(filtered[("contigA", "chr1P", "+")]) == 1

    # chr1A should have no hits left
    assert len(filtered.get(("contigA", "chr1A", "+"), [])) == 0
