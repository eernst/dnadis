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


# ----------------------------
# Telomere detection tests
# ----------------------------
def test_telomere_detection_plant_motif():
    """Test telomere detection with plant telomere motif TTTAGGG."""
    from final_finalizer.detection.telomere import detect_telomeres_single, _count_telomere_repeats

    # Test _count_telomere_repeats directly
    # Plant telomere motif: TTTAGGG (forward) / CCCTAAA (reverse)
    motif = "TTTAGGG"

    # Sequence with 5 forward repeats
    seq_with_telomere = "TTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGG"
    count = _count_telomere_repeats(seq_with_telomere, motif, min_repeats=3)
    assert count == 5, f"Expected 5 repeats, got {count}"

    # Sequence with reverse complement repeats
    seq_with_rc_telomere = "CCCTAAACCCTAAACCCTAAACCCTAAA"
    count_rc = _count_telomere_repeats(seq_with_rc_telomere, motif, min_repeats=3)
    assert count_rc == 4, f"Expected 4 RC repeats, got {count_rc}"

    # Sequence with only 2 repeats (below threshold)
    seq_short = "TTTAGGGTTTAGGG"
    count_short = _count_telomere_repeats(seq_short, motif, min_repeats=3)
    assert count_short == 0, f"Expected 0 (below threshold), got {count_short}"


def test_telomere_detection_at_contig_ends():
    """Test telomere detection at contig 5' and 3' ends."""
    from final_finalizer.detection.telomere import detect_telomeres_single

    # Create a synthetic contig with telomeres at both ends
    # 5' end: CCCTAAA repeats (reverse complement telomere on + strand)
    # 3' end: TTTAGGG repeats (forward telomere)
    telomere_5p = "CCCTAAA" * 10  # 70 bp
    middle = "A" * 500
    telomere_3p = "TTTAGGG" * 10  # 70 bp
    contig = telomere_5p + middle + telomere_3p

    result = detect_telomeres_single(
        sequence=contig,
        motif="TTTAGGG",
        window_size=100,
        min_repeats=3,
    )

    assert result.has_5p_telomere, "Should detect 5' telomere"
    assert result.has_3p_telomere, "Should detect 3' telomere"
    assert result.telomere_5p_count >= 3
    assert result.telomere_3p_count >= 3


def test_telomere_detection_no_telomeres():
    """Test that contigs without telomeres are correctly identified."""
    from final_finalizer.detection.telomere import detect_telomeres_single

    # Random sequence without telomeres
    contig = "ACGTACGTACGT" * 100

    result = detect_telomeres_single(
        sequence=contig,
        motif="TTTAGGG",
        window_size=100,
        min_repeats=3,
    )

    assert not result.has_5p_telomere
    assert not result.has_3p_telomere
    assert result.telomere_5p_count == 0
    assert result.telomere_3p_count == 0


# ----------------------------
# Full-length classification tests
# ----------------------------
def test_classify_full_length_both_telomeres():
    """Test that both telomeres = high confidence full-length."""
    from final_finalizer.classification.classifier import classify_full_length

    # Both telomeres - should be full-length with high confidence regardless of coverage
    is_full, confidence = classify_full_length(
        ref_coverage=0.50,  # Low coverage
        has_5p_telomere=True,
        has_3p_telomere=True,
        full_length_threshold=0.70,
    )
    assert is_full is True
    assert confidence == "high"


def test_classify_full_length_one_telomere_high_coverage():
    """Test that one telomere + high coverage = medium confidence full-length."""
    from final_finalizer.classification.classifier import classify_full_length

    is_full, confidence = classify_full_length(
        ref_coverage=0.80,  # Above threshold
        has_5p_telomere=True,
        has_3p_telomere=False,
        full_length_threshold=0.70,
    )
    assert is_full is True
    assert confidence == "medium"


def test_classify_full_length_one_telomere_low_coverage():
    """Test that one telomere + low coverage = medium confidence fragment."""
    from final_finalizer.classification.classifier import classify_full_length

    is_full, confidence = classify_full_length(
        ref_coverage=0.50,  # Below threshold
        has_5p_telomere=True,
        has_3p_telomere=False,
        full_length_threshold=0.70,
    )
    assert is_full is False
    assert confidence == "medium"


def test_classify_full_length_no_telomeres_high_coverage():
    """Test that no telomeres + high coverage = low confidence full-length."""
    from final_finalizer.classification.classifier import classify_full_length

    is_full, confidence = classify_full_length(
        ref_coverage=0.85,  # Above threshold
        has_5p_telomere=False,
        has_3p_telomere=False,
        full_length_threshold=0.70,
    )
    assert is_full is True
    assert confidence == "low"


def test_classify_full_length_no_telomeres_low_coverage():
    """Test that no telomeres + low coverage = low confidence fragment."""
    from final_finalizer.classification.classifier import classify_full_length

    is_full, confidence = classify_full_length(
        ref_coverage=0.40,  # Below threshold
        has_5p_telomere=False,
        has_3p_telomere=False,
        full_length_threshold=0.70,
    )
    assert is_full is False
    assert confidence == "low"


# ----------------------------
# Naming scheme tests
# ----------------------------
def test_generate_contig_names_full_length_single():
    """Test naming for single full-length contig."""
    from final_finalizer.classification.classifier import generate_contig_names
    from final_finalizer.models import ContigClassification

    clf = ContigClassification(
        original_name="contig_001",
        new_name="",
        classification="chrom_assigned",
        reversed=False,
        contaminant_taxid=None,
        contaminant_sci=None,
        assigned_ref_id="chr1A",
        ref_gene_proportion=0.5,
        contig_len=1000000,
        is_full_length=True,
        query_subgenome=None,
    )

    name_mapping = generate_contig_names(
        classifications=[clf],
        query_lengths={"contig_001": 1000000},
        add_subgenome_suffix=None,
    )

    assert name_mapping["contig_001"] == "chr1A"


def test_generate_contig_names_fragment():
    """Test naming for fragment contigs."""
    from final_finalizer.classification.classifier import generate_contig_names
    from final_finalizer.models import ContigClassification

    clfs = [
        ContigClassification(
            original_name="contig_001",
            new_name="",
            classification="chrom_assigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chr1A",
            ref_gene_proportion=0.3,
            contig_len=500000,
            is_full_length=False,  # Fragment
            query_subgenome=None,
        ),
        ContigClassification(
            original_name="contig_002",
            new_name="",
            classification="chrom_assigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chr1A",
            ref_gene_proportion=0.2,
            contig_len=300000,
            is_full_length=False,  # Fragment
            query_subgenome=None,
        ),
    ]

    name_mapping = generate_contig_names(
        classifications=clfs,
        query_lengths={"contig_001": 500000, "contig_002": 300000},
        add_subgenome_suffix=None,
    )

    # Both are fragments, should get _f1, _f2 suffix by length
    assert name_mapping["contig_001"] == "chr1A_f1"
    assert name_mapping["contig_002"] == "chr1A_f2"


def test_generate_contig_names_query_subgenome():
    """Test naming with query subgenome suffix."""
    from final_finalizer.classification.classifier import generate_contig_names
    from final_finalizer.models import ContigClassification

    clfs = [
        ContigClassification(
            original_name="contig_001",
            new_name="",
            classification="chrom_assigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chr1A",
            ref_gene_proportion=0.5,
            contig_len=1000000,
            is_full_length=True,
            query_subgenome=None,  # Primary
        ),
        ContigClassification(
            original_name="contig_002",
            new_name="",
            classification="chrom_assigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chr1A",
            ref_gene_proportion=0.4,
            contig_len=900000,
            is_full_length=True,
            query_subgenome="B",  # Query subgenome B
        ),
    ]

    name_mapping = generate_contig_names(
        classifications=clfs,
        query_lengths={"contig_001": 1000000, "contig_002": 900000},
        add_subgenome_suffix=None,
    )

    assert name_mapping["contig_001"] == "chr1A"  # Primary, no suffix
    assert name_mapping["contig_002"] == "chr1A_B"  # Query subgenome B


def test_generate_contig_names_multiple_full_length_copies():
    """Test naming when multiple full-length copies exist (unusual case)."""
    from final_finalizer.classification.classifier import generate_contig_names
    from final_finalizer.models import ContigClassification

    clfs = [
        ContigClassification(
            original_name="contig_001",
            new_name="",
            classification="chrom_assigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chr1A",
            ref_gene_proportion=0.5,
            contig_len=1000000,
            is_full_length=True,
            query_subgenome=None,
        ),
        ContigClassification(
            original_name="contig_002",
            new_name="",
            classification="chrom_assigned",
            reversed=False,
            contaminant_taxid=None,
            contaminant_sci=None,
            assigned_ref_id="chr1A",
            ref_gene_proportion=0.5,
            contig_len=900000,
            is_full_length=True,  # Also full-length!
            query_subgenome=None,
        ),
    ]

    name_mapping = generate_contig_names(
        classifications=clfs,
        query_lengths={"contig_001": 1000000, "contig_002": 900000},
        add_subgenome_suffix=None,
    )

    # Multiple full-length copies get _c1, _c2 suffix
    assert name_mapping["contig_001"] == "chr1A_c1"
    assert name_mapping["contig_002"] == "chr1A_c2"


# ----------------------------
# Rearrangement hypothesis detection tests (cluster-based)
# ----------------------------
def test_compute_largest_cluster_single_ref():
    """Test cluster computation with chains to single ref within gap tolerance."""
    from final_finalizer.classification.classifier import compute_largest_cluster_metrics_per_ref

    # All chains to same ref, close together - should form one cluster
    # Row format: (contig, contig_len, ref_id, chrom_id, subgenome, strand, chain_id,
    #              qstart, qend, qspan_bp, union_bp, ..., ref_start, ref_end)
    macro_block_rows = [
        ("contig1", 1000000, "chr1A", "chr1", "A", "+", 1, 100000, 300000, 200000, 150000, 0, 0, 0, 0, 0, 0, 100000, 300000),
        ("contig1", 1000000, "chr1A", "chr1", "A", "+", 2, 500000, 700000, 200000, 180000, 0, 0, 0, 0, 0, 0, 500000, 700000),
    ]

    result = compute_largest_cluster_metrics_per_ref(macro_block_rows)

    # chr1A: one cluster spanning 100000 to 700000 = 600000, with 330000 aligned bp
    span, aligned, genes = result[("contig1", "chr1A")]
    assert span == 600000
    assert aligned == 330000
    assert genes == 0  # gene_count fields are 0 in test data


def test_compute_largest_cluster_gap_splits_clusters():
    """Test that large gaps between same-ref chains create separate clusters."""
    from final_finalizer.classification.classifier import compute_largest_cluster_metrics_per_ref

    # chr12A chains with 1 Mb gap between them (> 500kb default max_gap)
    macro_block_rows = [
        ("contig1", 5000000, "chr12A", "chr12", "A", "-", 1, 100000, 300000, 200000, 150000, 0, 0, 0, 0, 0, 0, 100000, 300000),
        ("contig1", 5000000, "chr12A", "chr12", "A", "-", 2, 1500000, 1700000, 200000, 180000, 0, 0, 0, 0, 0, 0, 1500000, 1700000),
    ]

    result = compute_largest_cluster_metrics_per_ref(macro_block_rows)

    # chr12A: two separate clusters, each with span 200000
    # Largest cluster has span 200000
    span, aligned, genes = result[("contig1", "chr12A")]
    assert span == 200000


def test_compute_largest_cluster_ignores_intervening_chains():
    """Test that chains to other refs don't affect clustering of same-ref chains."""
    from final_finalizer.classification.classifier import compute_largest_cluster_metrics_per_ref

    # chr12A chains with chr1A chain between them, but gap < max_gap
    macro_block_rows = [
        ("contig1", 1000000, "chr12A", "chr12", "A", "-", 1, 100000, 300000, 200000, 150000, 0, 0, 0, 0, 0, 0, 100000, 300000),
        ("contig1", 1000000, "chr1A", "chr1", "A", "+", 1, 400000, 500000, 100000, 80000, 0, 0, 0, 0, 0, 0, 400000, 500000),
        ("contig1", 1000000, "chr12A", "chr12", "A", "-", 2, 600000, 800000, 200000, 180000, 0, 0, 0, 0, 0, 0, 600000, 800000),
    ]

    result = compute_largest_cluster_metrics_per_ref(macro_block_rows)

    # chr12A: gap between chains is 600000-300000=300000 < 500kb, so one cluster
    # spanning 100000 to 800000 = 700000
    span, aligned, genes = result[("contig1", "chr12A")]
    assert span == 700000
    assert aligned == 330000  # 150000 + 180000


def test_compute_largest_cluster_selects_largest():
    """Test that the largest cluster is selected when multiple exist."""
    from final_finalizer.classification.classifier import compute_largest_cluster_metrics_per_ref

    # chr12A has two clusters: small one then large one separated by big gap
    macro_block_rows = [
        ("contig1", 5000000, "chr12A", "chr12", "A", "-", 1, 50000, 100000, 50000, 40000, 0, 0, 0, 0, 0, 0, 50000, 100000),
        # 1.4 Mb gap > 500kb max_gap
        ("contig1", 5000000, "chr12A", "chr12", "A", "-", 2, 1500000, 1800000, 300000, 250000, 0, 0, 0, 0, 0, 0, 1500000, 1800000),
        ("contig1", 5000000, "chr12A", "chr12", "A", "-", 3, 1850000, 2000000, 150000, 120000, 0, 0, 0, 0, 0, 0, 1850000, 2000000),
    ]

    result = compute_largest_cluster_metrics_per_ref(macro_block_rows)

    # chr12A: first cluster = 50000, second cluster spans 1500000-2000000 = 500000
    span, aligned, genes = result[("contig1", "chr12A")]
    assert span == 500000
    assert aligned == 370000  # 250000 + 120000


def test_detect_rearrangement_candidates_no_offtarget():
    """Test that no rearrangement is detected when off-target is below span threshold."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    # (span, aligned_bp, gene_count) tuples
    qr_cluster_metrics = {
        ("contig1", "chr1A"): (900000, 700000, 200),  # 90% span
        ("contig1", "chr2A"): (50000, 40000, 20),      # 5% span (below 10% threshold)
    }
    contig_refs = {"chr1A", "chr2A"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
        ref_gene_density=30.0,
    )

    assert result is None  # No off-target >= 10%


def test_detect_rearrangement_candidates_single_offtarget():
    """Test detection of single off-target with gene density above calibrated threshold."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    # chr5B: 200kb span, 50 genes = 250 genes/Mb
    # ref_gene_density=30, density_frac=0.10 → threshold=3 genes/Mb → passes
    qr_cluster_metrics = {
        ("contig1", "chr1A"): (700000, 500000, 150),  # 70% span
        ("contig1", "chr5B"): (200000, 100000, 50),    # 20% span, 250 genes/Mb
    }
    contig_refs = {"chr1A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
        ref_gene_density=30.0,
    )

    assert result == "chr5B"


def test_detect_rearrangement_candidates_low_density_rejected_nucleotide():
    """Test that off-target with large span but low density is rejected in nucleotide mode."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    qr_cluster_metrics = {
        ("contig1", "chr1A"): (700000, 500000, 0),  # 70% span
        ("contig1", "chr5B"): (200000, 20000, 0),   # 20% span but only 10% density (< 15% min)
    }
    contig_refs = {"chr1A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
        min_density=0.15,
        synteny_mode="nucleotide",
    )

    assert result is None  # Density too low


def test_detect_rearrangement_candidates_low_gene_density_rejected_protein():
    """Test that off-target with low gene density is rejected in protein mode."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    # chr5B: 5 Mb span, 4 genes = 0.8 genes/Mb (paralogous noise)
    # ref_gene_density=30, density_frac=0.10 → threshold=3 genes/Mb → rejected
    qr_cluster_metrics = {
        ("contig1", "chr1A"): (7000000, 5000000, 150),  # 70% span
        ("contig1", "chr5B"): (5000000, 20000, 4),       # 50% span but 0.8 genes/Mb
    }
    contig_refs = {"chr1A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=10000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
        synteny_mode="protein",
        ref_gene_density=30.0,
    )

    assert result is None  # Gene density too low


def test_detect_rearrangement_candidates_below_gene_floor_rejected():
    """Test that off-target with < min_genes_floor is rejected even with OK density."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    # chr5B: 200kb span, 4 genes = 20 genes/Mb (passes density but fails floor of 5)
    qr_cluster_metrics = {
        ("contig1", "chr1A"): (700000, 500000, 150),
        ("contig1", "chr5B"): (200000, 20000, 4),
    }
    contig_refs = {"chr1A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
        synteny_mode="protein",
        ref_gene_density=30.0,
    )

    assert result is None  # Below min_genes_floor=5


def test_detect_rearrangement_candidates_multiple_offtarget():
    """Test detection of multiple off-target chromosomes, ordered by span fraction."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    qr_cluster_metrics = {
        ("contig1", "chr1A"): (500000, 400000, 100),  # 50% span
        ("contig1", "chr2A"): (150000, 50000, 30),     # 15% span, 200 genes/Mb
        ("contig1", "chr5B"): (250000, 100000, 50),    # 25% span, 200 genes/Mb
    }
    contig_refs = {"chr1A", "chr2A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
        ref_gene_density=30.0,
    )

    # Should be ordered by span fraction (highest first)
    assert result == "chr5B,chr2A"


def test_detect_rearrangement_candidates_no_assignment():
    """Test that no rearrangement is detected when contig has no assignment."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="",  # No assignment
        contig_len=1000000,
        qr_cluster_metrics={},
        contig_refs=set(),
        threshold=0.10,
    )

    assert result is None


def test_compute_reference_gene_density():
    """Test median on-target gene density computation."""
    from final_finalizer.classification.classifier import compute_reference_gene_density

    # 3 contigs: densities of 20, 30, 40 genes/Mb → median = 30
    qr_cluster_metrics = {
        ("contig1", "chr1"): (1_000_000, 500000, 20),   # 20 genes/Mb
        ("contig2", "chr2"): (2_000_000, 1000000, 60),   # 30 genes/Mb
        ("contig3", "chr3"): (500_000, 250000, 20),       # 40 genes/Mb
    }
    best_ref = {"contig1": "chr1", "contig2": "chr2", "contig3": "chr3"}

    density = compute_reference_gene_density(qr_cluster_metrics, best_ref)
    assert density == 30.0


def test_compute_reference_gene_density_empty():
    """Test that empty input returns 0.0."""
    from final_finalizer.classification.classifier import compute_reference_gene_density

    assert compute_reference_gene_density({}, {}) == 0.0


# ----------------------------
# Scaffolding tests
# ----------------------------
def test_builtin_scaffold_ordering():
    """Test that built-in scaffolder orders contigs by reference midpoint."""
    from final_finalizer.output.scaffolding import _builtin_scaffold

    sequences = {
        "ctgA": "ACGT" * 25,  # 100 bp
        "ctgB": "TGCA" * 25,  # 100 bp
    }
    # ctgB has earlier ref position than ctgA
    ref_ranges = {
        ("ctgA", "chr1"): (500000, 600000),
        ("ctgB", "chr1"): (100000, 200000),
    }
    orientations = {"ctgA": False, "ctgB": False}

    seq, agp = _builtin_scaffold(
        contig_names=["ctgA", "ctgB"],
        sequences=sequences,
        orientations=orientations,
        ref_ranges=ref_ranges,
        ref_id="chr1",
        scaffold_name="chr1",
        gap_size=10,
    )

    # ctgB should come first (earlier ref position)
    assert seq.startswith("TGCA")
    assert "N" * 10 in seq
    assert len(seq) == 100 + 10 + 100

    # AGP should have 3 lines: W, N, W
    assert len(agp) == 3
    assert "\tW\t" in agp[0]
    assert "ctgB" in agp[0]
    assert "\tN\t" in agp[1]
    assert "\tW\t" in agp[2]
    assert "ctgA" in agp[2]


def test_builtin_scaffold_reverse_complement():
    """Test that built-in scaffolder correctly orients reversed contigs."""
    from final_finalizer.output.scaffolding import _builtin_scaffold

    sequences = {"ctgA": "AAACCC"}
    ref_ranges = {("ctgA", "chr1"): (100, 200)}
    orientations = {"ctgA": True}  # Reversed

    seq, agp = _builtin_scaffold(
        contig_names=["ctgA"],
        sequences=sequences,
        orientations=orientations,
        ref_ranges=ref_ranges,
        ref_id="chr1",
        scaffold_name="chr1",
        gap_size=100,
    )

    assert seq == "GGGTTT"  # Reverse complement of AAACCC
    assert len(agp) == 1
    assert agp[0].endswith("\t-")


def test_scaffold_skip_full_length():
    """Test that full-length contigs with both telomeres get trivial AGP."""
    from final_finalizer.output.scaffolding import _group_contigs_by_haplotype
    from final_finalizer.models import ContigClassification

    clf = ContigClassification(
        original_name="ctg1",
        new_name="chr1A",
        classification="chrom_assigned",
        reversed=False,
        contaminant_taxid=None,
        contaminant_sci=None,
        assigned_ref_id="chr1A",
        ref_gene_proportion=0.9,
        contig_len=30000000,
        is_full_length=True,
        has_5p_telomere=True,
        has_3p_telomere=True,
        query_subgenome_grp=1,
    )

    groups = _group_contigs_by_haplotype(
        classifications=[clf],
        best_ref={"ctg1": "chr1A"},
        contig_refs={"ctg1": {"chr1A"}},
        qr_best_chain_ident={("ctg1", "chr1A"): 0.95},
    )

    assert ("chr1A", 1) in groups
    assert groups[("chr1A", 1)] == ["ctg1"]


def test_agp_format():
    """Test AGP 2.0 output format compliance."""
    from final_finalizer.output.scaffolding import write_agp, _agp_component_line, _agp_gap_line

    agp_lines = [
        _agp_component_line("chr1", 1, 1000, 1, "ctgA", 1, 1000, "+"),
        _agp_gap_line("chr1", 1001, 1100, 2, 100),
        _agp_component_line("chr1", 1101, 2100, 3, "ctgB", 1, 1000, "-"),
    ]

    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".agp", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        write_agp(agp_lines, tmp_path)
        content = tmp_path.read_text()
        lines = content.strip().split("\n")

        # First line should be AGP version header
        assert lines[0] == "##agp-version\t2.0"

        # Component line
        fields = lines[1].split("\t")
        assert fields[0] == "chr1"
        assert fields[4] == "W"
        assert fields[5] == "ctgA"
        assert fields[8] == "+"

        # Gap line
        fields = lines[2].split("\t")
        assert fields[4] == "N"
        assert fields[5] == "100"
        assert fields[6] == "scaffold"
        assert fields[7] == "yes"
        assert fields[8] == "align_genus"

        # Component line (reverse)
        fields = lines[3].split("\t")
        assert fields[4] == "W"
        assert fields[5] == "ctgB"
        assert fields[8] == "-"
    finally:
        tmp_path.unlink()


def test_scaffold_contig_grouping():
    """Test that contigs are correctly grouped by (ref, haplotype)."""
    from final_finalizer.output.scaffolding import _group_contigs_by_haplotype
    from final_finalizer.models import ContigClassification

    clfs = [
        ContigClassification(
            original_name="ctg1", new_name="chr1A", classification="chrom_assigned",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr1A", ref_gene_proportion=0.5, contig_len=1000000,
            query_subgenome_grp=1, seq_identity_vs_ref=0.95,
        ),
        ContigClassification(
            original_name="ctg2", new_name="chr1A_f2", classification="chrom_assigned",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr1A", ref_gene_proportion=0.3, contig_len=500000,
            query_subgenome_grp=1, seq_identity_vs_ref=0.94,
        ),
        ContigClassification(
            original_name="ctg3", new_name="chr2A", classification="chrom_assigned",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr2A", ref_gene_proportion=0.8, contig_len=2000000,
            query_subgenome_grp=1, seq_identity_vs_ref=0.92,
        ),
        # chrom_debris contig with evidence pointing to chr1A (scaffolding candidate)
        ContigClassification(
            original_name="ctg4", new_name="debris_1", classification="chrom_debris",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr1A", ref_gene_proportion=None, contig_len=50000,
        ),
        # Plain debris contig (should NOT be grouped - only chrom_debris allowed)
        ContigClassification(
            original_name="ctg4b", new_name="debris_2", classification="debris",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr1A", ref_gene_proportion=None, contig_len=30000,
        ),
        # Unclassified contig with no evidence (should not be grouped)
        ContigClassification(
            original_name="ctg5", new_name="unclass_1", classification="unclassified",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=None, ref_gene_proportion=None, contig_len=10000,
        ),
    ]

    groups = _group_contigs_by_haplotype(
        classifications=clfs,
        best_ref={"ctg1": "chr1A", "ctg2": "chr1A", "ctg3": "chr2A", "ctg4": "chr1A"},
        contig_refs={"ctg1": {"chr1A"}, "ctg2": {"chr1A"}, "ctg3": {"chr2A"}, "ctg4": {"chr1A"}},
        qr_best_chain_ident={
            ("ctg1", "chr1A"): 0.95, ("ctg2", "chr1A"): 0.94,
            ("ctg3", "chr2A"): 0.92, ("ctg4", "chr1A"): 0.93,
        },
    )

    # (chr1A, 1) should have ctg1, ctg2, and chrom_debris ctg4
    assert set(groups[("chr1A", 1)]) == {"ctg1", "ctg2", "ctg4"}
    # (chr2A, 1) should have ctg3
    assert groups[("chr2A", 1)] == ["ctg3"]
    # ctg4b (plain debris) and ctg5 should not appear anywhere
    all_grouped = set()
    for contigs in groups.values():
        all_grouped.update(contigs)
    assert "ctg4b" not in all_grouped
    assert "ctg5" not in all_grouped


# ----------------------------
# Collinearity score tests
# ----------------------------
def _make_macro_row(contig, ref_id, strand, qstart, qend, union_bp, ref_start, ref_end):
    """Helper to create a macro block row tuple with correct indices."""
    return (
        contig, 1000000, ref_id, "chr1", "A", strand, 1,
        qstart, qend, qend - qstart, union_bp, 0, 0, "0.95", "100.0",
        1, 0, ref_start, ref_end,
    )


def test_compute_collinearity_score_perfect():
    """All concordant chains → 1.0."""
    from final_finalizer.alignment.chain_parsing import compute_collinearity_score

    rows = [
        _make_macro_row("ctg1", "chr1A", "+", 0, 100000, 50000, 0, 100000),
        _make_macro_row("ctg1", "chr1A", "+", 200000, 300000, 60000, 200000, 300000),
        _make_macro_row("ctg1", "chr1A", "+", 400000, 500000, 70000, 400000, 500000),
    ]
    score = compute_collinearity_score(rows, "ctg1", "chr1A")
    assert score == 1.0


def test_compute_collinearity_score_discordant():
    """All discordant chains → 0.0."""
    from final_finalizer.alignment.chain_parsing import compute_collinearity_score

    # Query order: 0-100k, 200-300k, 400-500k
    # Ref order: 500k, 300k, 100k → all decreasing on + strand → all discordant
    rows = [
        _make_macro_row("ctg1", "chr1A", "+", 0, 100000, 50000, 500000, 600000),
        _make_macro_row("ctg1", "chr1A", "+", 200000, 300000, 60000, 300000, 400000),
        _make_macro_row("ctg1", "chr1A", "+", 400000, 500000, 70000, 100000, 200000),
    ]
    score = compute_collinearity_score(rows, "ctg1", "chr1A")
    assert score == 0.0


def test_compute_collinearity_score_single_chain():
    """Single chain → 1.0 (trivially collinear)."""
    from final_finalizer.alignment.chain_parsing import compute_collinearity_score

    rows = [
        _make_macro_row("ctg1", "chr1A", "+", 0, 100000, 50000, 0, 100000),
    ]
    score = compute_collinearity_score(rows, "ctg1", "chr1A")
    assert score == 1.0


def test_compute_collinearity_score_no_chains():
    """No chains → None."""
    from final_finalizer.alignment.chain_parsing import compute_collinearity_score

    score = compute_collinearity_score([], "ctg1", "chr1A")
    assert score is None


def test_compute_collinearity_score_bp_weighted():
    """Larger chains contribute more weight."""
    from final_finalizer.alignment.chain_parsing import compute_collinearity_score

    # First pair: concordant, small chains (10k bp each) → weight = min(10k, 10k) = 10k
    # Second pair: discordant, large chains (10k, 100k bp) → weight = min(10k, 100k) = 10k
    # concordant_weight=10k, total_weight=20k → score=0.5
    rows = [
        _make_macro_row("ctg1", "chr1A", "+", 0, 50000, 10000, 0, 50000),
        _make_macro_row("ctg1", "chr1A", "+", 100000, 200000, 10000, 100000, 200000),
        _make_macro_row("ctg1", "chr1A", "+", 300000, 400000, 100000, 50000, 100000),  # ref goes backward
    ]
    score = compute_collinearity_score(rows, "ctg1", "chr1A")
    assert abs(score - 0.5) < 0.01


def test_compute_collinearity_score_mixed_strands():
    """Inversion boundary (different strands) → discordant."""
    from final_finalizer.alignment.chain_parsing import compute_collinearity_score

    # Two + strand chains (concordant), then one - strand chain (inversion → discordant)
    rows = [
        _make_macro_row("ctg1", "chr1A", "+", 0, 100000, 50000, 0, 100000),
        _make_macro_row("ctg1", "chr1A", "+", 200000, 300000, 50000, 200000, 300000),
        _make_macro_row("ctg1", "chr1A", "-", 400000, 500000, 50000, 400000, 500000),
    ]
    score = compute_collinearity_score(rows, "ctg1", "chr1A")
    # Pair 1 (+/+, concordant): weight = 50k → concordant
    # Pair 2 (+/-, different strands): weight = 50k → discordant
    # score = 50k / 100k = 0.5
    assert abs(score - 0.5) < 0.01


def test_macro_block_rows_ref_coords():
    """Test that macro_block_rows from PAF parsing include ref_start/ref_end."""
    ev = ff.parse_paf_chain_evidence_and_segments(
        paf_gz_path=Path("/dev/null"),  # Will produce empty result
        contig_lengths={},
        assign_minlen=100,
        assign_minmapq=0,
        assign_tp="ALL",
        chain_q_gap=200000,
        chain_r_gap=400000,
        chain_diag_slop=150000,
        assign_min_ident=0.0,
        assign_chain_topk=3,
        assign_chain_score="matches",
        assign_chain_min_bp=0,
        assign_ref_score="all",
    )
    # Empty PAF: no rows, but qr_ref_ranges should exist
    assert ev.qr_ref_ranges is not None
    assert isinstance(ev.qr_ref_ranges, dict)


# ----------------------------
# 45S rDNA array detection tests
# ----------------------------

def _make_rdna_locus(contig, start, end, strand="+", identity=0.95,
                     consensus_coverage=0.90, copy_type="full"):
    """Helper to create RdnaLocus for testing."""
    from final_finalizer.models import RdnaLocus
    return RdnaLocus(
        contig=contig,
        start=start,
        end=end,
        strand=strand,
        identity=identity,
        consensus_coverage=consensus_coverage,
        copy_type=copy_type,
        sub_feature_loci=[],
    )


def test_detect_arrays_basic():
    """5 loci on same contig within gap -> 1 array, correct metrics."""
    from final_finalizer.detection.rdna_consensus import _detect_arrays

    loci = [
        _make_rdna_locus("ctg1", 1000, 11000, identity=0.95, copy_type="full"),
        _make_rdna_locus("ctg1", 12000, 22000, identity=0.96, copy_type="full"),
        _make_rdna_locus("ctg1", 23000, 33000, identity=0.94, copy_type="partial"),
        _make_rdna_locus("ctg1", 34000, 44000, identity=0.97, copy_type="full"),
        _make_rdna_locus("ctg1", 45000, 55000, identity=0.93, copy_type="fragment"),
    ]
    arrays = _detect_arrays(loci, min_tandem_copies=3, max_tandem_gap=50000)
    assert len(arrays) == 1

    arr = arrays[0]
    assert arr.array_id == "array_1"
    assert arr.contig == "ctg1"
    assert arr.start == 1000
    assert arr.end == 55000
    assert arr.n_total == 5
    assert arr.n_full == 3
    assert arr.n_partial == 1
    assert arr.n_fragment == 1
    assert arr.identity_min == 0.93
    assert arr.identity_max == 0.97
    assert arr.span == 54000


def test_detect_arrays_multiple_on_contig():
    """2 clusters separated by large gap -> 2 arrays."""
    from final_finalizer.detection.rdna_consensus import _detect_arrays

    loci = [
        # Cluster 1: positions 0-30k
        _make_rdna_locus("ctg1", 0, 10000),
        _make_rdna_locus("ctg1", 11000, 21000),
        _make_rdna_locus("ctg1", 22000, 32000),
        # Cluster 2: positions 500k-530k (gap > 50k)
        _make_rdna_locus("ctg1", 500000, 510000),
        _make_rdna_locus("ctg1", 511000, 521000),
        _make_rdna_locus("ctg1", 522000, 532000),
    ]
    arrays = _detect_arrays(loci, min_tandem_copies=3, max_tandem_gap=50000)
    assert len(arrays) == 2
    assert arrays[0].array_id == "array_1"
    assert arrays[1].array_id == "array_2"
    assert arrays[0].n_total == 3
    assert arrays[1].n_total == 3


def test_detect_arrays_below_min_copies():
    """2 loci -> no arrays, array_id stays None."""
    from final_finalizer.detection.rdna_consensus import _detect_arrays

    loci = [
        _make_rdna_locus("ctg1", 0, 10000),
        _make_rdna_locus("ctg1", 11000, 21000),
    ]
    arrays = _detect_arrays(loci, min_tandem_copies=3, max_tandem_gap=50000)
    assert len(arrays) == 0
    assert all(l.array_id is None for l in loci)


def test_detect_arrays_sets_array_id():
    """Verify constituent loci get array_id set."""
    from final_finalizer.detection.rdna_consensus import _detect_arrays

    loci = [
        _make_rdna_locus("ctg1", 0, 10000),
        _make_rdna_locus("ctg1", 11000, 21000),
        _make_rdna_locus("ctg1", 22000, 32000),
    ]
    arrays = _detect_arrays(loci, min_tandem_copies=3, max_tandem_gap=50000)
    assert len(arrays) == 1
    for locus in loci:
        assert locus.array_id == "array_1"


def test_rdna_locus_is_nor_candidate_compat():
    """Property returns True when array_id set, False when None."""
    locus_with = _make_rdna_locus("ctg1", 0, 10000)
    locus_with.array_id = "array_1"
    assert locus_with.is_nor_candidate is True

    locus_without = _make_rdna_locus("ctg1", 0, 10000)
    assert locus_without.is_nor_candidate is False


def test_detect_arrays_strand_majority():
    """Mixed strands -> array strand is majority."""
    from final_finalizer.detection.rdna_consensus import _detect_arrays

    loci = [
        _make_rdna_locus("ctg1", 0, 10000, strand="+"),
        _make_rdna_locus("ctg1", 11000, 21000, strand="+"),
        _make_rdna_locus("ctg1", 22000, 32000, strand="-"),
    ]
    arrays = _detect_arrays(loci, min_tandem_copies=3, max_tandem_gap=50000)
    assert len(arrays) == 1
    assert arrays[0].strand == "+"


# ===================================================================
# GMM 1D and subgenome inference tests
# ===================================================================

def test_gmm_1d_em_single_component():
    """GMM with k=1 should recover mean and std of a single Gaussian."""
    from final_finalizer.classification.classifier import _gmm_1d_em

    data = [0.60, 0.62, 0.58, 0.61, 0.59, 0.63, 0.60, 0.61]
    weights, means, stds, ll = _gmm_1d_em(data, k=1)

    assert len(weights) == 1
    assert abs(weights[0] - 1.0) < 1e-6
    assert abs(means[0] - 0.605) < 0.01
    assert stds[0] > 0


def test_gmm_1d_em_two_components():
    """GMM with k=2 should separate two distinct clusters."""
    from final_finalizer.classification.classifier import _gmm_1d_em

    low = [0.58, 0.59, 0.60, 0.61, 0.62]
    high = [0.70, 0.71, 0.72, 0.73, 0.74]
    data = low + high

    weights, means, stds, ll = _gmm_1d_em(data, k=2)

    assert len(means) == 2
    sorted_means = sorted(means)
    assert sorted_means[0] < 0.65
    assert sorted_means[1] > 0.65


def test_gmm_1d_bic_selects_k1_for_unimodal():
    """BIC should select k=1 for unimodal data."""
    from final_finalizer.classification.classifier import _gmm_1d_bic

    data = [0.60 + i * 0.005 for i in range(20)]
    best_k, weights, means, stds, labels = _gmm_1d_bic(data, max_k=3)

    assert best_k == 1
    assert all(l == 0 for l in labels)


def test_gmm_1d_bic_selects_k2_for_bimodal():
    """BIC should select k=2 for clearly bimodal data."""
    from final_finalizer.classification.classifier import _gmm_1d_bic

    # Two well-separated clusters (mimics real subgenome identity distributions)
    low = [0.58, 0.59, 0.60, 0.60, 0.61, 0.61, 0.62, 0.58, 0.59, 0.60,
           0.61, 0.62, 0.59, 0.60, 0.61, 0.60, 0.59, 0.61]
    high = [0.68, 0.69, 0.70, 0.70, 0.71, 0.72, 0.69, 0.70, 0.71, 0.68,
            0.69, 0.70, 0.71, 0.72, 0.70, 0.69, 0.71, 0.70]
    data = low + high

    best_k, weights, means, stds, labels = _gmm_1d_bic(data, max_k=3)

    assert best_k == 2
    # Each data point should be assigned to one of two clusters
    assert set(labels) == {0, 1}


def test_gmm_1d_bic_handles_small_data():
    """BIC should return k=1 for very small datasets."""
    from final_finalizer.classification.classifier import _gmm_1d_bic

    best_k, weights, means, stds, labels = _gmm_1d_bic([0.65], max_k=3)
    assert best_k == 1

    best_k, weights, means, stds, labels = _gmm_1d_bic([0.60, 0.70], max_k=3)
    assert best_k in (1, 2)  # Either is reasonable for n=2


def test_infer_query_subgenomes_single_contig_per_ref():
    """Single contig per ref chromosome -> no subgenome suffix."""
    from final_finalizer.classification.classifier import infer_query_subgenomes
    from final_finalizer.models import ContigClassification

    clfs = [
        ContigClassification(
            original_name="ctg1", new_name="", classification="chrom_assigned",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr1A", ref_gene_proportion=0.5, contig_len=10_000_000,
        ),
        ContigClassification(
            original_name="ctg2", new_name="", classification="chrom_assigned",
            reversed=False, contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id="chr2A", ref_gene_proportion=0.5, contig_len=10_000_000,
        ),
    ]
    idents = {("ctg1", "chr1A"): 0.65, ("ctg2", "chr2A"): 0.68}

    infer_query_subgenomes(clfs, idents)

    assert clfs[0].query_subgenome is None
    assert clfs[0].query_subgenome_grp == 1
    assert clfs[1].query_subgenome is None
    assert clfs[1].query_subgenome_grp == 1


def test_infer_query_subgenomes_two_subgenomes():
    """Two copies per chromosome with bimodal identities -> subgenome split."""
    from final_finalizer.classification.classifier import infer_query_subgenomes
    from final_finalizer.models import ContigClassification

    # Simulate an allotetraploid: each ref chromosome has 2 contigs,
    # one at ~0.69 identity (primary) and one at ~0.60 (secondary).
    clfs = []
    idents = {}
    for i in range(1, 11):
        ref_id = f"chr{i}A"
        # Primary (higher identity)
        prim = ContigClassification(
            original_name=f"ctg_{i}_high", new_name="",
            classification="chrom_assigned", reversed=False,
            contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=ref_id, ref_gene_proportion=0.5,
            contig_len=10_000_000,
        )
        # Secondary (lower identity)
        sec = ContigClassification(
            original_name=f"ctg_{i}_low", new_name="",
            classification="chrom_assigned", reversed=False,
            contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=ref_id, ref_gene_proportion=0.4,
            contig_len=9_000_000,
        )
        clfs.extend([prim, sec])
        idents[(f"ctg_{i}_high", ref_id)] = 0.68 + (i % 3) * 0.01
        idents[(f"ctg_{i}_low", ref_id)] = 0.59 + (i % 3) * 0.01

    infer_query_subgenomes(clfs, idents)

    # Check that contigs were split into two subgenomes
    for i in range(1, 11):
        high_clf = next(c for c in clfs if c.original_name == f"ctg_{i}_high")
        low_clf = next(c for c in clfs if c.original_name == f"ctg_{i}_low")
        assert high_clf.query_subgenome_grp != low_clf.query_subgenome_grp, (
            f"chr{i}A: high and low should be in different subgenome groups"
        )
        # Primary (higher identity) should have no suffix
        assert high_clf.query_subgenome is None
        # Secondary should have a letter suffix
        assert low_clf.query_subgenome is not None


def test_infer_query_subgenomes_no_split_similar_identities():
    """Two copies with very similar identities -> no subgenome split."""
    from final_finalizer.classification.classifier import infer_query_subgenomes
    from final_finalizer.models import ContigClassification

    clfs = []
    idents = {}
    for i in range(1, 11):
        ref_id = f"chr{i}A"
        c1 = ContigClassification(
            original_name=f"ctg_{i}_a", new_name="",
            classification="chrom_assigned", reversed=False,
            contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=ref_id, ref_gene_proportion=0.5,
            contig_len=10_000_000,
        )
        c2 = ContigClassification(
            original_name=f"ctg_{i}_b", new_name="",
            classification="chrom_assigned", reversed=False,
            contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=ref_id, ref_gene_proportion=0.5,
            contig_len=10_000_000,
        )
        clfs.extend([c1, c2])
        # Very similar identities — should NOT split
        idents[(f"ctg_{i}_a", ref_id)] = 0.65 + (i % 5) * 0.005
        idents[(f"ctg_{i}_b", ref_id)] = 0.64 + (i % 5) * 0.005

    infer_query_subgenomes(clfs, idents)

    # All should be in the same group (no subgenome split)
    for clf in clfs:
        assert clf.query_subgenome is None, (
            f"{clf.original_name}: expected no subgenome suffix, "
            f"got {clf.query_subgenome}"
        )


def test_infer_query_subgenomes_sets_identity():
    """Verify seq_identity_vs_ref is set for all contigs."""
    from final_finalizer.classification.classifier import infer_query_subgenomes
    from final_finalizer.models import ContigClassification

    clf = ContigClassification(
        original_name="ctg1", new_name="", classification="chrom_assigned",
        reversed=False, contaminant_taxid=None, contaminant_sci=None,
        assigned_ref_id="chr1A", ref_gene_proportion=0.5, contig_len=10_000_000,
    )
    idents = {("ctg1", "chr1A"): 0.72}

    infer_query_subgenomes([clf], idents)

    assert clf.seq_identity_vs_ref == 0.72


def test_infer_query_subgenomes_mixed_single_and_multi():
    """Mixture of single-copy and multi-copy chromosomes."""
    from final_finalizer.classification.classifier import infer_query_subgenomes
    from final_finalizer.models import ContigClassification

    clfs = []
    idents = {}
    # 8 chromosomes with 2 copies (bimodal)
    for i in range(1, 9):
        ref_id = f"chr{i}A"
        for label, id_offset in [("high", 0.10), ("low", 0.0)]:
            clf = ContigClassification(
                original_name=f"ctg_{i}_{label}", new_name="",
                classification="chrom_assigned", reversed=False,
                contaminant_taxid=None, contaminant_sci=None,
                assigned_ref_id=ref_id, ref_gene_proportion=0.5,
                contig_len=10_000_000,
            )
            clfs.append(clf)
            idents[(f"ctg_{i}_{label}", ref_id)] = 0.59 + id_offset

    # 3 chromosomes with single copy
    for i in range(9, 12):
        ref_id = f"chr{i}A"
        clf = ContigClassification(
            original_name=f"ctg_{i}_only", new_name="",
            classification="chrom_assigned", reversed=False,
            contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=ref_id, ref_gene_proportion=0.5,
            contig_len=10_000_000,
        )
        clfs.append(clf)
        idents[(f"ctg_{i}_only", ref_id)] = 0.65

    infer_query_subgenomes(clfs, idents)

    # Multi-copy chromosomes should be split
    for i in range(1, 9):
        high = next(c for c in clfs if c.original_name == f"ctg_{i}_high")
        low = next(c for c in clfs if c.original_name == f"ctg_{i}_low")
        assert high.query_subgenome_grp != low.query_subgenome_grp, (
            f"chr{i}A: expected different groups"
        )

    # Single-copy chromosomes should have group 1, no suffix
    for i in range(9, 12):
        only = next(c for c in clfs if c.original_name == f"ctg_{i}_only")
        assert only.query_subgenome is None
        assert only.query_subgenome_grp == 1


def test_infer_query_subgenomes_allotriploid():
    """Allotriploid: 2 query subgenomes map to ref subgenome A, 1 maps to P.

    The A channel should detect k=2 (bimodal identities) while the P channel
    should stay at k=1 (single-copy per chromosome).
    """
    from final_finalizer.classification.classifier import infer_query_subgenomes
    from final_finalizer.models import ContigClassification

    clfs = []
    idents = {}

    # 8 chromosomes in ref subgenome A, each with 2 query copies (bimodal)
    for i in range(1, 9):
        ref_id = f"chr{i}A"
        for label, id_offset in [("high", 0.10), ("low", 0.0)]:
            clf = ContigClassification(
                original_name=f"ctg_A{i}_{label}", new_name="",
                classification="chrom_assigned", reversed=False,
                contaminant_taxid=None, contaminant_sci=None,
                assigned_ref_id=ref_id, ref_gene_proportion=0.5,
                contig_len=10_000_000,
            )
            clfs.append(clf)
            idents[(f"ctg_A{i}_{label}", ref_id)] = 0.59 + id_offset

    # 6 chromosomes in ref subgenome P, each with 1 query copy
    for i in range(1, 7):
        ref_id = f"chr{i}P"
        clf = ContigClassification(
            original_name=f"ctg_P{i}", new_name="",
            classification="chrom_assigned", reversed=False,
            contaminant_taxid=None, contaminant_sci=None,
            assigned_ref_id=ref_id, ref_gene_proportion=0.5,
            contig_len=10_000_000,
        )
        clfs.append(clf)
        idents[(f"ctg_P{i}", ref_id)] = 0.75 + (i % 3) * 0.01

    infer_query_subgenomes(clfs, idents)

    # A-subgenome chromosomes should be split into two query subgenomes
    for i in range(1, 9):
        high = next(c for c in clfs if c.original_name == f"ctg_A{i}_high")
        low = next(c for c in clfs if c.original_name == f"ctg_A{i}_low")
        assert high.query_subgenome_grp != low.query_subgenome_grp, (
            f"chr{i}A: expected different subgenome groups"
        )
        assert high.query_subgenome is None  # primary
        assert low.query_subgenome is not None  # secondary (e.g. "B")

    # P-subgenome chromosomes should NOT be split (single copy each)
    for i in range(1, 7):
        p_clf = next(c for c in clfs if c.original_name == f"ctg_P{i}")
        assert p_clf.query_subgenome is None, (
            f"chr{i}P: expected no subgenome suffix, got {p_clf.query_subgenome}"
        )
        assert p_clf.query_subgenome_grp == 1
