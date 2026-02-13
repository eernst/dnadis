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
    span, aligned = result[("contig1", "chr1A")]
    assert span == 600000
    assert aligned == 330000


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
    span, aligned = result[("contig1", "chr12A")]
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
    span, aligned = result[("contig1", "chr12A")]
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
    span, aligned = result[("contig1", "chr12A")]
    assert span == 500000
    assert aligned == 370000  # 250000 + 120000


def test_detect_rearrangement_candidates_no_offtarget():
    """Test that no rearrangement is detected when off-target is below threshold."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    # (span, aligned_bp) tuples - need good density to pass
    qr_cluster_metrics = {
        ("contig1", "chr1A"): (900000, 700000),  # 90% span, 78% density
        ("contig1", "chr2A"): (50000, 40000),    # 5% span (below 10% threshold)
    }
    contig_refs = {"chr1A", "chr2A"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
    )

    assert result is None  # No off-target >= 10%


def test_detect_rearrangement_candidates_single_offtarget():
    """Test detection of single off-target chromosome with sufficient density."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    qr_cluster_metrics = {
        ("contig1", "chr1A"): (700000, 500000),  # 70% span
        ("contig1", "chr5B"): (200000, 100000),  # 20% span, 50% density (passes both)
    }
    contig_refs = {"chr1A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
    )

    assert result == "chr5B"


def test_detect_rearrangement_candidates_low_density_rejected():
    """Test that off-target with large span but low density is rejected."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    qr_cluster_metrics = {
        ("contig1", "chr1A"): (700000, 500000),  # 70% span
        ("contig1", "chr5B"): (200000, 20000),   # 20% span but only 10% density (< 15% min)
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
    )

    assert result is None  # Density too low


def test_detect_rearrangement_candidates_multiple_offtarget():
    """Test detection of multiple off-target chromosomes, ordered by span fraction."""
    from final_finalizer.classification.classifier import detect_rearrangement_candidates

    qr_cluster_metrics = {
        ("contig1", "chr1A"): (500000, 400000),  # 50% span
        ("contig1", "chr2A"): (150000, 50000),   # 15% span, 33% density
        ("contig1", "chr5B"): (250000, 100000),  # 25% span, 40% density (highest span)
    }
    contig_refs = {"chr1A", "chr2A", "chr5B"}

    result = detect_rearrangement_candidates(
        contig="contig1",
        assigned_ref_id="chr1A",
        contig_len=1000000,
        qr_cluster_metrics=qr_cluster_metrics,
        contig_refs=contig_refs,
        threshold=0.10,
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
