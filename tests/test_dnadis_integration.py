import gzip
import os
import shutil
import sys
import urllib.request
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import dnadis


def _path_from_env(var: str) -> Path | None:
    """Return a Path from the environment variable, or None if unset/empty."""
    val = os.environ.get(var)
    return Path(val) if val else None


# Integration tests requiring Arabidopsis TAIR10/Araport11 files read their paths
# from environment variables so they work on any machine. If unset, the tests
# using these fixtures are skipped. Typical download source: Phytozome
# (https://phytozome-next.jgi.doe.gov/) — requires login.
TAIR10_FASTA = _path_from_env("DNADIS_TAIR10_FASTA")
ARAPORT11_GFF3 = _path_from_env("DNADIS_ARAPORT11_GFF3")

_TAIR10_SKIP_MSG = (
    "Set DNADIS_TAIR10_FASTA to a TAIR10 Arabidopsis FASTA to run this test"
)
_ARAPORT11_SKIP_MSG = (
    "Set DNADIS_ARAPORT11_GFF3 to an Araport11 Arabidopsis GFF3 to run this test"
)
LER_URL = "https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.chr.all.v2.0.fasta.gz"
LER_FASTA = Path("tests/data/Ler.chr.all.v2.0.fasta.gz")
CMV_FASTA = Path("tests/data/GCF_000864745.1_ViralMultiSegProj15470_genomic.fna.gz")


def _download_if_missing(url: str, dest: Path) -> None:
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url) as resp, dest.open("wb") as out:
            shutil.copyfileobj(resp, out)
    except Exception as exc:
        pytest.skip(f"failed to download {url}: {exc}")


def _first_transcript_record(gff3_path: Path):
    opener = gzip.open if gff3_path.suffix in (".gz", ".bgz") else open
    with opener(gff3_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = parts
            if ftype.lower() not in ("mrna", "transcript"):
                continue
            attrs_dict = {}
            for kv in attrs.split(";"):
                kv = kv.strip()
                if not kv or "=" not in kv:
                    continue
                key, value = kv.split("=", 1)
                attrs_dict[key] = value
            tid = attrs_dict.get("ID")
            parent = attrs_dict.get("Parent", "")
            gene_id = parent.split(",")[0] if parent else tid
            if not tid:
                continue
            return tid, seqid, int(start), int(end), strand, gene_id
    return None


@pytest.mark.integration
def test_arabidopsis_full_chromosomes(tmp_path):
    if TAIR10_FASTA is None or not TAIR10_FASTA.exists():
        pytest.skip(_TAIR10_SKIP_MSG)
    if ARAPORT11_GFF3 is None or not ARAPORT11_GFF3.exists():
        pytest.skip(_ARAPORT11_SKIP_MSG)

    _download_if_missing(LER_URL, LER_FASTA)
    if not CMV_FASTA.exists():
        pytest.skip(f"missing CMV FASTA: {CMV_FASTA}")

    ref_lengths = dnadis.read_fasta_lengths(TAIR10_FASTA)
    for chrom in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"):
        assert chrom in ref_lengths
        assert ref_lengths[chrom] > 0

    assert dnadis.get_min_nuclear_chrom_length(ref_lengths) > 0

    tx2loc, tx2gene = dnadis.parse_gff3_transcript_coords(ARAPORT11_GFF3)
    assert tx2loc
    first_record = _first_transcript_record(ARAPORT11_GFF3)
    assert first_record is not None
    tid, seqid, start, end, strand, gene_id = first_record
    assert tid in tx2loc
    assert tx2loc[tid] == (seqid, start - 1, end, strand if strand in ("+", "-") else "+")
    assert tx2gene.get(tid) == gene_id

    ler_lengths = dnadis.read_fasta_lengths(LER_FASTA)
    assert len(ler_lengths) >= 5
    top5 = sorted(ler_lengths.values(), reverse=True)[:5]
    assert all(length > 1_000_000 for length in top5)

    cmv_lengths = dnadis.read_fasta_lengths(CMV_FASTA)
    assert cmv_lengths
    assert all(length > 0 for length in cmv_lengths.values())

    combined_fasta = tmp_path / "ler_plus_cmv.fasta.gz"
    with gzip.open(combined_fasta, "wb") as out_fh:
        with gzip.open(LER_FASTA, "rb") as in_fh:
            shutil.copyfileobj(in_fh, out_fh)
        with gzip.open(CMV_FASTA, "rb") as in_fh:
            shutil.copyfileobj(in_fh, out_fh)

    combined_lengths = dnadis.read_fasta_lengths(combined_fasta)
    assert any(name in combined_lengths for name in cmv_lengths.keys())


@pytest.mark.integration
def test_arabidopsis_normalization_pipeline():
    """Test that normalization correctly handles Chr -> chr conversion."""
    if TAIR10_FASTA is None or not TAIR10_FASTA.exists():
        pytest.skip(_TAIR10_SKIP_MSG)
    if ARAPORT11_GFF3 is None or not ARAPORT11_GFF3.exists():
        pytest.skip(_ARAPORT11_SKIP_MSG)

    # Test read_fasta_lengths_with_map normalization
    ref_lengths_norm, orig_to_norm, norm_to_orig = dnadis.read_fasta_lengths_with_map(TAIR10_FASTA)

    # Original IDs use "Chr" prefix (Phytozome convention)
    assert "Chr1" in orig_to_norm
    assert orig_to_norm["Chr1"] == "chr1"
    assert orig_to_norm["ChrC"] == "chrC"
    assert orig_to_norm["ChrM"] == "chrM"

    # Normalized lengths use "chr" prefix
    assert "chr1" in ref_lengths_norm
    assert "chrC" in ref_lengths_norm
    assert "chrM" in ref_lengths_norm

    # Reverse mapping preserves original style
    assert norm_to_orig["chr1"] == "Chr1"
    assert norm_to_orig["chrC"] == "ChrC"
    assert norm_to_orig["chrM"] == "ChrM"

    # Test parse_gff3_transcript_coords with normalization
    tx2loc_norm, tx2gene = dnadis.parse_gff3_transcript_coords(ARAPORT11_GFF3, ref_id_map=orig_to_norm)
    first_record = _first_transcript_record(ARAPORT11_GFF3)
    assert first_record is not None
    tid, seqid, start, end, strand, gene_id = first_record

    # With ref_id_map, seqid should be normalized
    assert tid in tx2loc_norm
    normalized_seqid = orig_to_norm.get(seqid, seqid)
    assert tx2loc_norm[tid][0] == normalized_seqid

    # Test split_chrom_subgenome preserves style
    assert dnadis.split_chrom_subgenome("Chr1") == ("Chr1", "NA")
    assert dnadis.split_chrom_subgenome("chr1") == ("chr1", "NA")
    assert dnadis.split_chrom_subgenome("ChrC") == ("ChrC", "NA")
    assert dnadis.split_chrom_subgenome("chrM") == ("chrM", "NA")

    # Test normalize_organelle_id with various conventions
    assert dnadis.normalize_organelle_id("ChrC") == "chrC"
    assert dnadis.normalize_organelle_id("ChrM") == "chrM"
    assert dnadis.normalize_organelle_id("Mt") == "chrM"
    assert dnadis.normalize_organelle_id("Pt") == "chrC"
    assert dnadis.normalize_organelle_id("Chr1") is None  # Not an organelle


@pytest.mark.integration
def test_contig_orientation_detection(tmp_path):
    """Test that contig orientation detection correctly identifies reverse-complemented contigs.

    This test creates a synthetic query assembly with:
    1. A forward-oriented chromosome (copied from reference)
    2. A reverse-complemented chromosome (reverse complement of reference)
    3. Runs the full pipeline to verify orientation detection
    4. Checks that the reversed field is correctly set in the output

    The test helps verify that the orientation detection logic works correctly,
    which is critical for proper chromosome assembly orientation.

    CURRENT STATUS: This test is marked as xfail because there is a known bug
    where miniprot does not report the correct strand orientation for reverse-
    complemented contigs. All alignments are reported as "+" strand even when
    the contig is reverse-complemented relative to the reference.
    """
    if TAIR10_FASTA is None or not TAIR10_FASTA.exists():
        pytest.skip(_TAIR10_SKIP_MSG)
    if ARAPORT11_GFF3 is None or not ARAPORT11_GFF3.exists():
        pytest.skip(_ARAPORT11_SKIP_MSG)

    # Read reference sequences - we'll use Chr4 and Chr5 for testing
    # These are smaller chromosomes (~18-27 Mbp) which makes the test faster
    ref_seqs = dnadis.read_fasta_sequences(TAIR10_FASTA)
    assert "Chr4" in ref_seqs, "Chr4 not found in TAIR10 reference"
    assert "Chr5" in ref_seqs, "Chr5 not found in TAIR10 reference"

    chr4_seq = ref_seqs["Chr4"]
    chr5_seq = ref_seqs["Chr5"]

    # Create a subset reference FASTA with only Chr4 and Chr5
    # This works around gffread's requirement for consistent line lengths
    subset_ref_fasta = tmp_path / "subset_ref.fasta"
    subset_ref_seqs = {
        "Chr4": chr4_seq,
        "Chr5": chr5_seq,
    }
    dnadis.write_fasta(subset_ref_seqs, subset_ref_fasta)

    # Create a subset GFF3 with only Chr4 and Chr5 annotations
    subset_gff3 = tmp_path / "subset_ref.gff3"
    with gzip.open(ARAPORT11_GFF3, "rt") as in_f, subset_gff3.open("w") as out_f:
        for line in in_f:
            if line.startswith("#"):
                out_f.write(line)
            elif line.startswith("Chr4\t") or line.startswith("Chr5\t"):
                out_f.write(line)

    # Create synthetic query assembly:
    # - contig_1: Copy of Chr4 in forward orientation (should not be reversed)
    # - contig_2: Reverse complement of Chr5 (should be detected and reversed)
    chr5_revcomp = dnadis.reverse_complement(chr5_seq)

    query_fasta = tmp_path / "test_query.fasta"
    query_seqs = {
        "contig_forward": chr4_seq,
        "contig_reversed": chr5_revcomp,
    }
    dnadis.write_fasta(query_seqs, query_fasta)

    # Verify the query assembly was created correctly
    assert query_fasta.exists()
    query_lengths = dnadis.read_fasta_lengths(query_fasta)
    assert "contig_forward" in query_lengths
    assert "contig_reversed" in query_lengths
    assert query_lengths["contig_forward"] == len(chr4_seq)
    assert query_lengths["contig_reversed"] == len(chr5_revcomp)

    # Create output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)

    # Assembly name is derived from query filename: "test_query"
    asm_name = "test_query"
    output_prefix = output_dir / asm_name / asm_name

    # Run the full pipeline
    # We use a simplified set of arguments focusing on the core alignment and orientation detection
    import subprocess
    import sys

    cmd = [
        sys.executable, "-m", "dnadis",
        "--query", str(query_fasta),
        "--ref", str(subset_ref_fasta),
        "--ref-gff3", str(subset_gff3),
        "--output-dir", str(output_dir),
        "-t", "4",
        "--skip-contaminants",  # Skip contaminant detection to speed up test
        "--skip-rdna",  # Skip rDNA detection to speed up test
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=600,  # 10 minute timeout
    )

    # Check that the pipeline succeeded
    assert result.returncode == 0, f"Pipeline failed with stderr:\n{result.stderr}"

    # Verify that output files were created
    contig_summary = output_prefix.with_suffix(".contig_summary.tsv")
    assert contig_summary.exists(), "contig_summary.tsv not created"

    macro_blocks = output_prefix.with_suffix(".macro_blocks.tsv")
    assert macro_blocks.exists(), "macro_blocks.tsv not created"

    chrs_fasta = Path(f"{output_prefix}.chrs.fasta")
    assert chrs_fasta.exists(), "chrs.fasta not created"

    # Parse the contig_summary.tsv to check orientation detection
    classifications = {}
    with contig_summary.open("r") as f:
        header = f.readline().strip().split("\t")
        # Find column indices
        contig_idx = header.index("contig")
        original_name_idx = header.index("original_name")
        classification_idx = header.index("classification")
        reversed_idx = header.index("reversed")
        assigned_ref_idx = header.index("assigned_ref_id")

        for line in f:
            fields = line.strip().split("\t")
            original_name = fields[original_name_idx]
            classifications[original_name] = {
                "new_name": fields[contig_idx],
                "classification": fields[classification_idx],
                "reversed": fields[reversed_idx],
                "assigned_ref_id": fields[assigned_ref_idx],
            }

    # Verify both contigs are present and classified as chromosomes
    assert "contig_forward" in classifications, "contig_forward not found in output"
    assert "contig_reversed" in classifications, "contig_reversed not found in output"

    # Check classifications
    fwd_class = classifications["contig_forward"]
    rev_class = classifications["contig_reversed"]

    assert fwd_class["classification"] == "chrom_assigned", \
        f"contig_forward not classified as chromosome: {fwd_class['classification']}"
    assert rev_class["classification"] == "chrom_assigned", \
        f"contig_reversed not classified as chromosome: {rev_class['classification']}"

    # Check reference assignments
    # Normalized reference IDs use lowercase (chr4, chr5)
    assert fwd_class["assigned_ref_id"] == "chr4" or fwd_class["assigned_ref_id"] == "Chr4", \
        f"contig_forward assigned to wrong reference: {fwd_class['assigned_ref_id']}"
    assert rev_class["assigned_ref_id"] == "chr5" or rev_class["assigned_ref_id"] == "Chr5", \
        f"contig_reversed assigned to wrong reference: {rev_class['assigned_ref_id']}"

    # CRITICAL TEST: Check orientation detection
    # contig_forward should NOT be reversed (it's already in correct orientation)
    assert fwd_class["reversed"] == "no", \
        f"contig_forward incorrectly marked for reversal: {fwd_class['reversed']}"

    # contig_reversed SHOULD be reversed (it's reverse-complemented)
    # NOTE: This assertion will fail until the orientation detection bug is fixed.
    # The bug is that miniprot aligns to reverse-complemented sequences but reports
    # strand as "+" relative to the query sequence, not relative to the reference.
    # This causes determine_contig_orientations() to see all "+" strands and not
    # detect that the contig needs reversal.
    #
    # Expected behavior after fix:
    # - Macro blocks should show "-" strand for reverse-complemented contigs
    # - determine_contig_orientations() should detect rev > fwd and set reversed=True
    # - Output FASTA should be reverse-complemented back to match reference orientation
    assert rev_class["reversed"] == "yes", \
        f"contig_reversed NOT marked for reversal (this is the bug!): {rev_class['reversed']}\n" \
        f"Check macro_blocks.tsv - all strands are likely '+' even for the reverse-complemented contig.\n" \
        f"This indicates the alignment step is not detecting relative orientation correctly."

    # Verify the output sequences are correctly oriented
    output_seqs = dnadis.read_fasta_sequences(chrs_fasta)

    # The output should have renamed contigs (e.g., chr4, chr5 or Chr4, Chr5)
    # Find the contigs by their content
    fwd_output_name = fwd_class["new_name"]
    rev_output_name = rev_class["new_name"]

    assert fwd_output_name in output_seqs, f"{fwd_output_name} not in output FASTA"
    assert rev_output_name in output_seqs, f"{rev_output_name} not in output FASTA"

    # Check that contig_forward sequence matches reference chr4 (forward orientation)
    fwd_output_seq = output_seqs[fwd_output_name]
    assert fwd_output_seq == chr4_seq, \
        "contig_forward sequence doesn't match reference (should be unchanged)"

    # Check that contig_reversed was reverse-complemented back to match reference chr5
    rev_output_seq = output_seqs[rev_output_name]
    assert rev_output_seq == chr5_seq, \
        "contig_reversed was not properly reverse-complemented (should match reference chr5)"

    # Additional verification: Check macro_blocks.tsv for strand information
    forward_strands = 0
    reverse_strands = 0
    with macro_blocks.open("r") as f:
        header = f.readline().strip().split("\t")
        contig_idx = header.index("contig")
        strand_idx = header.index("strand")

        for line in f:
            fields = line.strip().split("\t")
            contig = fields[contig_idx]
            strand = fields[strand_idx]

            if contig == "contig_forward":
                if strand == "+":
                    forward_strands += 1
            elif contig == "contig_reversed":
                if strand == "-":
                    reverse_strands += 1

    # Verify that we detected the expected strand patterns
    assert forward_strands > 0, \
        "No forward (+) strands detected for contig_forward in macro_blocks"
    assert reverse_strands > 0, \
        "No reverse (-) strands detected for contig_reversed in macro_blocks"
