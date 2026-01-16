import gzip
import shutil
import sys
import urllib.request
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import final_finalizer as ff


TAIR10_FASTA = Path(
    "/home/eernst/reference/phytozome/PhytozomeV13/Athaliana/Araport11/assembly/Athaliana_447_TAIR10.fa.gz"
)
ARAPORT11_GFF3 = Path(
    "/home/eernst/reference/phytozome/PhytozomeV13/Athaliana/Araport11/annotation/Athaliana_447_Araport11.gene.gff3.gz"
)
LER_URL = "https://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/Ler/Ler.chr.all.v2.0.fasta.gz"
LER_FASTA = Path("data/Ler.chr.all.v2.0.fasta.gz")


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
def test_arabidopsis_full_chromosomes():
    if not TAIR10_FASTA.exists():
        pytest.skip(f"missing TAIR10 FASTA: {TAIR10_FASTA}")
    if not ARAPORT11_GFF3.exists():
        pytest.skip(f"missing Araport11 GFF3: {ARAPORT11_GFF3}")

    _download_if_missing(LER_URL, LER_FASTA)

    ref_lengths = ff.read_fasta_lengths(TAIR10_FASTA)
    for chrom in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"):
        assert chrom in ref_lengths
        assert ref_lengths[chrom] > 0

    assert ff.get_min_nuclear_chrom_length(ref_lengths) > 0

    tx2loc, tx2gene = ff.parse_gff3_transcript_coords(ARAPORT11_GFF3)
    assert tx2loc
    first_record = _first_transcript_record(ARAPORT11_GFF3)
    assert first_record is not None
    tid, seqid, start, end, strand, gene_id = first_record
    assert tid in tx2loc
    assert tx2loc[tid] == (seqid, start - 1, end, strand if strand in ("+", "-") else "+")
    assert tx2gene.get(tid) == gene_id

    ler_lengths = ff.read_fasta_lengths(LER_FASTA)
    assert len(ler_lengths) >= 5
    top5 = sorted(ler_lengths.values(), reverse=True)[:5]
    assert all(length > 1_000_000 for length in top5)


@pytest.mark.integration
def test_arabidopsis_normalization_pipeline():
    """Test that normalization correctly handles Chr -> chr conversion."""
    if not TAIR10_FASTA.exists():
        pytest.skip(f"missing TAIR10 FASTA: {TAIR10_FASTA}")
    if not ARAPORT11_GFF3.exists():
        pytest.skip(f"missing Araport11 GFF3: {ARAPORT11_GFF3}")

    # Test read_fasta_lengths_with_map normalization
    ref_lengths_norm, orig_to_norm, norm_to_orig = ff.read_fasta_lengths_with_map(TAIR10_FASTA)

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
    tx2loc_norm, tx2gene = ff.parse_gff3_transcript_coords(ARAPORT11_GFF3, ref_id_map=orig_to_norm)
    first_record = _first_transcript_record(ARAPORT11_GFF3)
    assert first_record is not None
    tid, seqid, start, end, strand, gene_id = first_record

    # With ref_id_map, seqid should be normalized
    assert tid in tx2loc_norm
    normalized_seqid = orig_to_norm.get(seqid, seqid)
    assert tx2loc_norm[tid][0] == normalized_seqid

    # Test split_chrom_subgenome preserves style
    assert ff.split_chrom_subgenome("Chr1") == ("Chr1", "NA")
    assert ff.split_chrom_subgenome("chr1") == ("chr1", "NA")
    assert ff.split_chrom_subgenome("ChrC") == ("ChrC", "NA")
    assert ff.split_chrom_subgenome("chrM") == ("chrM", "NA")

    # Test normalize_organelle_id with various conventions
    assert ff.normalize_organelle_id("ChrC") == "chrC"
    assert ff.normalize_organelle_id("ChrM") == "chrM"
    assert ff.normalize_organelle_id("Mt") == "chrM"
    assert ff.normalize_organelle_id("Pt") == "chrC"
    assert ff.normalize_organelle_id("Chr1") is None  # Not an organelle
