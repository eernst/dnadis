"""Checks on the R Markdown report templates and their shared R helpers."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_DIR))  # refresh_reports.py lives at the repo root

REPORTS_DIR = REPO_DIR / "dnadis" / "output" / "reports"
TEMPLATES = sorted(REPORTS_DIR.glob("*.tmpl.Rmd"))

_CHUNK_RE = re.compile(r"^```\{r\s*(?P<header>[^}]*)\}\s*$(?P<body>.*?)^```\s*$",
                       re.MULTILINE | re.DOTALL)
_FIG_OPT_RE = re.compile(r"\b(?:fig\.(?:height|width)|out\.(?:height|width))\s*=\s*"
                         r"(?P<expr>[^,()]*(?:\([^()]*\))?[^,]*)")
_IDENT_RE = re.compile(r"(?<![\w.])([A-Za-z][\w.]*)\b(?!\s*\()")


def _chunks(text: str) -> list[tuple[str, str]]:
    return [(m.group("header"), m.group("body")) for m in _CHUNK_RE.finditer(text)]


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_fig_option_variables_defined_in_unconditional_chunks(template):
    """knitr evaluates fig.* chunk options even when the chunk's eval is FALSE,
    so a variable used there must be assigned in a chunk that always runs."""
    chunks = _chunks(template.read_text())
    assert chunks, f"no R chunks parsed from {template.name}"

    always_assigned: set[str] = set()
    for header, body in chunks:
        if re.search(r"\beval\s*=", header):
            continue
        always_assigned.update(re.findall(r"^\s*([A-Za-z][\w.]*)\s*<-", body, re.MULTILINE))

    missing = []
    for header, _ in chunks:
        for m in _FIG_OPT_RE.finditer(header):
            for name in _IDENT_RE.findall(m.group("expr")):
                if name in {"TRUE", "FALSE", "NA", "NULL"}:
                    continue
                if name not in always_assigned:
                    label = header.split(",")[0].strip()
                    missing.append(f"{label}: {name}")
    assert not missing, (
        "chunk options use variables assigned only in conditional chunks: "
        + ", ".join(missing)
    )


def _r_packages_available() -> bool:
    if shutil.which("Rscript") is None:
        return False
    rc = subprocess.run(
        ["Rscript", "-e",
         "for (p in c('pacman', 'readr', 'dplyr')) "
         "if (!requireNamespace(p, quietly = TRUE)) quit(status = 1)"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    ).returncode
    return rc == 0


@pytest.mark.skipif(not _r_packages_available(), reason="Rscript with readr/dplyr not available")
def test_read_dnadis_tsv_keeps_digit_only_names_as_character(tmp_path):
    """Name columns are character in every file, so per-assembly tables bind
    even when one file has only digit-only names (reference chromosomes
    named "1".."10") and another has "1_f1"-style names."""
    a = tmp_path / "a.contig_summary.tsv"
    b = tmp_path / "b.contig_summary.tsv"
    a.write_text(
        "contig\tassigned_ref_id\tassigned_subgenome\tlength\tredundancy_target\n"
        "1\t1\tNA\t1000\t9\n"
        "ptg000002l\tNA\tNA\t200\t3\n"
    )
    b.write_text(
        "contig\tassigned_ref_id\tassigned_subgenome\tlength\tredundancy_target\n"
        "1_f1\t1\tNA\t800\t\n"
        "ptg000005l\tNA\tNA\t300\t1_f1\n"
    )
    script = f"""
suppressMessages(source("{REPORTS_DIR / 'common.R'}"))
a <- read_dnadis_tsv("{a}")
b <- read_dnadis_tsv("{b}",
                     col_types = cols(assigned_subgenome = col_character(),
                                      length = col_integer()))
stopifnot(is.character(a$contig), is.character(a$assigned_ref_id),
          is.character(a$redundancy_target), is.character(b$redundancy_target),
          is.numeric(a$length), is.integer(b$length))
df <- dplyr::bind_rows(a, b)
stopifnot(identical(df$redundancy_target, c("9", "3", NA, "1_f1")),
          identical(df$contig, c("1", "ptg000002l", "1_f1", "ptg000005l")))
cat("OK\\n")
"""
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True,
                          timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "OK" in proc.stdout


def _contiguity_cases():
    """(ref_id, [(contig, ref_copy, has_5p, has_3p, n_gaps), ...]) for each state rule."""
    return [
        ("chr1", [("t2t", 1, True, True, 0)]),
        ("chr2", [("one", 1, True, False, 0)]),
        ("chr3", [("p", 1, True, False, 0), ("q", 1, False, True, 0)]),
        ("chr4", [("h1", 1, True, True, 0), ("h2", 2, True, True, 0)]),
        ("chr5", [("h1b", 1, True, True, 0), ("h2a", 2, False, False, 0), ("h2b", 2, True, False, 0)]),
        ("chr6", [("t2t6", 1, True, True, 0), ("unplaced", None, False, False, 0)]),
        ("chr7", [("u1", None, False, False, 0), ("u2", None, True, True, 0)]),
        # Scaffold gaps: telomere-to-telomere with gaps is "gapped"; a gapped
        # contig without both telomeres stays "single"; the least contiguous
        # copy sets the state.
        ("chr8", [("gapped", 1, True, True, 3)]),
        ("chr9", [("gap_single", 1, True, False, 2)]),
        ("chr10", [("t2t10", 1, True, True, 0), ("gapped10", 2, True, True, 1)]),
        ("chr11", [("gapped11", 1, True, True, 1), ("one11", 2, False, False, 0)]),
    ]


@pytest.mark.skipif(not _r_packages_available(), reason="Rscript with readr/dplyr not available")
def test_contiguity_state_python_and_r_agree(tmp_path):
    """comparison.py and common.R compute the same per-chromosome state,
    copy count and telomere flag from the same contigs."""
    from dnadis.models import ContigClassification
    from dnadis.output.comparison import _chrom_contiguity_state

    rows = ["contig\toriginal_name\tclassification\tassigned_ref_id\tlength\t"
            "ref_copy\thas_5p_telomere\thas_3p_telomere\tn_gaps"]
    expected = {}
    for ref_id, contigs in _contiguity_cases():
        clfs = []
        for name, copy, t5, t3, gaps in contigs:
            name = f"{ref_id}_{name}"
            c = ContigClassification(
                original_name=name, new_name=name, classification="chrom_assigned",
                reversed=False, cobiont_taxid=None, cobiont_sci=None,
                assigned_ref_id=ref_id, ref_gene_proportion=None, contig_len=1_000_000,
            )
            c.ref_copy, c.has_5p_telomere, c.has_3p_telomere = copy, t5, t3
            c.n_gaps = gaps
            clfs.append(c)
            yn = lambda v: "yes" if v else "no"
            rows.append(f"{name}\t{name}\tchrom_assigned\t{ref_id}\t1000000\t"
                        f"{'' if copy is None else copy}\t{yn(t5)}\t{yn(t3)}\t{gaps}")
        state, n_copies = _chrom_contiguity_state(clfs, 0)
        both = any(c.has_5p_telomere and c.has_3p_telomere for c in clfs)
        expected[ref_id] = f"{state}|{n_copies}|{'TRUE' if both else 'FALSE'}"
    summary = tmp_path / "x.contig_summary.tsv"
    summary.write_text("\n".join(rows) + "\n")
    refs = tmp_path / "ref_lengths.tsv"
    refs.write_text("ref_id\tref_len\n" + "".join(
        f"{r}\t5000000\n" for r, _ in _contiguity_cases()))

    script = f"""
suppressMessages(source("{REPORTS_DIR / 'common.R'}"))
s <- read_dnadis_tsv("{summary}")
r <- read_dnadis_tsv("{refs}")
ct <- compute_chrom_contiguity(s, r)
cat(paste0(ct$ref_id, "=", ct$state, "|", ct$n_copies, "|", ct$both_telo), sep = "\\n")
"""
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True,
                          timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    got = dict(line.split("=", 1) for line in proc.stdout.split())
    assert got == expected
    assert expected["chr8"].startswith("gapped|") and expected["chr9"].startswith("single|")
    assert expected["chr10"].startswith("gapped|") and expected["chr11"].startswith("single|")


@pytest.mark.skipif(not _r_packages_available(), reason="Rscript with readr/dplyr not available")
def test_compute_chrom_contiguity_names_missing_columns(tmp_path):
    """A contig summary from a different dnadis version stops with a message
    naming the missing column and the fix, instead of a silent fallback."""
    summary = tmp_path / "x.contig_summary.tsv"
    summary.write_text("contig\toriginal_name\tclassification\tassigned_ref_id\tlength\t"
                       "has_5p_telomere\thas_3p_telomere\n"
                       "a\ta\tchrom_assigned\tchr1\t10\tyes\tno\n")
    script = f"""
suppressMessages(source("{REPORTS_DIR / 'common.R'}"))
s <- read_dnadis_tsv("{summary}")
compute_chrom_contiguity(s, tibble::tibble(ref_id = "chr1", ref_len = 10))
"""
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True,
                          timeout=300)
    assert proc.returncode != 0
    assert "ref_copy" in proc.stderr and "Rerun dnadis" in proc.stderr


def test_generated_report_inlines_common_r():
    """Generated .Rmd files carry the helpers instead of sourcing common.R."""
    from dnadis.output.plotting import inline_common_r

    for template in TEMPLATES:
        text = template.read_text()
        assert 'source("__COMMON_R__")' in text
        filled = inline_common_r(text)
        assert "__COMMON_R__" not in filled
        assert "compute_chrom_contiguity <- function" in filled


def test_templates_declare_output_schema():
    import refresh_reports

    for template in TEMPLATES:
        assert refresh_reports.output_schema(template.read_text()) >= 1, template.name


def test_assembly_report_requires_columns_dnadis_writes(tmp_path):
    """Columns the assembly report requires from rearrangements.tsv are ones
    write_rearrangements_tsv() writes, even with no calls."""
    from dnadis.output.tsv_output import write_rearrangements_tsv

    text = (REPORTS_DIR / "assembly_report.tmpl.Rmd").read_text()
    calls = re.findall(r'require_columns\(\s*rearr_df\s*,\s*c\(([^)]*)\)', text)
    assert calls, "assembly report does not check rearrangements.tsv columns"
    required = {c for call in calls for c in re.findall(r'"([^"]+)"', call)}

    out = tmp_path / "rearrangements.tsv"
    write_rearrangements_tsv(out, [])
    header = out.read_text().splitlines()[0].split("\t")
    assert required <= set(header), required - set(header)


def test_refresh_skips_reports_with_older_outputs(tmp_path, capsys):
    """A report generated before the template's output schema is left alone."""
    import refresh_reports

    template_dir = tmp_path / "templates"
    template_dir.mkdir()
    (template_dir / "assembly_report.tmpl.Rmd").write_text(
        "params:\n  summary_tsv: \"__SUMMARY__\"\n  output_schema: 2\n")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    old = run_dir / "asm.assembly_report.Rmd"
    old_text = "params:\n  summary_tsv: \"/x/asm.contig_summary.tsv\"\n  output_schema: 1\n"
    old.write_text(old_text)

    result = refresh_reports.refresh_script(
        "assembly_report", old, template_dir, run_dir / "asm")
    assert result is None
    assert old.read_text() == old_text
    assert "Rerun dnadis" in capsys.readouterr().err

    old.write_text(old_text.replace("output_schema: 1", "output_schema: 2"))
    assert refresh_reports.refresh_script(
        "assembly_report", old, template_dir, run_dir / "asm") is not None


def test_refresh_prefers_current_agp_name(tmp_path):
    """With both AGP names in the folder, refresh uses the current one."""
    import refresh_reports

    prefix = tmp_path / "asm" / "asm"
    prefix.parent.mkdir()
    for sfx in (".chrs.scaffolds.agp", ".scaffolded.agp"):
        Path(str(prefix) + sfx).write_text("")
    template = 'params:\n  agp_tsv: "__AGP__"\n'
    values = refresh_reports.build_placeholder_values(
        "assembly_report", prefix, template, template)
    assert values["__AGP__"].endswith("asm.chrs.scaffolds.agp")

    Path(str(prefix) + ".chrs.scaffolds.agp").unlink()
    values = refresh_reports.build_placeholder_values(
        "assembly_report", prefix, template, template)
    assert values["__AGP__"].endswith("asm.scaffolded.agp")
