"""Checks on the R Markdown report templates and their shared R helpers."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPORTS_DIR = Path(__file__).resolve().parent.parent / "dnadis" / "output" / "reports"
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
