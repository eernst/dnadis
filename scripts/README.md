# Developer scripts

These scripts support development and CI, not the main `dnadis` pipeline.

- **`update_readme_conda_versions.py`** — invoked by the scheduled CI workflow
  (`.github/workflows/ci-scheduled.yml`) to refresh the "Latest tested conda
  package versions" block in `README.md` based on a `conda list` snapshot.
- **`update_readme_python.py`** — companion to the above; updates the latest
  tested Python version reported in `README.md`.
- **`bundle_single_file.py`** — experimental utility that concatenates the
  `dnadis` package into a single script for environments that cannot install a
  proper Python package.
- **`bundle_reports.py`** — end-user utility that packs the HTML reports from a
  dnadis output folder into a compressed tarball, preserving the per-assembly
  subdirectories so the comparison report's links still resolve. Point it at an
  output directory: `python scripts/bundle_reports.py OUTPUT_DIR`.

The top-level `refresh_reports.py` is a separate developer utility for
re-rendering HTML reports from updated R/Rmd templates; see its `--help` output
and `docs/output_formats.md` for usage.
