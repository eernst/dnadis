# Bundled reference data

This directory ships small reference files that `dnadis` uses by default when
corresponding user-supplied inputs are not provided.

## 45S rDNA reference (Arabidopsis thaliana)

- **`athal-45s-ref.fa`** — the canonical 45S ribosomal DNA repeat from *Arabidopsis
  thaliana* (Col-0), spanning `CP002686.1:14195483-14204860` (9378 bp). Used as
  the default rDNA reference by `dnadis/detection/rdna.py` when the user does
  not pass `--rdna-ref`.

- **`athal-45s-features.tsv`** — sub-feature coordinates (18S, ITS1, 5.8S, ITS2,
  25S) for the reference above. 1-based, inclusive. Coordinates are based on
  Gruendler et al. 1989; Unfried et al. 1989; Unfried and Gruendler 1990;
  Cokus et al. 2008.

For non-Arabidopsis assemblies, supply a species-specific rDNA reference via
`--rdna-ref`. The bundled file is a sensible default for eukaryotic rDNA
detection because the 45S locus is highly conserved, but a closer reference
improves sensitivity.

Additional bundled reference data lives under `dnadis/data/rfam/` (eukaryotic
rRNA covariance models from Rfam 15.0, used by the rDNA consensus
sub-feature annotation step).
