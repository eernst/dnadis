# Rfam rRNA Covariance Models

This directory contains a minimal subset of Rfam covariance models for
annotating eukaryotic ribosomal RNA genes.

## Source

Models extracted from Rfam 15.0 (September 2024):
- https://rfam.org/
- https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/

## Included Models

| Accession | Name | Description | Length |
|-----------|------|-------------|--------|
| RF00001 | 5S_rRNA | 5S ribosomal RNA | 120 bp |
| RF00002 | 5_8S_rRNA | 5.8S ribosomal RNA | 154 bp |
| RF01960 | SSU_rRNA_eukarya | 18S ribosomal RNA | 1831 bp |
| RF02543 | LSU_rRNA_eukarya | 28S/25S ribosomal RNA | 3401 bp |

## Usage

These models are used by `final_finalizer` to annotate rRNA sub-features
on the consensus 45S sequence with accurate, structure-based boundaries.

Requires Infernal (available via conda: `conda install -c bioconda infernal`).

## License

Rfam is freely available under the Creative Commons Zero (CC0) license.
See: https://rfam.org/about#license

## Citation

If you use these models, please cite:

Kalvari I, Nawrocki EP, Ontiveros-Palacios N, et al. (2021)
Rfam 14: expanded coverage of metagenomic, viral and microRNA families.
Nucleic Acids Research, 49(D1):D192-D200.
https://doi.org/10.1093/nar/gkaa1047
