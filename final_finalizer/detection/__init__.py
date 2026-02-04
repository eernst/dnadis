# final_finalizer/detection/__init__.py
"""Detection modules for organelles, rDNA, debris, and contaminants."""

from final_finalizer.detection.blast import (
    parse_blast_coverage,
    run_blastn_megablast,
    run_makeblastdb,
)
from final_finalizer.detection.contaminant import detect_contaminants
from final_finalizer.detection.debris import detect_chromosome_debris
from final_finalizer.detection.organelle import (
    detect_organelles,
    extract_organelle_from_ref,
    prepare_organelle_references,
)
from final_finalizer.detection.rdna import detect_rdna_contigs, prepare_rdna_reference
from final_finalizer.detection.rdna_consensus import build_rdna_consensus

__all__ = [
    # blast.py
    "run_makeblastdb",
    "run_blastn_megablast",
    "parse_blast_coverage",
    # organelle.py
    "extract_organelle_from_ref",
    "prepare_organelle_references",
    "detect_organelles",
    # rdna.py
    "prepare_rdna_reference",
    "detect_rdna_contigs",
    # rdna_consensus.py
    "build_rdna_consensus",
    # debris.py
    "detect_chromosome_debris",
    # contaminant.py
    "detect_contaminants",
]
