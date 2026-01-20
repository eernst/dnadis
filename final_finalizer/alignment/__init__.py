# final_finalizer/alignment/__init__.py
"""Alignment tools and chain parsing modules."""

from final_finalizer.alignment.chain_parsing import (
    _filter_overlapping_hits_by_identity,
    parse_miniprot_synteny_evidence_and_segments,
    parse_paf_chain_evidence_and_segments,
)
from final_finalizer.alignment.external_tools import (
    get_minimap2_exe,
    run_gffread_extract_proteins,
    run_minimap2,
    run_minimap2_synteny,
    run_miniprot,
)
from final_finalizer.alignment.stats import (
    parse_paf_primary,
    write_stats_tsv,
)

__all__ = [
    # external_tools.py
    "get_minimap2_exe",
    "run_minimap2",
    "run_minimap2_synteny",
    "run_gffread_extract_proteins",
    "run_miniprot",
    # chain_parsing.py
    "parse_paf_chain_evidence_and_segments",
    "parse_miniprot_synteny_evidence_and_segments",
    "_filter_overlapping_hits_by_identity",
    # stats.py
    "parse_paf_primary",
    "write_stats_tsv",
]
