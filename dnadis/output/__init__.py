# dnadis/output/__init__.py
"""Output modules for FASTA, TSV, and plotting."""

from dnadis.output.fasta_output import write_classified_fastas
from dnadis.output.tsv_output import (
    write_contig_summary_tsv,
    write_macro_blocks_tsv,
    write_chain_segments_tsv,
    write_chain_summary_tsv,
    compute_summary,
    build_segment_support_from_rows,
)
from dnadis.output.plotting import (
    have_rscript,
    run_assembly_report,
)

__all__ = [
    "write_classified_fastas",
    "write_contig_summary_tsv",
    "write_macro_blocks_tsv",
    "write_chain_segments_tsv",
    "write_chain_summary_tsv",
    "compute_summary",
    "build_segment_support_from_rows",
    "have_rscript",
    "run_assembly_report",
]
