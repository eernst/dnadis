# final_finalizer/output/__init__.py
"""Output modules for FASTA, TSV, and plotting."""

from final_finalizer.output.fasta_output import write_classified_fastas
from final_finalizer.output.tsv_output import (
    write_contig_summary_tsv,
    write_macro_blocks_tsv,
    write_chain_segments_tsv,
    write_chain_summary_tsv,
    compute_summary,
    build_segment_support_from_rows,
)
from final_finalizer.output.plotting import (
    have_rscript,
    run_plot,
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
    "run_plot",
]
