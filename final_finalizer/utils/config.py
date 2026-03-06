#!/usr/bin/env python3
"""
TOML configuration file support for final_finalizer.

Provides functions for loading configuration from TOML files,
merging with command-line arguments, and generating config templates.

Configuration files allow users to:
- Store complex parameter sets for reproducible analyses
- Share analysis parameters across projects
- Override specific parameters via CLI while keeping base config

The CONFIG_SCHEMA defines the mapping between TOML sections/keys
and argparse argument names, ensuring proper validation and type conversion.

Usage:
    # Generate template
    ./final_finalizer.py --dump-config > config.toml

    # Use config file (CLI args override config values)
    ./final_finalizer.py --config config.toml -t 64
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Optional

import tomllib


# Schema mapping TOML sections/keys to argparse argument names
# Format: {"section": {"toml_key": "argparse_dest"}}
CONFIG_SCHEMA: Dict[str, Dict[str, str]] = {
    "required": {
        "ref": "ref",
        "query": "query",
        "output_dir": "output_dir",
        "ref_gff3": "ref_gff3",
    },
    "common": {
        "threads": "threads",
        "skip_plot": "skip_plot",
        "chr_like_minlen": "chr_like_minlen",
        "add_subgenome_suffix": "add_subgenome_suffix",
        "ref_id_pattern": "ref_id_pattern",
    },
    "read_depth": {
        "reads": "reads",
        "reads_type": "reads_type",
        "skip_depth": "skip_depth",
        "depth_window_size": "depth_window_size",
        "depth_target_coverage": "depth_target_coverage",
        "keep_depth_bam": "keep_depth_bam",
    },
    "synteny": {
        "synteny_mode": "synteny_mode",
    },
    "pipeline_toggles": {
        "skip_organelles": "skip_organelles",
        "skip_rdna": "skip_rdna",
        "skip_contaminants": "skip_contaminants",
    },
    "external_tools": {
        "miniprot": "miniprot",
        "miniprot_args": "miniprot_args",
    },
    "reference_inputs": {
        "chrC_ref": "chrC_ref",
        "chrM_ref": "chrM_ref",
        "rdna_ref": "rdna_ref",
        "centrifuger_idx": "centrifuger_idx",
    },
    "thresholds_chromosome": {
        "assign_min_frac": "assign_min_frac",
        "assign_min_ratio": "assign_min_ratio",
        "chimera_primary_frac": "chimera_primary_frac",
        "chimera_secondary_frac": "chimera_secondary_frac",
        "low_ref_span_threshold": "low_ref_span_threshold",
    },
    "thresholds_synteny": {
        "assign_minlen": "assign_minlen",
        "assign_minlen_protein": "assign_minlen_protein",
        "assign_minmapq": "assign_minmapq",
        "assign_min_ident": "assign_min_ident",
        "assign_tp": "assign_tp",
    },
    "thresholds_chain": {
        "chain_q_gap": "chain_q_gap",
        "chain_r_gap": "chain_r_gap",
        "chain_diag_slop": "chain_diag_slop",
        "assign_chain_min_bp": "assign_chain_min_bp",
        "assign_chain_score": "assign_chain_score",
        "assign_chain_topk": "assign_chain_topk",
        "assign_ref_score": "assign_ref_score",
    },
    "thresholds_protein": {
        "miniprot_min_genes": "miniprot_min_genes",
        "miniprot_min_segments": "miniprot_min_segments",
        "min_span_frac": "min_span_frac",
        "min_span_bp": "min_span_bp",
    },
    "thresholds_organelle": {
        "organelle_min_cov": "organelle_min_cov",
        "chrC_len_tolerance": "chrC_len_tolerance",
        "chrM_len_tolerance": "chrM_len_tolerance",
    },
    "thresholds_rdna": {
        "rdna_min_cov": "rdna_min_cov",
        "skip_rdna_consensus": "skip_rdna_consensus",
        "rdna_ref_features": "rdna_ref_features",
    },
    "thresholds_chromosome_debris": {
        "chr_debris_min_cov": "chr_debris_min_cov",
        "chr_debris_min_identity": "chr_debris_min_identity",
    },
    "thresholds_contaminant": {
        "contaminant_min_score": "contaminant_min_score",
        "contaminant_min_coverage": "contaminant_min_coverage",
    },
    "thresholds_debris": {
        "debris_min_cov": "debris_min_cov",
        "debris_min_protein_hits": "debris_min_protein_hits",
    },
    "minimap2_tuning": {
        "preset": "preset",
        "kmer": "kmer",
        "window": "window",
        "aln_minlen": "aln_minlen",
    },
    "full_length_detection": {
        "full_length_ref_coverage": "full_length_ref_coverage",
        "skip_telomeres": "skip_telomeres",
        "telomere_motif": "telomere_motif",
        "telomere_window": "telomere_window",
        "telomere_min_repeats": "telomere_min_repeats",
        "subgenome_k": "subgenome_k",
    },
    "scaffolding": {
        "scaffold": "scaffold",
        "scaffold_gap_size": "scaffold_gap_size",
    },
    "compleasm": {
        "compleasm_lineage": "compleasm_lineage",
        "compleasm_library": "compleasm_library",
        "compleasm_path": "compleasm_path",
        "skip_compleasm": "skip_compleasm",
    },
    "multi_assembly": {
        "fofn": "fofn",
        "assembly_dir": "assembly_dir",
    },
    "distributed": {
        "cluster": "cluster",
        "max_threads_dist": "max_threads_dist",
        "max_mem_dist": "max_mem_dist",
        "max_time_dist": "max_time_dist",
        "partition": "partition",
        "qos": "qos",
    },
}


def load_config(path: Path) -> Dict[str, Any]:
    """Load configuration from a TOML file.

    Args:
        path: Path to the TOML configuration file

    Returns:
        Flat dictionary of configuration values with argparse dest names as keys

    Raises:
        ValueError: If the configuration file is invalid
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Configuration file not found: {path}")

    with open(path, "rb") as f:
        raw_config = tomllib.load(f)

    # Flatten TOML structure to argparse dest names
    flat_config: Dict[str, Any] = {}

    for section, keys in CONFIG_SCHEMA.items():
        if section not in raw_config:
            continue
        section_data = raw_config[section]
        for toml_key, argparse_dest in keys.items():
            if toml_key in section_data:
                value = section_data[toml_key]
                # Convert paths to Path objects for path-like arguments
                if argparse_dest in ("ref", "query", "ref_gff3", "reads", "chrC_ref", "chrM_ref", "rdna_ref", "fofn", "assembly_dir"):
                    if value is not None:
                        value = Path(value)
                flat_config[argparse_dest] = value

    return flat_config


def merge_config_with_args(config: Dict[str, Any], args: argparse.Namespace) -> None:
    """Merge configuration from TOML file with command-line arguments.

    Command-line arguments take precedence over config file values.
    Only sets config values if the argparse argument has its default value.

    Args:
        config: Configuration dictionary from load_config()
        args: Parsed command-line arguments (modified in place)
    """
    # Get default values by creating a new parser and getting defaults
    # This is a simplified approach - we check if the current value differs from None
    # or common defaults
    for key, value in config.items():
        if not hasattr(args, key):
            continue

        current_value = getattr(args, key)

        # Skip if CLI explicitly set a value (non-None for optional args)
        # This is a heuristic - works for most cases
        if current_value is not None and key not in (
            "threads", "skip_plot",
            "max_threads_dist", "max_mem_dist", "max_time_dist",
            "partition", "qos",
        ):
            continue

        # For boolean flags that default to False, only override if not set via CLI
        if isinstance(value, bool) and current_value is False:
            setattr(args, key, value)
        # For numeric defaults, be more careful
        elif key == "threads" and current_value == 8:  # Default threads
            setattr(args, key, value)
        elif current_value is None:
            setattr(args, key, value)


def dump_config_template(args: argparse.Namespace) -> str:
    """Generate a TOML configuration template based on current argument values.

    Args:
        args: Parsed command-line arguments

    Returns:
        TOML-formatted configuration template string
    """
    lines = [
        "# final_finalizer configuration file",
        "# Generated template - customize values as needed",
        "",
    ]

    def format_value(val: Any) -> str:
        """Format a Python value for TOML output."""
        if val is None:
            return '""  # Not set'
        if isinstance(val, bool):
            return "true" if val else "false"
        if isinstance(val, str):
            return f'"{val}"'
        if isinstance(val, Path):
            return f'"{val}"'
        if isinstance(val, (int, float)):
            return str(val)
        if isinstance(val, list):
            formatted = ", ".join(f'"{v}"' for v in val)
            return f"[{formatted}]"
        return f'"{val}"'

    for section, keys in CONFIG_SCHEMA.items():
        section_title = section.replace("_", " ").title()
        lines.append(f"# {section_title}")
        lines.append(f"[{section}]")

        for toml_key, argparse_dest in keys.items():
            if hasattr(args, argparse_dest):
                value = getattr(args, argparse_dest)
                formatted = format_value(value)
                # Use underscores in TOML keys to match schema
                lines.append(f"{toml_key} = {formatted}")

        lines.append("")

    return "\n".join(lines)
