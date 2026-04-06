#!/usr/bin/env python3
"""
Contaminant detection for dnadis.

Contains functions for identifying contaminated contigs using centrifuger
taxonomic classification.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from dnadis.models import ContaminantHit, ContaminantHitExtended
from dnadis.utils.io_utils import have_exe, open_maybe_gzip
from dnadis.utils.logging_config import get_logger

logger = get_logger("contaminant")


# NCBI taxonomy rank names used by taxonkit
# Note: "kingdom" here will be populated with domain (Bacteria, Archaea) or
# Cavalier-Smith kingdoms for eukaryotes (Animalia, Plantae, Fungi, Chromista, Protozoa)
TAXONOMY_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Mapping of NCBI eukaryotic superkingdoms to Cavalier-Smith kingdoms
CAVALIER_SMITH_KINGDOMS = {
    # Animalia
    "Metazoa": "Animalia",
    # Fungi
    "Fungi": "Fungi",
    # Plantae (green plants + red/glaucophyte algae)
    "Viridiplantae": "Plantae",
    "Rhodophyta": "Plantae",
    "Glaucophyta": "Plantae",
    # Chromista (stramenopiles, alveolates, rhizaria, haptophytes, cryptophytes)
    "Stramenopiles": "Chromista",
    "Alveolata": "Chromista",
    "Rhizaria": "Chromista",
    "Haptista": "Chromista",
    "Cryptophyceae": "Chromista",
    "Sar": "Chromista",
    # Protozoa (other single-celled eukaryotes)
    "Amoebozoa": "Protozoa",
    "Excavata": "Protozoa",
    "Opisthokonta": "Protozoa",  # If not Metazoa or Fungi
}


def _fasta_to_fastq_stream(fasta_path: Path, fastq_path: Path) -> None:
    """Convert FASTA to FASTQ with dummy quality scores.

    Centrifuger requires FASTQ input. For contigs, we use a constant high
    quality score (I = Phred 40) since these are assembled sequences.
    """
    with open_maybe_gzip(fasta_path, "rt") as fh_in, fastq_path.open("w") as fh_out:
        name: Optional[str] = None
        seq_parts: List[str] = []

        for line in fh_in:
            if line.startswith(">"):
                # Write previous record
                if name is not None and seq_parts:
                    seq = "".join(seq_parts)
                    fh_out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

        # Write last record
        if name is not None and seq_parts:
            seq = "".join(seq_parts)
            fh_out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def _validate_centrifuger_index(idx_prefix: str) -> bool:
    """Check if centrifuger index files exist."""
    idx_path = Path(idx_prefix)
    # Centrifuger index consists of .1.cfr, .2.cfr, .3.cfr files
    required_exts = [".1.cfr", ".2.cfr", ".3.cfr"]
    for ext in required_exts:
        if not Path(str(idx_path) + ext).exists():
            return False
    return True


def _get_centrifuger_name_table(idx_prefix: str, work_dir: Path) -> Dict[int, str]:
    """Get taxid -> scientific name mapping from centrifuger index.

    Runs centrifuger-inspect --name-table and caches the result.
    Returns dict mapping taxid (int) -> scientific name (str).
    """
    cache_path = work_dir / "centrifuger_names.tsv"

    # Use cached result if available
    if cache_path.exists():
        name_table: Dict[int, str] = {}
        with cache_path.open("r") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    try:
                        name_table[int(parts[0])] = parts[1]
                    except ValueError:
                        continue
        return name_table

    # Run centrifuger-inspect
    if not have_exe("centrifuger-inspect"):
        logger.warning("centrifuger-inspect not found, cannot look up scientific names")
        return {}

    cmd = ["centrifuger-inspect", "-x", idx_prefix, "--name-table"]
    try:
        ret = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if ret.returncode != 0:
            logger.warning(f"centrifuger-inspect failed: {ret.stderr}")
            return {}
    except Exception as e:
        logger.warning(f"centrifuger-inspect error: {e}")
        return {}

    # Parse output: "taxid | name | scientific name |"
    name_table = {}
    for line in ret.stdout.splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) >= 2:
            try:
                taxid = int(parts[0])
                sci_name = parts[1]
                name_table[taxid] = sci_name
            except ValueError:
                continue

    # Cache result
    work_dir.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w") as fh:
        for taxid, name in sorted(name_table.items()):
            fh.write(f"{taxid}\t{name}\n")

    logger.info(f"Loaded {len(name_table)} taxid -> name mappings from centrifuger index")
    return name_table


def _get_taxonomic_lineage(
    taxids: List[int],
    work_dir: Path,
) -> Dict[int, Dict[str, Optional[str]]]:
    """Query TaxonKit for taxonomic lineage of given taxids.

    Uses 'taxonkit lineage' and 'taxonkit reformat' to get full lineage
    with standardized ranks (kingdom, phylum, class, order, family, genus, species).

    Args:
        taxids: List of NCBI taxonomy IDs to look up
        work_dir: Working directory for intermediate files

    Returns:
        Dict mapping taxid -> {kingdom, phylum, class, order, family, genus, species}.
        Falls back to empty dict if taxonkit not available.
    """
    if not taxids:
        return {}

    if not have_exe("taxonkit"):
        logger.warning("taxonkit not found; taxonomic lineage unavailable for alluvial plot")
        return {}

    work_dir.mkdir(parents=True, exist_ok=True)
    cache_path = work_dir / "taxonkit_lineage.tsv"

    # Check for cached results
    if cache_path.exists():
        lineages: Dict[int, Dict[str, Optional[str]]] = {}
        try:
            with cache_path.open("r") as fh:
                header = fh.readline().strip().split("\t")
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    try:
                        tid = int(parts[0])
                        lineage_dict: Dict[str, Optional[str]] = {}
                        for i, rank in enumerate(TAXONOMY_RANKS):
                            if i + 1 < len(parts) and parts[i + 1]:
                                lineage_dict[rank] = parts[i + 1]
                            else:
                                lineage_dict[rank] = None
                        lineages[tid] = lineage_dict
                    except ValueError:
                        continue
            # Check if cache has all requested taxids
            if all(tid in lineages for tid in taxids):
                logger.info(f"Using cached taxonomic lineage for {len(lineages)} taxids")
                return lineages
        except Exception as e:
            logger.warning(f"Failed to read lineage cache: {e}")

    # Write taxids to temp file
    taxid_file = work_dir / "taxids.txt"
    with taxid_file.open("w") as f:
        for tid in set(taxids):
            f.write(f"{tid}\n")

    # Run taxonkit lineage -n -R to get lineage with rank names
    # Then reformat with specific ranks
    try:
        # Step 1: Get lineage with taxid
        lineage_cmd = ["taxonkit", "lineage", "-n", str(taxid_file)]
        lineage_result = subprocess.run(lineage_cmd, capture_output=True, text=True, check=False)

        if lineage_result.returncode != 0:
            logger.warning(f"taxonkit lineage failed: {lineage_result.stderr}")
            return {}

        # Write intermediate lineage output for reformat
        lineage_tmp = work_dir / "lineage_raw.tsv"
        with lineage_tmp.open("w") as f:
            f.write(lineage_result.stdout)

        # Step 2: Reformat to standardized ranks
        # Use {d} for domain (Bacteria, Archaea, Eukaryota, Viruses)
        # Use {K} for eukaryotic superkingdom (Metazoa, Fungi, Viridiplantae, etc.)
        # Format: domain;superkingdom;phylum;class;order;family;genus;species
        reformat_cmd = [
            "taxonkit", "reformat",
            "-I", "1",  # taxid column
            "-f", "{d};{K};{p};{c};{o};{f};{g};{s}",
            str(lineage_tmp)
        ]
        reformat_result = subprocess.run(reformat_cmd, capture_output=True, text=True, check=False)

        if reformat_result.returncode != 0:
            logger.warning(f"taxonkit reformat failed: {reformat_result.stderr}")
            return {}

    except Exception as e:
        logger.warning(f"taxonkit error: {e}")
        return {}

    # Parse reformat output
    # Format: taxid \t name \t lineage \t formatted_lineage
    # formatted_lineage: domain;superkingdom;phylum;class;order;family;genus;species
    lineages = {}
    for line in reformat_result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            tid = int(parts[0])
            formatted = parts[3] if len(parts) > 3 else ""
            ranks = formatted.split(";") if formatted else []

            # ranks[0] = domain (Bacteria, Archaea, Eukaryota, Viruses)
            # ranks[1] = superkingdom (for eukaryotes: Metazoa, Fungi, Viridiplantae, etc.)
            # ranks[2:] = phylum, class, order, family, genus, species
            domain = ranks[0] if len(ranks) > 0 and ranks[0] else None
            superkingdom = ranks[1] if len(ranks) > 1 and ranks[1] else None

            # Determine kingdom value
            if domain in ("Bacteria", "Archaea", "Viruses"):
                kingdom = domain
            elif domain == "Eukaryota" and superkingdom:
                # Map to Cavalier-Smith kingdoms
                kingdom = CAVALIER_SMITH_KINGDOMS.get(superkingdom, "Protozoa")
            else:
                kingdom = domain or superkingdom

            lineage_dict: Dict[str, Optional[str]] = {"kingdom": kingdom}
            # Map remaining ranks (phylum onwards starts at index 2)
            rank_names = ["phylum", "class", "order", "family", "genus", "species"]
            for i, rank in enumerate(rank_names):
                idx = i + 2  # offset by domain and superkingdom
                if idx < len(ranks) and ranks[idx] and ranks[idx] != "":
                    lineage_dict[rank] = ranks[idx]
                else:
                    lineage_dict[rank] = None
            lineages[tid] = lineage_dict
        except ValueError:
            continue

    # Cache results
    try:
        with cache_path.open("w") as fh:
            fh.write("taxid\t" + "\t".join(TAXONOMY_RANKS) + "\n")
            for tid, lineage_dict in sorted(lineages.items()):
                row = [str(tid)]
                for rank in TAXONOMY_RANKS:
                    row.append(lineage_dict.get(rank) or "")
                fh.write("\t".join(row) + "\n")
        logger.info(f"Cached taxonomic lineage for {len(lineages)} taxids")
    except Exception as e:
        logger.warning(f"Failed to cache lineage: {e}")

    return lineages


def _get_genus_lineages(
    genera: Set[str],
    work_dir: Path,
) -> Dict[str, Dict[str, Optional[str]]]:
    """Look up taxonomic lineage for genus names using taxonkit name2taxid.

    This is a fallback for GTDB-specific taxids that don't exist in NCBI taxonomy.
    By looking up the genus name, we can still get kingdom/phylum/family information.

    Args:
        genera: Set of genus names to look up
        work_dir: Working directory for intermediate files

    Returns:
        Dict mapping genus name -> {kingdom, phylum, class, order, family, genus, species}
    """
    if not genera:
        return {}

    cache_path = work_dir / "genus_lineage_cache.tsv"

    # Check cache
    cached: Dict[str, Dict[str, Optional[str]]] = {}
    if cache_path.exists():
        try:
            with cache_path.open("r") as fh:
                header = fh.readline()
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        genus_name = parts[0]
                        lineage_dict: Dict[str, Optional[str]] = {}
                        for i, rank in enumerate(TAXONOMY_RANKS):
                            if i + 1 < len(parts) and parts[i + 1]:
                                lineage_dict[rank] = parts[i + 1]
                            else:
                                lineage_dict[rank] = None
                        cached[genus_name] = lineage_dict
        except Exception:
            pass

    # Find genera not in cache
    to_lookup = genera - set(cached.keys())
    if not to_lookup:
        return cached

    # Write genus names to temp file
    genus_file = work_dir / "genera_to_lookup.txt"
    with genus_file.open("w") as f:
        for g in to_lookup:
            f.write(f"{g}\n")

    try:
        # Step 1: Convert genus names to taxids
        name2taxid_cmd = ["taxonkit", "name2taxid", str(genus_file)]
        name2taxid_result = subprocess.run(name2taxid_cmd, capture_output=True, text=True, check=False)

        if name2taxid_result.returncode != 0:
            logger.warning(f"taxonkit name2taxid failed: {name2taxid_result.stderr}")
            return cached

        # Parse name -> taxid mapping
        genus_taxids: Dict[str, int] = {}
        for line in name2taxid_result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1]:
                try:
                    genus_taxids[parts[0]] = int(parts[1])
                except ValueError:
                    pass

        if not genus_taxids:
            return cached

        # Step 2: Get lineage for these taxids
        taxid_file = work_dir / "genus_taxids.txt"
        with taxid_file.open("w") as f:
            for tid in genus_taxids.values():
                f.write(f"{tid}\n")

        lineage_cmd = ["taxonkit", "lineage", "-n", str(taxid_file)]
        lineage_result = subprocess.run(lineage_cmd, capture_output=True, text=True, check=False)

        if lineage_result.returncode != 0:
            return cached

        lineage_tmp = work_dir / "genus_lineage_raw.tsv"
        with lineage_tmp.open("w") as f:
            f.write(lineage_result.stdout)

        # Step 3: Reformat to standardized ranks (same format as main function)
        reformat_cmd = [
            "taxonkit", "reformat",
            "-I", "1",
            "-f", "{d};{K};{p};{c};{o};{f};{g};{s}",
            str(lineage_tmp)
        ]
        reformat_result = subprocess.run(reformat_cmd, capture_output=True, text=True, check=False)

        if reformat_result.returncode != 0:
            return cached

        # Parse reformatted output and map back to genus names
        taxid_to_lineage: Dict[int, Dict[str, Optional[str]]] = {}
        for line in reformat_result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                tid = int(parts[0])
                formatted = parts[3] if len(parts) > 3 else ""
                ranks = formatted.split(";") if formatted else []

                # Same parsing logic as main function
                domain = ranks[0] if len(ranks) > 0 and ranks[0] else None
                superkingdom = ranks[1] if len(ranks) > 1 and ranks[1] else None

                if domain in ("Bacteria", "Archaea", "Viruses"):
                    kingdom = domain
                elif domain == "Eukaryota" and superkingdom:
                    kingdom = CAVALIER_SMITH_KINGDOMS.get(superkingdom, "Protozoa")
                else:
                    kingdom = domain or superkingdom

                lineage_dict: Dict[str, Optional[str]] = {"kingdom": kingdom}
                # Map remaining ranks (phylum onwards starts at index 2)
                rank_names = ["phylum", "class", "order", "family", "genus", "species"]
                for i, rank in enumerate(rank_names):
                    idx = i + 2  # offset by domain and superkingdom
                    if idx < len(ranks) and ranks[idx]:
                        lineage_dict[rank] = ranks[idx]
                    else:
                        lineage_dict[rank] = None
                taxid_to_lineage[tid] = lineage_dict
            except ValueError:
                continue

        # Map back to genus names
        for genus_name, tid in genus_taxids.items():
            if tid in taxid_to_lineage:
                cached[genus_name] = taxid_to_lineage[tid]

        # Update cache
        with cache_path.open("w") as fh:
            fh.write("genus\t" + "\t".join(TAXONOMY_RANKS) + "\n")
            for genus_name, lineage_dict in sorted(cached.items()):
                row = [genus_name]
                for rank in TAXONOMY_RANKS:
                    row.append(lineage_dict.get(rank) or "")
                fh.write("\t".join(row) + "\n")

        logger.info(f"Resolved {len(genus_taxids)} genus names to NCBI taxonomy")

    except Exception as e:
        logger.warning(f"Genus lineage lookup failed: {e}")

    return cached


def detect_contaminants(
    query_fasta: Path,
    query_lengths: Dict[str, int],
    centrifuger_idx: str,
    work_dir: Path,
    threads: int,
    min_score: int,
    exclude_contigs: Set[str],
    fetch_lineage: bool = True,
) -> Dict[str, ContaminantHitExtended]:
    """Identify contaminant contigs using centrifuger taxonomic classification.

    Centrifuger is much faster than BLAST for taxonomic classification.
    Any contig that gets a significant classification hit is considered
    a potential contaminant.

    Args:
        query_fasta: Query FASTA file
        query_lengths: Dict of contig name -> length
        centrifuger_idx: Path prefix to centrifuger index
        work_dir: Working directory for intermediate files
        threads: Number of threads
        min_score: Minimum centrifuger score to consider a hit significant
        exclude_contigs: Set of contigs to exclude from screening
        fetch_lineage: Whether to fetch taxonomic lineage via taxonkit

    Returns:
        Dict mapping contig_name -> ContaminantHitExtended with taxid, name, coverage, score, and lineage
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Check for centrifuger executable
    if not have_exe("centrifuger"):
        logger.warning("centrifuger not found in PATH, skipping contaminant detection")
        return {}

    # Validate index
    if not _validate_centrifuger_index(centrifuger_idx):
        logger.warning(f"centrifuger index not found: {centrifuger_idx}")
        return {}

    logger.info(f"Using centrifuger index: {centrifuger_idx}")

    # Convert FASTA to FASTQ (centrifuger requires FASTQ)
    fastq_path = work_dir / "contigs.fq"
    if not fastq_path.exists():
        logger.info("Converting FASTA to FASTQ for centrifuger")
        _fasta_to_fastq_stream(query_fasta, fastq_path)

    # Run centrifuger
    output_path = work_dir / "centrifuger_results.tsv"
    err_path = work_dir / "centrifuger.err"

    if not output_path.exists():
        cmd = [
            "centrifuger",
            "-x", centrifuger_idx,
            "-u", str(fastq_path),
            "-t", str(threads),
        ]

        logger.info(f"Running centrifuger -> {output_path}")

        err_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w") as out_fh, err_path.open("w") as err_fh:
            ret = subprocess.run(cmd, stdout=out_fh, stderr=err_fh, check=False)

        if ret.returncode != 0:
            logger.warning(f"centrifuger failed with return code {ret.returncode}")
            return {}

    # Parse results
    contaminants: Dict[str, ContaminantHitExtended] = {}

    if not output_path.exists() or output_path.stat().st_size == 0:
        return contaminants

    # Get taxid -> scientific name mapping
    name_table = _get_centrifuger_name_table(centrifuger_idx, work_dir)

    # Temporary storage for basic hit info before lineage lookup
    hits_basic: List[Tuple[str, int, str, float, int]] = []  # (contig, taxid, sci_name, coverage, score)

    # Centrifuger output columns:
    # 1: Read ID, 2: Sequence ID, 3: Taxid, 4: Score, 5: Second-best score,
    # 6: Matching bp, 7: Read length, 8: Number of classifications
    with output_path.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue

            contig_name = fields[0]
            if contig_name in exclude_contigs:
                continue

            try:
                taxid = int(fields[2])
                score = int(fields[3])
                matching_bp = int(fields[5])
                read_len = int(fields[6])
            except ValueError:
                continue

            # Skip low-score hits
            if score < min_score:
                continue

            # Calculate coverage
            coverage = matching_bp / read_len if read_len > 0 else 0.0

            # Look up scientific name
            sci_name = name_table.get(taxid, "")

            hits_basic.append((contig_name, taxid, sci_name, coverage, score))
            contig_len = query_lengths.get(contig_name, read_len)
            logger.info(f"Contaminant: {contig_name} ({contig_len:,} bp, taxid={taxid}, {sci_name}, score={score}, cov={coverage:.2f})")

    # Fetch taxonomic lineage for all taxids
    lineages: Dict[int, Dict[str, Optional[str]]] = {}
    if fetch_lineage and hits_basic:
        unique_taxids = list(set(taxid for _, taxid, _, _, _ in hits_basic))
        lineages = _get_taxonomic_lineage(unique_taxids, work_dir)

        # For GTDB-specific taxids that returned empty lineage, try looking up genus
        # Extract genus names from sci_name for entries with missing lineage
        genera_to_lookup: Set[str] = set()
        for contig_name, taxid, sci_name, coverage, score in hits_basic:
            lineage = lineages.get(taxid, {})
            # Check if lineage is empty (no phylum = probably GTDB-specific taxid)
            if not lineage.get("phylum") and sci_name:
                # Parse genus from sci_name
                genus = sci_name.replace("_", " ").split()[0] if sci_name else None
                if genus and genus not in genera_to_lookup:
                    genera_to_lookup.add(genus)

        # Look up genus names in NCBI taxonomy using taxonkit name2taxid
        if genera_to_lookup and have_exe("taxonkit"):
            genus_lineages = _get_genus_lineages(genera_to_lookup, work_dir)
            # Store for later use when building extended hits
            lineages["__genus_fallback__"] = genus_lineages  # type: ignore

    # Get genus fallback lineages if available
    genus_fallback: Dict[str, Dict[str, Optional[str]]] = {}
    if "__genus_fallback__" in lineages:
        genus_fallback = lineages.pop("__genus_fallback__")  # type: ignore

    # Build extended contaminant hits
    for contig_name, taxid, sci_name, coverage, score in hits_basic:
        lineage = lineages.get(taxid, {})

        # Fallback: parse genus from sci_name if taxonkit didn't provide it
        # sci_name format is typically "Genus_species" or "Genus species"
        genus = lineage.get("genus")
        species = lineage.get("species")
        parsed_genus = None
        if sci_name:
            parts = sci_name.replace("_", " ").split()
            if parts:
                parsed_genus = parts[0]
                if not genus:
                    genus = parsed_genus
                if not species and len(parts) > 1:
                    species = " ".join(parts[1:])

        # Use genus-level lineage as fallback for GTDB-specific taxids
        # (when species-level lookup returned empty)
        kingdom = lineage.get("kingdom")
        phylum = lineage.get("phylum")
        tax_class = lineage.get("class")
        order = lineage.get("order")
        family = lineage.get("family")

        if not phylum and parsed_genus and parsed_genus in genus_fallback:
            genus_lineage = genus_fallback[parsed_genus]
            kingdom = kingdom or genus_lineage.get("kingdom")
            phylum = phylum or genus_lineage.get("phylum")
            tax_class = tax_class or genus_lineage.get("class")
            order = order or genus_lineage.get("order")
            family = family or genus_lineage.get("family")

        contaminants[contig_name] = ContaminantHitExtended(
            taxid=taxid,
            sci_name=sci_name,
            coverage=coverage,
            score=score,
            kingdom=kingdom,
            phylum=phylum,
            tax_class=tax_class,
            order=order,
            family=family,
            genus=genus,
            species=species,
        )

    logger.info(f"Contaminant contigs: {len(contaminants)}")
    return contaminants
