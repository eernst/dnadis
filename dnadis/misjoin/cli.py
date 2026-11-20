"""Command line for ``dnadis misjoin {identify,break}``."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Set

from dnadis.misjoin import breaker
from dnadis.misjoin import identify as idf
from dnadis.misjoin import sisters as sis
from dnadis.utils.logging_config import get_logger, setup_logging

logger = get_logger("misjoin")


def _name_path_pairs(values: Optional[List[str]], what: str) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    for v in values or []:
        if "=" not in v:
            raise SystemExit(f"--{what} expects NAME=PATH, got {v!r}")
        name, path = v.split("=", 1)
        out[name] = Path(path)
    return out


def _fofn_paths(fofn: Optional[Path]) -> Dict[str, Path]:
    if fofn is None:
        return {}
    from dnadis.utils.multi_assembly import parse_fofn

    return {name: Path(path) for path, name, _ in parse_fofn(fofn)}


def _identify_parser(sub) -> None:
    p = sub.add_parser(
        "identify",
        help="find contigs that switch reference chromosome, and call misjoins",
        description=(
            "Find junctions where a contig switches from one reference chromosome to "
            "another (macro blocks), test each against the other assemblies in the "
            "results folder (sisters) and, optionally, against reads, and write one "
            "row per candidate break to <assembly>.misjoins.tsv."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("results_dir", type=Path, help="dnadis output folder (-o of the run)")
    p.add_argument("-a", "--assembly", action="append", default=None,
                   help="assembly to examine (repeatable; default: every assembly in the folder)")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="write <assembly>.misjoins.tsv here (default: each assembly's folder)")
    p.add_argument("--work-dir", type=Path, default=None,
                   help="alignment cache (default: RESULTS_DIR/misjoin_work)")
    p.add_argument("--fasta", action="append", metavar="NAME=PATH",
                   help="input FASTA of an assembly (default: read from its class FASTAs)")
    p.add_argument("--fofn", type=Path, default=None,
                   help="the run's FOFN, for every assembly's input FASTA")
    p.add_argument("--sister-results", type=Path, action="append", default=[],
                   help="another dnadis output folder whose assemblies serve as sisters (repeatable)")
    p.add_argument("--exclude-sister", action="append", default=[],
                   help="assembly not to use as a sister (repeatable), e.g. one built from "
                        "the same contigs")
    p.add_argument("--derived", action="append", default=[], metavar="NAME=OF",
                   help="NAME was built from OF's contigs (e.g. scaffolds of them), or OF from "
                        "NAME's; neither is a sister of the other (repeatable)")
    p.add_argument("--no-sisters", action="store_true", help="reference evidence only")
    p.add_argument("--bam", action="append", metavar="NAME=PATH",
                   help="reads aligned to an assembly's input contigs (indexed BAM/CRAM)")
    p.add_argument("--cram-reference", type=Path, default=None, help="reference FASTA for CRAM input")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--jobs", type=int, default=1,
                   help="sister alignments run at once (threads are shared among them)")

    d = idf.IdentifyParams()
    s = d.sisters
    r = d.reads
    g = p.add_argument_group("junctions")
    g.add_argument("--min-block-bp", type=int, default=d.min_block_bp,
                   help="ignore macro blocks with less aligned sequence")
    g.add_argument("--min-side-bp", type=int, default=d.min_side_bp,
                   help="aligned bp a reference chromosome needs on one side of a junction")
    g = p.add_argument_group("sister evidence")
    g.add_argument("--flank", type=int, default=s.flank,
                   help="bp of each side aligned to the sisters")
    g.add_argument("--min-identity", type=float, default=s.min_identity,
                   help="minimum gap-compressed identity of a sister hit")
    g.add_argument("--same-genotype-identity", type=float, default=s.same_genotype_identity,
                   help="a carrier this identical is the same genotype and does not count")
    g.add_argument("--min-extend", type=int, default=s.min_extend,
                   help="sister sequence beyond the junction needed to contradict it")
    g.add_argument("--min-contradict", type=int, default=d.min_contradict,
                   help="contradicting sister units needed to call a misjoin")
    g.add_argument("--repeat-min-copies", type=int, default=d.repeat_min_copies,
                   help="copies per sister that make the segment between the sides a repeat")
    g.add_argument("--min-chrom-frac", type=float, default=s.min_chrom_frac,
                   help="fraction of each junction chromosome a sister unit's chromosome "
                        "contigs must hold for the unit to count")
    g.add_argument("--min-genotype-carriers", type=int, default=d.min_genotype_carriers,
                   help="same-genotype carriers (with no same-genotype contradiction) that "
                        "make a join genotype_specific instead of a misjoin")
    g.add_argument("--min-side-density", type=float, default=d.min_side_density,
                   help="a junction whose two sides are both aligned over less than this "
                        "fraction of their extent is reviewed, not broken")
    g.add_argument("--max-insert-bp", type=int, default=d.max_insert_bp,
                   help="a shorter run between two runs of one chromosome is an insertion, "
                        "not a misjoin")
    g = p.add_argument_group("telomeres")
    g.add_argument("--telomere-motif", default=d.telomere_motif,
                   help="telomere repeat motif ('' to skip the telomere check)")
    g.add_argument("--telomere-min-bp", type=int, default=d.telomere_min_bp,
                   help="bp of telomere repeat runs that mark a chromosome end at a break")
    g.add_argument("--telomere-window", type=int, default=d.telomere_window,
                   help="distance from the break interval searched for telomere arrays")
    g = p.add_argument_group("read evidence")
    g.add_argument("--min-mapq", type=int, default=r.min_mapq)
    g.add_argument("--min-clip", type=int, default=r.min_clip)
    g.add_argument("--min-split", type=int, default=r.min_split,
                   help="split reads that confirm a break")
    g.add_argument("--min-spanning", type=int, default=r.min_spanning,
                   help="spanning reads that refute a break")
    p.add_argument("-v", "--verbose", action="store_true")


def _break_parser(sub) -> None:
    p = sub.add_parser(
        "break",
        help="split contigs at the breaks of an identify TSV",
        description=(
            "Split contigs of an input FASTA at the break positions of a misjoins TSV "
            "(rows with action 'break' unless --actions says otherwise). Pieces are "
            "named <contig>_b1, _b2, ...; the AGP maps them to the input contigs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fasta", type=Path, help="input assembly FASTA (the one dnadis read)")
    p.add_argument("tsv", type=Path, help="<assembly>.misjoins.tsv, or a filtered copy")
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="output FASTA (.gz: bgzip-compressed)")
    p.add_argument("--agp", type=Path, default=None,
                   help="output AGP (default: output name with .agp)")
    p.add_argument("--assembly", default=None, help="use only this assembly's rows")
    p.add_argument("--actions", default="break",
                   help="comma-separated actions to apply, or 'all' for every row")
    p.add_argument("-t", "--threads", type=int, default=1, help="compression threads")
    p.add_argument("-v", "--verbose", action="store_true")


def run_identify(args) -> int:
    found = sis.discover_assemblies(args.results_dir)
    if not found:
        raise SystemExit(f"No dnadis assembly outputs under {args.results_dir}")
    if args.assembly:
        unknown = sorted(set(args.assembly) - set(found))
        if unknown:
            raise SystemExit(f"Not in {args.results_dir}: {', '.join(unknown)}")
        targets = {a: found[a] for a in sorted(args.assembly)}
    else:
        targets = dict(sorted(found.items()))
    sisters = dict(found)
    for extra in args.sister_results:
        for name, a in sis.discover_assemblies(extra).items():
            if name in sisters:
                logger.warning(f"Sister {name} from {extra} shadows one of the same name; skipped")
                continue
            sisters[name] = a
    for name in args.exclude_sister:
        sisters.pop(name, None)

    fastas = _fofn_paths(args.fofn)
    fastas.update(_name_path_pairs(args.fasta, "fasta"))
    derived: Dict[str, Set[str]] = {}
    for item in args.derived:
        if "=" not in item:
            raise SystemExit(f"--derived expects NAME=OF, got {item}")
        a, b = item.split("=", 1)
        derived.setdefault(a, set()).add(b)
        derived.setdefault(b, set()).add(a)
    inputs = idf.IdentifyInputs(
        targets=targets, sisters=sisters, fastas=fastas,
        bams=_name_path_pairs(args.bam, "bam"), cram_reference=args.cram_reference,
        derived=derived,
    )
    params = idf.IdentifyParams(
        min_block_bp=args.min_block_bp, min_side_bp=args.min_side_bp,
        min_contradict=args.min_contradict, repeat_min_copies=args.repeat_min_copies,
        telomere_motif=args.telomere_motif, telomere_min_bp=args.telomere_min_bp,
        telomere_window=args.telomere_window,
        min_genotype_carriers=args.min_genotype_carriers, max_insert_bp=args.max_insert_bp,
        min_side_density=args.min_side_density,
    )
    params.sisters.min_chrom_frac = args.min_chrom_frac
    params.sisters.flank = args.flank
    params.sisters.min_identity = args.min_identity
    params.sisters.same_genotype_identity = args.same_genotype_identity
    params.sisters.min_extend = args.min_extend
    params.reads.min_mapq = args.min_mapq
    params.reads.min_clip = args.min_clip
    params.reads.min_split = args.min_split
    params.reads.min_spanning = args.min_spanning

    work = args.work_dir or (args.results_dir / "misjoin_work")
    calls = idf.identify(inputs, work, params, threads=args.threads, jobs=args.jobs,
                         use_sisters=not args.no_sisters)
    for name, jcs in calls.items():
        rows = [row for jc in jcs for row in idf.call_rows(jc)]
        out_dir = args.out_dir or targets[name].directory
        path = out_dir / f"{name}.misjoins.tsv"
        idf.write_tsv(rows, path)
        n_break = sum(1 for jc in jcs if jc.action == "break")
        n_review = sum(1 for jc in jcs if jc.action == "review")
        logger.done(f"{name}: {len(jcs)} junctions, {n_break} to break, {n_review} to review: {path}")
    return 0


def run_break(args) -> int:
    actions = None if args.actions == "all" else tuple(a for a in args.actions.split(",") if a)
    breaks = breaker.read_breaks(args.tsv, actions, args.assembly)
    agp = args.agp
    if agp is None:
        base = args.output.name
        for ext in (".gz", ".bgz"):
            if base.endswith(ext):
                base = base[: -len(ext)]
        for ext in (".fasta", ".fa", ".fna"):
            if base.endswith(ext):
                base = base[: -len(ext)]
        agp = args.output.with_name(base + ".agp")
    breaker.break_fasta(args.fasta, breaks, args.output, agp, threads=args.threads)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="dnadis misjoin",
        description="Find contigs that join two reference chromosomes, call misjoins "
                    "with sister-assembly and read evidence, and split them.",
    )
    sub = p.add_subparsers(dest="command", required=True)
    _identify_parser(sub)
    _break_parser(sub)
    args = p.parse_args(argv)
    setup_logging(verbose=args.verbose)
    if args.command == "identify":
        return run_identify(args)
    return run_break(args)
