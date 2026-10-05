# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from viral_read_quantifier.star.index import _STAR_GENOME_CHR_BIN_NBITS_DEFAULT

logger = logging.getLogger(__name__)

_DEFAULT_OUTDIR = Path("output_dir")
_DEFAULT_STAR_THREADS = 4


def _add_star_index_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--genome-sa-index-nbases",
        type=int,
        default=None,
        help=(
            "STAR --genomeSAindexNbases (SA pre-indexing string length). "
            "Auto-computed as min(14, log2(GenomeLength)/2-1) when omitted. "
            "See https://github.com/alexdobin/STAR/blob/master/source/parametersDefault"
        ),
    )
    parser.add_argument(
        "--genome-chr-bin-nbits",
        type=int,
        default=None,
        help=(f"STAR --genomeChrBinNbits value (default: {_STAR_GENOME_CHR_BIN_NBITS_DEFAULT})"),
    )


def _add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=_DEFAULT_OUTDIR,
        help=f"Output directory (default: {_DEFAULT_OUTDIR})",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=_DEFAULT_STAR_THREADS,
        help=f"Number of threads (default: {_DEFAULT_STAR_THREADS})",
    )
    parser.add_argument(
        "--star-bin",
        type=Path,
        default=None,
        help="Path to STAR binary (default: search PATH via shutil.which)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if output already exists",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viral_read_quantifier",
        description="Virome analysis — RNA-seq read mapping and strand-specific quantification",
    )
    sub = parser.add_subparsers(dest="command")

    # --- index ---
    p_index = sub.add_parser("index", help="Build STAR genome index from FASTA")
    p_index.add_argument("--genome", type=Path, required=True, help="Genome FASTA file")
    _add_common_options(p_index)
    _add_star_index_options(p_index)

    # --- map ---
    p_map = sub.add_parser(
        "map",
        help="Run fastp QC, STAR mapping, strand separation, and quantification",
    )
    p_map.add_argument(
        "--sra-info",
        type=Path,
        required=True,
        help="SraRunInfo.csv (Run, LibraryLayout columns)",
    )
    p_map.add_argument(
        "--strandedness-info",
        type=Path,
        required=True,
        help="Strandedness CSV (Run, strandedness columns)",
    )
    p_map.add_argument("--index-dir", type=Path, default=None, help="STAR index directory")
    p_map.add_argument(
        "--fastq-dir",
        type=Path,
        default=None,
        help="Directory containing FASTQ files (default: <output-dir>/fastq)",
    )
    p_map.add_argument(
        "--alignment-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing pre-aligned BAM/SAM files "
            "(BAM must be coordinate-sorted; SAM is sorted internally. "
            "Mutually exclusive with --fastq-dir)"
        ),
    )
    _add_common_options(p_map)
    p_map.add_argument(
        "--limit-bam-sort-ram",
        type=int,
        default=None,
        help="STAR --limitBAMsortRAM (bytes). Set when BAM sorting runs out of memory",
    )

    # --- visualize ---
    p_viz = sub.add_parser(
        "visualize",
        help="Generate per-sample images (SVG + PDF) from strand depth/counts TSVs",
    )
    p_viz.add_argument(
        "--depth-tsv",
        type=Path,
        required=True,
        help="Merged strand_depth.tsv (output of 'map')",
    )
    p_viz.add_argument(
        "--counts-tsv",
        type=Path,
        required=True,
        help="Merged strand_counts.tsv (output of 'map')",
    )
    p_viz.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=_DEFAULT_OUTDIR,
        help=f"Output directory (default: {_DEFAULT_OUTDIR})",
    )

    # --- run ---
    p_run = sub.add_parser("run", help="Run full pipeline: index + map + visualize")
    p_run.add_argument("--genome", type=Path, default=None, help="Genome FASTA file")
    p_run.add_argument(
        "--sra-info",
        type=Path,
        required=True,
        help="SraRunInfo.csv (Run, LibraryLayout columns)",
    )
    p_run.add_argument(
        "--strandedness-info",
        type=Path,
        required=True,
        help="Strandedness CSV (Run, strandedness columns)",
    )
    p_run.add_argument(
        "--fastq-dir",
        type=Path,
        default=None,
        help="Directory containing FASTQ files (default: <output-dir>/fastq)",
    )
    p_run.add_argument(
        "--alignment-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing pre-aligned BAM/SAM files "
            "(BAM must be coordinate-sorted; SAM is sorted internally. "
            "Mutually exclusive with --fastq-dir)"
        ),
    )
    _add_common_options(p_run)
    _add_star_index_options(p_run)
    p_run.add_argument(
        "--limit-bam-sort-ram",
        type=int,
        default=None,
        help="STAR --limitBAMsortRAM (bytes). Set when BAM sorting runs out of memory",
    )

    return parser


def _resolve_fastq_dir(args: argparse.Namespace) -> Path:
    if args.fastq_dir is not None:
        return args.fastq_dir
    return args.output_dir / "fastq"


def _validate_input_mode(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    alignment_dir = getattr(args, "alignment_dir", None)
    fastq_dir = getattr(args, "fastq_dir", None)

    if alignment_dir is not None and fastq_dir is not None:
        parser.error("--alignment-dir and --fastq-dir are mutually exclusive")

    if alignment_dir is None:
        if args.command == "map" and getattr(args, "index_dir", None) is None:
            parser.error("--index-dir is required in FASTQ mode")
        if args.command == "run" and getattr(args, "genome", None) is None:
            parser.error("--genome is required in FASTQ mode")


def _dispatch_index(args: argparse.Namespace) -> None:
    from viral_read_quantifier.cmd_index import run_index

    run_index(
        genome_fasta=args.genome,
        output_dir=args.output_dir,
        star_bin=args.star_bin,
        threads=args.threads,
        force=args.force,
        genome_sa_index_nbases=args.genome_sa_index_nbases,
        genome_chr_bin_nbits=args.genome_chr_bin_nbits,
    )


def _dispatch_map(args: argparse.Namespace) -> None:
    if args.alignment_dir is not None:
        from viral_read_quantifier.cmd_map import run_map_from_alignments

        run_map_from_alignments(
            sra_info_csv=args.sra_info,
            strandedness_info_csv=args.strandedness_info,
            alignment_dir=args.alignment_dir,
            output_dir=args.output_dir,
            threads=args.threads,
            force=args.force,
        )
    else:
        from viral_read_quantifier.cmd_map import run_map

        run_map(
            sra_info_csv=args.sra_info,
            strandedness_info_csv=args.strandedness_info,
            index_dir=args.index_dir,
            fastq_dir=_resolve_fastq_dir(args),
            output_dir=args.output_dir,
            star_bin=args.star_bin,
            threads=args.threads,
            force=args.force,
            limit_bam_sort_ram_bytes=args.limit_bam_sort_ram,
        )


def _dispatch_visualize(args: argparse.Namespace) -> None:
    from viral_read_quantifier.cmd_visualize import run_visualize

    run_visualize(
        depth_tsv=args.depth_tsv,
        counts_tsv=args.counts_tsv,
        output_dir=args.output_dir,
    )


def _dispatch_run(args: argparse.Namespace) -> None:
    from viral_read_quantifier.cmd_map import (
        MERGED_STRAND_COUNTS_BASENAME,
        MERGED_STRAND_DEPTH_BASENAME,
    )
    from viral_read_quantifier.cmd_visualize import run_visualize

    if args.alignment_dir is None:
        from viral_read_quantifier.cmd_index import run_index

        logger.info("=== Step 1/3: index ===")
        args.index_dir = run_index(
            genome_fasta=args.genome,
            output_dir=args.output_dir,
            star_bin=args.star_bin,
            threads=args.threads,
            force=args.force,
            genome_sa_index_nbases=args.genome_sa_index_nbases,
            genome_chr_bin_nbits=args.genome_chr_bin_nbits,
        )
        logger.info("=== Step 2/3: map ===")
    else:
        logger.info("=== Step 1/2: map (alignment mode) ===")

    _dispatch_map(args)

    summary_dir = args.output_dir / "summary"
    step_label = "Step 3/3" if args.alignment_dir is None else "Step 2/2"
    logger.info("=== %s: visualize ===", step_label)
    run_visualize(
        depth_tsv=summary_dir / MERGED_STRAND_DEPTH_BASENAME,
        counts_tsv=summary_dir / MERGED_STRAND_COUNTS_BASENAME,
        output_dir=args.output_dir,
    )


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command in ("map", "run"):
        _validate_input_mode(args, parser)

    match args.command:
        case "index":
            _dispatch_index(args)
        case "map":
            _dispatch_map(args)
        case "visualize":
            _dispatch_visualize(args)
        case "run":
            _dispatch_run(args)
        case _:
            raise NotImplementedError(f"Subcommand '{args.command}' is not implemented")
