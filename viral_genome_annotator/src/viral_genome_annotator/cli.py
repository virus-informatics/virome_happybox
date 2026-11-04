# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

_DEFAULT_OUTDIR = "output_dir"
_DEFAULT_UPDATE_REFSEQ_OUTDIR = f"{_DEFAULT_OUTDIR}/update_refseq"
_DEFAULT_BLASTDB_PATH = f"{_DEFAULT_UPDATE_REFSEQ_OUTDIR}/blast_db/refseq_viral_nucl"


def _add_update_refseq_subcommands(
    parent: argparse._SubParsersAction,
) -> None:
    p_update = parent.add_parser(
        "update_refseq",
        help="Fetch and build reference data from NCBI FTP",
    )
    update_sub = p_update.add_subparsers(dest="update_refseq_command")

    outdir_help = f"Output directory (default: {_DEFAULT_UPDATE_REFSEQ_OUTDIR})"

    # download-refseq
    p_dl = update_sub.add_parser("download-refseq", help="Download files from NCBI FTP")
    p_dl.add_argument("-o", "--output-dir", default=_DEFAULT_UPDATE_REFSEQ_OUTDIR, help=outdir_help)
    p_dl.add_argument(
        "--update-taxdump",
        action="store_true",
        help="Also download/update NCBI taxdump (skipped by default)",
    )

    # download-assembly
    p_da = update_sub.add_parser(
        "download-assembly",
        help="Download GCF→NC assembly mapping from NCBI Datasets API",
    )
    p_da.add_argument("-o", "--output-dir", default=_DEFAULT_UPDATE_REFSEQ_OUTDIR, help=outdir_help)
    p_da.add_argument(
        "--batch-size",
        type=int,
        default=500,
        help="Accessions per API request (default: 500)",
    )

    # gen_meta
    p_gm = update_sub.add_parser("gen_meta", help="Generate GFF3 + taxonomy.tsv from GBFF")
    p_gm.add_argument("-o", "--output-dir", default=_DEFAULT_UPDATE_REFSEQ_OUTDIR, help=outdir_help)
    p_gm.add_argument(
        "--taxdump-dir",
        default=None,
        help="Extracted taxdump directory for taxonkit (default: <output-dir>/taxdump)",
    )

    # gen_db
    p_gd = update_sub.add_parser("gen_db", help="Build BLAST DB from FASTA")
    p_gd.add_argument("-o", "--output-dir", default=_DEFAULT_UPDATE_REFSEQ_OUTDIR, help=outdir_help)

    # run
    p_run_all = update_sub.add_parser("run", help="Run full pipeline: download + gen_meta + gen_db")
    p_run_all.add_argument(
        "-o",
        "--output-dir",
        default=_DEFAULT_UPDATE_REFSEQ_OUTDIR,
        help=outdir_help,
    )
    p_run_all.add_argument(
        "--taxdump-dir",
        default=None,
        help="Extracted taxdump directory for taxonkit (default: <output-dir>/taxdump)",
    )
    p_run_all.add_argument(
        "--update-taxdump",
        action="store_true",
        help="Also download/update NCBI taxdump (skipped by default)",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viral_genome_annotator",
        description="Virome analysis — genome annotation visualizer",
    )
    sub = parser.add_subparsers(dest="command")

    _outdir_help = f"Working directory (default: {_DEFAULT_OUTDIR})"

    # --- run ---
    p_run = sub.add_parser("run", help="Run blast + visualize in one step")
    p_run.add_argument("-i", "--input", required=True, help="Query FASTA file")
    p_run.add_argument("-o", "--outdir", default=_DEFAULT_OUTDIR, help=_outdir_help)
    p_run.add_argument(
        "--db",
        default=_DEFAULT_BLASTDB_PATH,
        help=f"BLAST DB path (default: {_DEFAULT_BLASTDB_PATH})",
    )
    p_run.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for tblastx (default: 4)",
    )
    p_run.add_argument(
        "--force",
        action="store_true",
        help="Re-run tblastx even if output TSV exists",
    )
    p_run.add_argument(
        "--refseq-dir",
        default=_DEFAULT_UPDATE_REFSEQ_OUTDIR,
        help=(
            "Reference data directory containing annotations.gff3,"
            " taxonomy.tsv, and assembly/assembly_map.json"
            f" (default: {_DEFAULT_UPDATE_REFSEQ_OUTDIR})"
        ),
    )
    p_run.add_argument(
        "--title",
        help="Figure title (default: auto-generated from BLAST TSV filename)",
    )
    p_run.add_argument(
        "--min-covered-bases",
        type=int,
        default=500,
        help="Minimum covered bases to include (default: 500)",
    )
    p_run.add_argument(
        "--max-track-count",
        type=int,
        default=15,
        help="Maximum tracks per figure (default: 15)",
    )

    # --- blast ---
    p_blast = sub.add_parser("blast", help="Run tblastx (TSV output)")
    p_blast.add_argument("-i", "--input", required=True, help="Query FASTA file")
    p_blast.add_argument("-o", "--outdir", default=_DEFAULT_OUTDIR, help=_outdir_help)
    p_blast.add_argument("--db", required=True, help="BLAST DB path")
    p_blast.add_argument(
        "--output-tsv",
        help="Output TSV path (default: auto-generated from FASTA header)",
    )
    p_blast.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    p_blast.add_argument("--force", action="store_true", help="Re-run even if output TSV exists")

    # --- visualize ---
    p_vis = sub.add_parser(
        "visualize",
        help="Prepare metadata, compute coverage, and render figures",
    )
    p_vis.add_argument(
        "-i",
        "--input",
        nargs="+",
        help="BLAST TSV files (multiple allowed; default: auto-detect from outdir/tblastx/)",
    )
    p_vis.add_argument("-o", "--outdir", default=_DEFAULT_OUTDIR, help=_outdir_help)
    p_vis.add_argument(
        "--refseq-dir",
        default=_DEFAULT_UPDATE_REFSEQ_OUTDIR,
        help=(
            "Reference data directory containing annotations.gff3,"
            " taxonomy.tsv, and assembly/assembly_map.json"
            f" (default: {_DEFAULT_UPDATE_REFSEQ_OUTDIR})"
        ),
    )
    p_vis.add_argument(
        "--title",
        help="Figure title (default: auto-generated from BLAST TSV filename)",
    )
    p_vis.add_argument(
        "--min-covered-bases",
        type=int,
        default=500,
        help="Minimum covered bases to include (default: 500)",
    )
    p_vis.add_argument(
        "--max-track-count",
        type=int,
        default=15,
        help="Maximum tracks per figure (default: 15)",
    )

    # --- update_refseq ---
    _add_update_refseq_subcommands(sub)

    return parser


def _dispatch_run(args: argparse.Namespace) -> None:
    from viral_genome_annotator.blast import run_tblastx
    from viral_genome_annotator.external import require_command
    from viral_genome_annotator.visualize import run_visualize

    require_command("tblastx")

    outdir = Path(args.outdir)
    blast_tsv = run_tblastx(
        query_fasta=Path(args.input),
        db_path=Path(args.db),
        workdir=outdir,
        num_threads=args.threads,
        force=args.force,
    )

    html_paths = run_visualize(
        blast_tsv,
        Path(args.refseq_dir),
        outdir,
        min_covered_bases=args.min_covered_bases,
        max_track_count=args.max_track_count,
        title=args.title,
    )
    logger.info("Generated %d HTML figure(s)", len(html_paths))


def _dispatch_blast(args: argparse.Namespace) -> None:
    from viral_genome_annotator.blast import run_tblastx
    from viral_genome_annotator.external import require_command

    require_command("tblastx")

    run_tblastx(
        query_fasta=Path(args.input),
        db_path=Path(args.db),
        workdir=Path(args.outdir),
        output_tsv=Path(args.output_tsv) if args.output_tsv else None,
        num_threads=args.threads,
        force=args.force,
    )


def _dispatch_visualize(args: argparse.Namespace) -> None:
    from viral_genome_annotator.visualize import discover_blast_tsv_files, run_visualize

    outdir = Path(args.outdir)
    refseq_dir = Path(args.refseq_dir)

    if args.input:
        blast_tsv_paths = [Path(p) for p in args.input]
    else:
        blast_tsv_paths = discover_blast_tsv_files(outdir)

    all_html_paths: list[Path] = []
    for blast_tsv_path in blast_tsv_paths:
        html_paths = run_visualize(
            blast_tsv_path,
            refseq_dir,
            outdir,
            min_covered_bases=args.min_covered_bases,
            max_track_count=args.max_track_count,
            title=args.title,
        )
        all_html_paths.extend(html_paths)
    logger.info("Generated %d HTML figure(s)", len(all_html_paths))


def _dispatch_update_refseq(args: argparse.Namespace) -> None:
    from viral_genome_annotator.update_refseq.cmd_download_assembly import run_download_assembly
    from viral_genome_annotator.update_refseq.cmd_download_refseq import run_download
    from viral_genome_annotator.update_refseq.cmd_gen_db import run_gen_db
    from viral_genome_annotator.update_refseq.cmd_gen_meta import run_gen_meta
    from viral_genome_annotator.update_refseq.cmd_run import run_all

    if args.update_refseq_command is None:
        raise SystemExit(
            "usage: viral_genome_annotator update_refseq"
            " {download-refseq,download-assembly,gen_meta,gen_db,run}",
        )

    output_dir = Path(args.output_dir)
    taxdump_dir = Path(args.taxdump_dir) if getattr(args, "taxdump_dir", None) else None

    match args.update_refseq_command:
        case "download-refseq":
            run_download(output_dir, update_taxdump=args.update_taxdump)
        case "download-assembly":
            run_download_assembly(output_dir, batch_size=args.batch_size)
        case "gen_meta":
            run_gen_meta(output_dir, taxdump_dir or output_dir / "taxdump")
        case "gen_db":
            run_gen_db(output_dir)
        case "run":
            run_all(output_dir, taxdump_dir, update_taxdump=args.update_taxdump)


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

    if args.command == "update_refseq":
        _dispatch_update_refseq(args)
        return

    if args.command == "run":
        _dispatch_run(args)
        return

    if args.command == "blast":
        _dispatch_blast(args)
        return

    if args.command == "visualize":
        _dispatch_visualize(args)
        return

    raise NotImplementedError(f"Subcommand '{args.command}' is not implemented")
