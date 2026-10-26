# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# match's value patterns need a dotted name: a bare constant would be a capture pattern
from viral_phylogeny_builder import _commands
from viral_phylogeny_builder._cli_parsers import (
    build_reftree_parser,
    query_parser,
    update_reference_parser,
)
from viral_phylogeny_builder._tool_config import load_tool_config

logger = logging.getLogger(__name__)

_PROG_NAME = "viral_phylogeny_builder"
_DEFAULT_TOOL_CONFIG_FILENAME = "config.toml"


def _load_tool_config(args: argparse.Namespace) -> dict[str, str]:
    path = args.tool_config
    if path is None:
        default = Path(_DEFAULT_TOOL_CONFIG_FILENAME)
        path = default if default.exists() else None
    return load_tool_config(path)


def _subcommand_usage(command: str, subcommands: tuple[str, ...]) -> str:
    return f"usage: {_PROG_NAME} {command} {{{','.join(subcommands)}}}"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=_PROG_NAME,
        description="Virome analysis — phylogenetic tree visualization",
    )
    sub = parser.add_subparsers(dest="command")

    update_reference_parser.add_subcommands(sub)
    build_reftree_parser.add_subcommands(sub)
    query_parser.add_subcommands(sub)

    return parser


def _run_update_reference_pipeline(args: argparse.Namespace) -> None:
    from viral_phylogeny_builder.update_reference.fetch_genome import (
        run_update_reference_fetch_genome,
    )
    from viral_phylogeny_builder.update_reference.fetch_metadata import (
        run_update_reference_fetch_metadata,
    )

    tool_config = _load_tool_config(args)

    run_update_reference_fetch_metadata(
        virus_family=args.virus_family,
        output_dir=args.output_dir,
        exclude_virus_taxid=args.exclude_virus_taxid,
        taxdump_dir=args.taxdump_dir,
        update_taxdump=args.update_taxdump,
        force=args.force,
        tool_config=tool_config,
    )
    genome_fasta = run_update_reference_fetch_genome(
        virus_family=args.virus_family,
        output_dir=args.output_dir,
        force=args.force,
        tool_config=tool_config,
    )
    logger.info("Update reference complete: %s", genome_fasta)


def _dispatch_update_reference(args: argparse.Namespace) -> None:
    if args.update_reference_command is None:
        raise SystemExit(
            _subcommand_usage(
                _commands.UPDATE_REFERENCE_COMMAND,
                _commands.UPDATE_REFERENCE_SUBCOMMANDS,
            ),
        )

    tool_config = _load_tool_config(args)

    match args.update_reference_command:
        case _commands.UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND:
            from viral_phylogeny_builder.update_reference.fetch_metadata import (
                run_update_reference_fetch_metadata,
            )

            run_update_reference_fetch_metadata(
                virus_family=args.virus_family,
                output_dir=args.output_dir,
                exclude_virus_taxid=args.exclude_virus_taxid,
                taxdump_dir=args.taxdump_dir,
                update_taxdump=args.update_taxdump,
                force=args.force,
                tool_config=tool_config,
            )
        case _commands.UPDATE_REFERENCE_FETCH_GENOME_SUBCOMMAND:
            from viral_phylogeny_builder.update_reference.fetch_genome import (
                run_update_reference_fetch_genome,
            )

            run_update_reference_fetch_genome(
                virus_family=args.virus_family,
                output_dir=args.output_dir,
                force=args.force,
                tool_config=tool_config,
            )
        case _commands.UPDATE_REFERENCE_RUN_SUBCOMMAND:
            _run_update_reference_pipeline(args)


def _run_build_reftree_pipeline(args: argparse.Namespace) -> None:
    from viral_phylogeny_builder._pipeline import step_workdir
    from viral_phylogeny_builder.build_reftree.align import (
        build_reftree_align_final_filename,
        run_build_reftree_align,
    )
    from viral_phylogeny_builder.build_reftree.search import (
        build_reftree_search_rdrp_filename,
        build_reftree_search_rdrp_metadata_filename,
        run_build_reftree_search,
    )
    from viral_phylogeny_builder.build_reftree.visualize import run_build_reftree_visualize

    virus_family = args.virus_family
    tool_config = _load_tool_config(args)

    run_build_reftree_search(
        virus_family=virus_family,
        genome_fasta=args.genome_fasta,
        reference_metadata_csv=args.reference_metadata_csv,
        rdrp_db=args.rdrp_db,
        output_dir=args.output_dir,
        threads=args.threads,
        force=args.force,
        tool_config=tool_config,
    )

    search_workdir = step_workdir(
        args.output_dir,
        _commands.BUILD_REFTREE_COMMAND,
        _commands.BUILD_REFTREE_SEARCH_SUBCOMMAND,
        virus_family,
    )
    rdrp_fasta = search_workdir / build_reftree_search_rdrp_filename(virus_family)
    rdrp_metadata_csv = search_workdir / build_reftree_search_rdrp_metadata_filename(virus_family)

    run_build_reftree_align(
        virus_family=virus_family,
        rdrp_fasta=rdrp_fasta,
        rdrp_metadata_csv=rdrp_metadata_csv,
        output_dir=args.output_dir,
        threads=args.threads,
        force=args.force,
        tool_config=tool_config,
    )

    align_workdir = step_workdir(
        args.output_dir,
        _commands.BUILD_REFTREE_COMMAND,
        _commands.BUILD_REFTREE_ALIGN_SUBCOMMAND,
        virus_family,
    )
    alignment_fasta = align_workdir / build_reftree_align_final_filename(virus_family)

    run_build_reftree_visualize(
        virus_family=virus_family,
        alignment_fasta=alignment_fasta,
        rdrp_metadata_csv=rdrp_metadata_csv,
        output_dir=args.output_dir,
        threads=args.threads,
        force=args.force,
        midpoint_rooting=args.midpoint_rooting,
        tool_config=tool_config,
    )
    logger.info("Build-reftree pipeline complete for %s", virus_family)


def _dispatch_build_reftree(args: argparse.Namespace) -> None:
    if args.build_reftree_command is None:
        raise SystemExit(
            _subcommand_usage(
                _commands.BUILD_REFTREE_COMMAND,
                _commands.BUILD_REFTREE_SUBCOMMANDS,
            ),
        )

    tool_config = _load_tool_config(args)

    match args.build_reftree_command:
        case _commands.BUILD_REFTREE_SEARCH_SUBCOMMAND:
            from viral_phylogeny_builder.build_reftree.search import run_build_reftree_search

            run_build_reftree_search(
                virus_family=args.virus_family,
                genome_fasta=args.genome_fasta,
                reference_metadata_csv=args.reference_metadata_csv,
                rdrp_db=args.rdrp_db,
                output_dir=args.output_dir,
                threads=args.threads,
                force=args.force,
                tool_config=tool_config,
            )
        case _commands.BUILD_REFTREE_ALIGN_SUBCOMMAND:
            from viral_phylogeny_builder.build_reftree.align import run_build_reftree_align

            run_build_reftree_align(
                virus_family=args.virus_family,
                rdrp_fasta=args.rdrp_fasta,
                rdrp_metadata_csv=args.rdrp_metadata_csv,
                output_dir=args.output_dir,
                threads=args.threads,
                force=args.force,
                tool_config=tool_config,
            )
        case _commands.BUILD_REFTREE_VISUALIZE_SUBCOMMAND:
            from viral_phylogeny_builder.build_reftree.visualize import (
                run_build_reftree_visualize,
            )

            run_build_reftree_visualize(
                virus_family=args.virus_family,
                alignment_fasta=args.alignment_fasta,
                rdrp_metadata_csv=args.rdrp_metadata_csv,
                output_dir=args.output_dir,
                threads=args.threads,
                force=args.force,
                midpoint_rooting=args.midpoint_rooting,
                tool_config=tool_config,
            )
        case _commands.BUILD_REFTREE_RUN_SUBCOMMAND:
            _run_build_reftree_pipeline(args)


def _run_query_pipeline(args: argparse.Namespace) -> None:
    from viral_phylogeny_builder.query.align import run_query_align
    from viral_phylogeny_builder.query.search import run_query_search
    from viral_phylogeny_builder.query.visualize import run_query_visualize

    tool_config = _load_tool_config(args)

    renamed_fasta = run_query_search(
        virus_family=args.virus_family,
        input_fasta=args.input_fasta,
        rdrp_db=args.rdrp_db,
        output_dir=args.output_dir,
        blastp_evalue=args.blastp_evalue,
        threads=args.threads,
        force=args.force,
        tool_config=tool_config,
    )

    name_map_tsv = renamed_fasta.with_name(
        renamed_fasta.name.replace(".renamed.fasta", ".name_map.tsv"),
    )

    final_alignment = run_query_align(
        virus_family=args.virus_family,
        query_rdrp_fasta=renamed_fasta,
        reference_alignment=args.reference_alignment,
        output_dir=args.output_dir,
        threads=args.threads,
        force=args.force,
        gap_threshold_pct=args.gap_threshold_pct,
        tool_config=tool_config,
    )

    run_query_visualize(
        virus_family=args.virus_family,
        alignment_fasta=final_alignment,
        name_map_tsv=name_map_tsv,
        host_csv=args.host_csv,
        reference_taxon_csv=args.reference_taxon_csv,
        output_dir=args.output_dir,
        threads=args.threads,
        force=args.force,
        midpoint_rooting=args.midpoint_rooting,
        tool_config=tool_config,
    )


def _dispatch_query(args: argparse.Namespace) -> None:
    if args.query_command is None:
        raise SystemExit(
            _subcommand_usage(_commands.QUERY_COMMAND, _commands.QUERY_SUBCOMMANDS),
        )

    tool_config = _load_tool_config(args)

    match args.query_command:
        case _commands.QUERY_SEARCH_SUBCOMMAND:
            from viral_phylogeny_builder.query.search import run_query_search

            run_query_search(
                virus_family=args.virus_family,
                input_fasta=args.input_fasta,
                rdrp_db=args.rdrp_db,
                output_dir=args.output_dir,
                blastp_evalue=args.blastp_evalue,
                threads=args.threads,
                force=args.force,
                tool_config=tool_config,
            )
        case _commands.QUERY_ALIGN_SUBCOMMAND:
            from viral_phylogeny_builder.query.align import run_query_align

            run_query_align(
                virus_family=args.virus_family,
                query_rdrp_fasta=args.query_rdrp_fasta,
                reference_alignment=args.reference_alignment,
                output_dir=args.output_dir,
                threads=args.threads,
                force=args.force,
                gap_threshold_pct=args.gap_threshold_pct,
                tool_config=tool_config,
            )
        case _commands.QUERY_VISUALIZE_SUBCOMMAND:
            from viral_phylogeny_builder.query.visualize import run_query_visualize

            run_query_visualize(
                virus_family=args.virus_family,
                alignment_fasta=args.alignment_fasta,
                name_map_tsv=args.name_map_tsv,
                host_csv=args.host_csv,
                reference_taxon_csv=args.reference_taxon_csv,
                output_dir=args.output_dir,
                threads=args.threads,
                force=args.force,
                midpoint_rooting=args.midpoint_rooting,
                tool_config=tool_config,
            )
        case _commands.QUERY_RUN_SUBCOMMAND:
            _run_query_pipeline(args)


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

    match args.command:
        case _commands.UPDATE_REFERENCE_COMMAND:
            _dispatch_update_reference(args)
        case _commands.BUILD_REFTREE_COMMAND:
            _dispatch_build_reftree(args)
        case _commands.QUERY_COMMAND:
            _dispatch_query(args)
        case _:
            raise NotImplementedError(f"Subcommand '{args.command}' is not implemented")
