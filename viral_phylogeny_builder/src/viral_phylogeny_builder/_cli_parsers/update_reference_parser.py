# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""argparse definitions for the update-reference command."""

from __future__ import annotations

import argparse
from pathlib import Path

from viral_phylogeny_builder import _commands
from viral_phylogeny_builder._cli_parsers import options
from viral_phylogeny_builder.update_reference.taxdump import TAXDUMP_DIR_NAME

# NCBI registers far more SARS-CoV-2 records than any other Coronaviridae species,
# which would otherwise swamp representative-sequence selection.
_DEFAULT_EXCLUDE_VIRUS_TAXID = 2697049  # SARS-CoV-2


def _add_exclude_virus_taxid_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--exclude-virus-taxid",
        type=int,
        default=_DEFAULT_EXCLUDE_VIRUS_TAXID,
        help=f"NCBI taxonomy ID to exclude (default: {_DEFAULT_EXCLUDE_VIRUS_TAXID}, SARS-CoV-2)",
    )


def _add_update_taxdump_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--update-taxdump",
        action="store_true",
        help="Also download/update NCBI taxdump (skipped by default)",
    )


def _add_taxdump_dir_option(parser: argparse.ArgumentParser) -> None:
    # Left unset by default so that the pipeline can tell a supplied directory from the
    # managed one: the managed one is downloaded when missing or on --update-taxdump, a
    # supplied one is only checked. Its default value cannot be built here either, as it
    # follows --output-dir.
    parser.add_argument(
        "--taxdump-dir",
        type=Path,
        default=None,
        help=(
            "Directory holding the NCBI taxdump that taxonkit reads"
            f" (default: <output-dir>/{TAXDUMP_DIR_NAME}, downloaded on first run)"
        ),
    )


def add_subcommands(parent: argparse._SubParsersAction) -> None:
    p_cmd = parent.add_parser(
        _commands.UPDATE_REFERENCE_COMMAND,
        help="Fetch NCBI metadata and download reference genomes",
    )
    sub = p_cmd.add_subparsers(dest="update_reference_command")

    # fetch-taxdump
    p_ft = sub.add_parser(
        _commands.UPDATE_REFERENCE_FETCH_TAXDUMP_SUBCOMMAND,
        help=f"Download/update NCBI taxdump into <output-dir>/{TAXDUMP_DIR_NAME}",
    )
    options.add_output_dir_option(p_ft, default=options.DEFAULT_REFERENCE_OUTDIR)
    options.add_tool_config_option(p_ft)

    # fetch-metadata
    p_fm = sub.add_parser(
        _commands.UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
        help="Fetch NCBI virus metadata and select representative sequences",
    )
    options.add_virus_family_option(p_fm)
    options.add_output_dir_option(p_fm, default=options.DEFAULT_REFERENCE_OUTDIR)
    options.add_force_option(p_fm)
    options.add_tool_config_option(p_fm)
    _add_exclude_virus_taxid_option(p_fm)
    _add_taxdump_dir_option(p_fm)
    _add_update_taxdump_option(p_fm)

    # fetch-genome
    p_fg = sub.add_parser(
        _commands.UPDATE_REFERENCE_FETCH_GENOME_SUBCOMMAND,
        help="Download reference genomes using accessions from fetch-metadata",
    )
    options.add_virus_family_option(p_fg)
    options.add_output_dir_option(p_fg, default=options.DEFAULT_REFERENCE_OUTDIR)
    options.add_force_option(p_fg)
    options.add_tool_config_option(p_fg)

    # run
    p_run = sub.add_parser(
        _commands.UPDATE_REFERENCE_RUN_SUBCOMMAND,
        help="Run full update-reference pipeline",
    )
    options.add_virus_family_option(p_run)
    options.add_output_dir_option(p_run, default=options.DEFAULT_REFERENCE_OUTDIR)
    options.add_force_option(p_run)
    options.add_tool_config_option(p_run)
    _add_exclude_virus_taxid_option(p_run)
    _add_taxdump_dir_option(p_run)
    _add_update_taxdump_option(p_run)
