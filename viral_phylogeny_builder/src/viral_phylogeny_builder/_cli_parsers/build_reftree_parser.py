# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""argparse definitions for the build-reftree command."""

from __future__ import annotations

import argparse
from pathlib import Path

from viral_phylogeny_builder import _commands
from viral_phylogeny_builder._cli_parsers import options


def add_subcommands(parent: argparse._SubParsersAction) -> None:
    p_cmd = parent.add_parser(
        _commands.BUILD_REFTREE_COMMAND,
        help="Build reference phylogenetic tree from RdRp sequences",
    )
    sub = p_cmd.add_subparsers(dest="build_reftree_command")

    # search
    p_s = sub.add_parser(
        _commands.BUILD_REFTREE_SEARCH_SUBCOMMAND,
        help="Extract RdRp sequences via getorf + MMseqs2",
    )
    options.add_virus_family_option(p_s)
    p_s.add_argument(
        "--genome-fasta",
        type=Path,
        required=True,
        help="Reference genome FASTA",
    )
    p_s.add_argument(
        "--reference-metadata-csv",
        type=Path,
        required=True,
        help="Reference metadata CSV (from update-reference fetch-metadata)",
    )
    p_s.add_argument(
        "--rdrp-db",
        type=Path,
        required=True,
        help="Curated RdRp reference FASTA",
    )
    options.add_output_dir_option(p_s)
    options.add_threads_option(p_s)
    options.add_force_option(p_s)
    options.add_tool_config_option(p_s)

    # align
    p_a = sub.add_parser(
        _commands.BUILD_REFTREE_ALIGN_SUBCOMMAND,
        help="Multiple sequence alignment via mafft + trimal",
    )
    options.add_virus_family_option(p_a)
    p_a.add_argument(
        "--rdrp-fasta",
        type=Path,
        required=True,
        help="Extracted RdRp FASTA (from search)",
    )
    p_a.add_argument(
        "--rdrp-metadata-csv",
        type=Path,
        required=True,
        help="RdRp metadata CSV (from search)",
    )
    options.add_output_dir_option(p_a)
    options.add_threads_option(p_a)
    options.add_force_option(p_a)
    options.add_tool_config_option(p_a)

    # visualize
    p_v = sub.add_parser(
        _commands.BUILD_REFTREE_VISUALIZE_SUBCOMMAND,
        help="Build phylogenetic tree (IQ-TREE) and render visualization",
    )
    options.add_virus_family_option(p_v)
    p_v.add_argument(
        "--alignment-fasta",
        type=Path,
        required=True,
        help="Trimmed alignment FASTA (from align)",
    )
    p_v.add_argument(
        "--rdrp-metadata-csv",
        type=Path,
        required=True,
        help="RdRp metadata CSV (from search)",
    )
    options.add_output_dir_option(p_v)
    options.add_threads_option(p_v)
    options.add_force_option(p_v)
    options.add_midpoint_rooting_option(p_v)
    options.add_tool_config_option(p_v)

    # run
    p_run = sub.add_parser(
        _commands.BUILD_REFTREE_RUN_SUBCOMMAND,
        help="Run full build-reftree pipeline",
    )
    options.add_virus_family_option(p_run)
    p_run.add_argument(
        "--genome-fasta",
        type=Path,
        required=True,
        help="Reference genome FASTA",
    )
    p_run.add_argument(
        "--reference-metadata-csv",
        type=Path,
        required=True,
        help="Reference metadata CSV (from update-reference fetch-metadata)",
    )
    p_run.add_argument(
        "--rdrp-db",
        type=Path,
        required=True,
        help="Curated RdRp reference FASTA",
    )
    options.add_output_dir_option(p_run)
    options.add_threads_option(p_run)
    options.add_force_option(p_run)
    options.add_midpoint_rooting_option(p_run)
    options.add_tool_config_option(p_run)
