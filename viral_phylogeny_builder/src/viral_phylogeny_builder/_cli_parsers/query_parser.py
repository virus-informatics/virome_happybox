# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""argparse definitions for the query command."""

from __future__ import annotations

import argparse
from pathlib import Path

from viral_phylogeny_builder import _commands
from viral_phylogeny_builder._cli_parsers import options

_DEFAULT_BLASTP_EVALUE = 1e-5
_DEFAULT_GAP_THRESHOLD_PCT = 70


def _add_blastp_evalue_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--blastp-evalue",
        type=float,
        default=_DEFAULT_BLASTP_EVALUE,
        help=f"BLASTP E-value threshold (default: {_DEFAULT_BLASTP_EVALUE})",
    )


def _add_gap_threshold_pct_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--gap-threshold-pct",
        type=int,
        default=_DEFAULT_GAP_THRESHOLD_PCT,
        help=(
            f"Maximum gap percentage for sequence filtering (default: {_DEFAULT_GAP_THRESHOLD_PCT})"
        ),
    )


def add_subcommands(parent: argparse._SubParsersAction) -> None:
    p_cmd = parent.add_parser(
        _commands.QUERY_COMMAND,
        help="Add query sequences to reference tree",
    )
    sub = p_cmd.add_subparsers(dest="query_command")

    # search
    p_s = sub.add_parser(
        _commands.QUERY_SEARCH_SUBCOMMAND,
        help="RdRp homology search for query sequences (getorf + BLASTP)",
    )
    options.add_virus_family_option(p_s)
    p_s.add_argument(
        "--input-fasta",
        type=Path,
        required=True,
        help="Query nucleotide FASTA",
    )
    p_s.add_argument(
        "--rdrp-db",
        type=Path,
        required=True,
        help="RdRp protein DB FASTA",
    )
    options.add_output_dir_option(p_s)
    options.add_threads_option(p_s)
    options.add_force_option(p_s)
    options.add_tool_config_option(p_s)
    _add_blastp_evalue_option(p_s)

    # align
    p_a = sub.add_parser(
        _commands.QUERY_ALIGN_SUBCOMMAND,
        help="Add query RdRp to reference alignment (mafft --add + trimal + gap filter)",
    )
    options.add_virus_family_option(p_a)
    p_a.add_argument(
        "--query-rdrp-fasta",
        type=Path,
        required=True,
        help="Query RdRp amino acid FASTA (from search)",
    )
    p_a.add_argument(
        "--reference-alignment",
        type=Path,
        required=True,
        help="Reference MSA FASTA (from build-reftree align)",
    )
    options.add_output_dir_option(p_a)
    options.add_threads_option(p_a)
    options.add_force_option(p_a)
    options.add_tool_config_option(p_a)
    _add_gap_threshold_pct_option(p_a)

    # visualize
    p_v = sub.add_parser(
        _commands.QUERY_VISUALIZE_SUBCOMMAND,
        help="Build tree (IQ-TREE) and render with query + reference metadata",
    )
    options.add_virus_family_option(p_v)
    p_v.add_argument(
        "--alignment-fasta",
        type=Path,
        required=True,
        help="Final alignment FASTA (from align)",
    )
    p_v.add_argument(
        "--name-map-tsv",
        type=Path,
        required=True,
        help="Name map TSV (node/contig/orf_id, from search)",
    )
    p_v.add_argument(
        "--host-csv",
        type=Path,
        required=True,
        help="Per-contig host annotation CSV (user-provided)",
    )
    p_v.add_argument(
        "--reference-taxon-csv",
        type=Path,
        required=True,
        help="Reference taxon metadata CSV (from build-reftree)",
    )
    options.add_output_dir_option(p_v)
    options.add_threads_option(p_v)
    options.add_force_option(p_v)
    options.add_midpoint_rooting_option(p_v)
    options.add_tool_config_option(p_v)

    # run
    p_run = sub.add_parser(_commands.QUERY_RUN_SUBCOMMAND, help="Run full query pipeline")
    options.add_virus_family_option(p_run)
    p_run.add_argument(
        "--input-fasta",
        type=Path,
        required=True,
        help="Query nucleotide FASTA",
    )
    p_run.add_argument(
        "--host-csv",
        type=Path,
        required=True,
        help="Query sample host metadata CSV",
    )
    p_run.add_argument(
        "--reference-alignment",
        type=Path,
        required=True,
        help="Reference MSA FASTA (from build-reftree align)",
    )
    p_run.add_argument(
        "--reference-taxon-csv",
        type=Path,
        required=True,
        help="Reference taxon metadata CSV (from build-reftree)",
    )
    p_run.add_argument(
        "--rdrp-db",
        type=Path,
        required=True,
        help="Curated RdRp protein DB FASTA",
    )
    options.add_output_dir_option(p_run)
    options.add_threads_option(p_run)
    options.add_force_option(p_run)
    options.add_midpoint_rooting_option(p_run)
    options.add_tool_config_option(p_run)
    _add_blastp_evalue_option(p_run)
    _add_gap_threshold_pct_option(p_run)
