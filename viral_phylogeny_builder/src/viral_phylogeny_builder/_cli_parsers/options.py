# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Option defaults and the option definitions shared by more than one command."""

from __future__ import annotations

import argparse
from pathlib import Path

_DEFAULT_OUTDIR = Path("output_dir")
_DEFAULT_THREADS = 4


def add_output_dir_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=_DEFAULT_OUTDIR,
        help=f"Output directory (default: {_DEFAULT_OUTDIR})",
    )


def add_threads_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--threads",
        type=int,
        default=_DEFAULT_THREADS,
        help=f"Number of threads (default: {_DEFAULT_THREADS})",
    )


def add_tool_config_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--tool-config",
        type=Path,
        default=None,
        help="TOML file with external tool paths (default: config.toml if present)",
    )


def add_force_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if output already exists",
    )


def add_virus_family_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--virus-family",
        required=True,
        help="Virus family name (e.g. Coronaviridae)",
    )


def add_midpoint_rooting_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--midpoint-rooting",
        action="store_true",
        help="Apply midpoint rooting to the tree before visualization",
    )
