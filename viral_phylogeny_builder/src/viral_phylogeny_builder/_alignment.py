# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Alignment tool wrappers shared by build-reftree and query pipelines."""

from __future__ import annotations

from pathlib import Path

from viral_phylogeny_builder._pipeline import run_subprocess

ALIGNMENT_REQUIRED_TOOLS = ("mafft", "trimal", "seqkit")


def run_trimal(input_fasta: Path, output_fasta: Path, *, trimal_cmd: str = "trimal") -> None:
    run_subprocess(
        [trimal_cmd, "-in", str(input_fasta), "-out", str(output_fasta), "-gappyout"],
        "trimal -gappyout",
    )


def run_seqkit_rmdup(input_fasta: Path, output_fasta: Path, *, seqkit_cmd: str = "seqkit") -> None:
    run_subprocess(
        [seqkit_cmd, "rmdup", "-n", str(input_fasta), "-o", str(output_fasta)],
        "seqkit rmdup",
    )
