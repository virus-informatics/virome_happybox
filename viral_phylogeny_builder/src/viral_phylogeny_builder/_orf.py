# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""ORF extraction wrapper shared by build-reftree and query pipelines."""

from __future__ import annotations

from pathlib import Path

from viral_phylogeny_builder._pipeline import run_subprocess

ORF_REQUIRED_TOOLS = ("getorf",)

# getorf measures -minsize on the nucleotide ORF, not on the translated protein. Its own
# default is 30, which is far below the length of an RdRp.
GETORF_MIN_ORF_NUCLEOTIDES = 500


def run_getorf(input_fasta: Path, output_fasta: Path, *, getorf_cmd: str = "getorf") -> None:
    """Translate every ORF of input_fasta into output_fasta as amino acid sequences.

    The record IDs getorf writes are "{input record ID}_{1-based ORF number}".
    """
    run_subprocess(
        [
            getorf_cmd,
            "-sequence",
            str(input_fasta),
            "-outseq",
            str(output_fasta),
            "-minsize",
            str(GETORF_MIN_ORF_NUCLEOTIDES),
        ],
        "getorf",
    )
