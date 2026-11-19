# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""BLASTP wrapper shared by query pipelines (makeblastdb + blastp)."""

from __future__ import annotations

import re
from pathlib import Path

from viral_phylogeny_builder._pipeline import run_subprocess

BLASTP_REQUIRED_TOOLS = ("makeblastdb", "blastp")

BLASTP_OUTFMT_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

_BLASTP_OUTFMT = f"6 {' '.join(BLASTP_OUTFMT_COLUMNS)}"

_MAX_TARGET_SEQS = 10


def sanitize_evalue_tag(evalue: float) -> str:
    """Replace non-alphanumeric characters with underscores for use in filenames."""
    return re.sub(r"[^A-Za-z0-9]", "_", str(evalue))


def run_makeblastdb(
    input_fasta: Path,
    db_prefix: Path,
    *,
    makeblastdb_cmd: str = "makeblastdb",
) -> None:
    run_subprocess(
        [
            makeblastdb_cmd,
            "-in",
            str(input_fasta),
            "-dbtype",
            "prot",
            "-out",
            str(db_prefix),
        ],
        "makeblastdb",
    )


def run_blastp(
    query_fasta: Path,
    db_prefix: Path,
    output_tsv: Path,
    *,
    evalue: float,
    threads: int,
    blastp_cmd: str = "blastp",
) -> None:
    run_subprocess(
        [
            blastp_cmd,
            "-query",
            str(query_fasta),
            "-db",
            str(db_prefix),
            "-evalue",
            str(evalue),
            "-outfmt",
            _BLASTP_OUTFMT,
            "-max_target_seqs",
            str(_MAX_TARGET_SEQS),
            "-num_threads",
            str(threads),
            "-out",
            str(output_tsv),
        ],
        "blastp",
    )
