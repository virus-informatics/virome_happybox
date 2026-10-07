# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Run tblastx against a BLAST nucleotide DB."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

_OUTPUT_FILENAME_MAX_BYTES = 255

_TBLASTX_EVALUE_THRESHOLD = "1e-4"
_TBLASTX_MAX_TARGET_SEQS = 10
_TBLASTX_DEFAULT_NUM_THREADS = 4

_TBLASTX_OUTFMT_COLUMNS = (
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


def infer_sra_id_from_fasta(fasta_path: Path) -> str:
    """Read the first FASTA header and extract the SRA ID prefix (before '_')."""
    with fasta_path.open() as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                prefix = header.split("_")[0]
                if prefix:
                    return prefix
    raise ValueError(f"Cannot infer SRA ID from {fasta_path}: no FASTA header found")


def _build_tblastx_command(
    query_fasta: Path,
    db_path: Path,
    output_tsv: Path,
    num_threads: int,
) -> list[str]:
    outfmt = "6 " + " ".join(_TBLASTX_OUTFMT_COLUMNS)
    return [
        "tblastx",
        "-query",
        str(query_fasta),
        "-db",
        str(db_path),
        "-evalue",
        _TBLASTX_EVALUE_THRESHOLD,
        "-max_target_seqs",
        str(_TBLASTX_MAX_TARGET_SEQS),
        "-num_threads",
        str(num_threads),
        "-outfmt",
        outfmt,
        "-out",
        str(output_tsv),
    ]


def run_tblastx(
    query_fasta: Path,
    db_path: Path,
    workdir: Path,
    *,
    output_tsv: Path | None = None,
    num_threads: int = _TBLASTX_DEFAULT_NUM_THREADS,
    force: bool = False,
) -> Path:
    if not query_fasta.exists():
        raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")
    if query_fasta.stat().st_size == 0:
        raise ValueError(f"Query FASTA is empty: {query_fasta}")

    tblastx_dir = workdir / "tblastx"
    tblastx_dir.mkdir(parents=True, exist_ok=True)

    if output_tsv is None:
        stem = infer_sra_id_from_fasta(query_fasta)
        logger.info("Inferred output name from FASTA header: %s", stem)
        output_tsv = tblastx_dir / f"{stem}.tblastx.tsv"

    filename_bytes = len(output_tsv.name.encode())
    if filename_bytes > _OUTPUT_FILENAME_MAX_BYTES:
        raise ValueError(
            f"Output filename too long ({filename_bytes} bytes,"
            f" max {_OUTPUT_FILENAME_MAX_BYTES}): {output_tsv.name}",
        )

    if output_tsv.exists() and output_tsv.stat().st_size > 0:
        if not force:
            logger.info("Skipping tblastx (cached): %s", output_tsv)
            return output_tsv
        logger.warning("Overwriting existing output: %s", output_tsv)

    cmd = _build_tblastx_command(query_fasta, db_path, output_tsv, num_threads)
    logger.info("Running tblastx: %s -> %s", query_fasta.name, output_tsv)

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            cmd,
            result.stdout,
            result.stderr,
        )

    if result.stderr:
        logger.debug("tblastx stderr:\n%s", result.stderr.rstrip())

    logger.info("tblastx complete: %s", output_tsv)
    return output_tsv
