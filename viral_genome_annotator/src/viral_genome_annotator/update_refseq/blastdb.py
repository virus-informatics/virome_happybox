# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Build a BLAST nucleotide DB from a gzipped FASTA via makeblastdb."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

_BLASTDB_NAME = "refseq_viral_nucl"
_BLASTDB_TITLE = "RefSeq Viral Nucleotide"

# makeblastdb generates multiple sidecar files (.ndb, .nhr, .nin, ...).
# The stem (dbtype=nucl) always starts with ".n".
_BLASTDB_NUCL_EXTENSIONS = (
    ".ndb",
    ".nhr",
    ".nin",
    ".njs",
    ".not",
    ".nsq",
    ".ntf",
    ".nto",
)


def _run_makeblastdb(fna_gz_path: Path, db_path: Path) -> subprocess.CompletedProcess[str]:
    """Pipe gzip-decompressed FASTA into makeblastdb via stdin."""
    decompress = subprocess.Popen(
        ["gzip", "-dc", str(fna_gz_path)],
        stdout=subprocess.PIPE,
    )
    result = subprocess.run(
        [
            "makeblastdb",
            "-in",
            "-",
            "-dbtype",
            "nucl",
            "-out",
            str(db_path),
            "-title",
            _BLASTDB_TITLE,
        ],
        stdin=decompress.stdout,
        capture_output=True,
        text=True,
    )
    if decompress.stdout is not None:
        decompress.stdout.close()
    decompress.wait()

    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            result.args,
            result.stdout,
            result.stderr,
        )
    if decompress.returncode != 0:
        raise subprocess.SubprocessError(f"gzip -dc failed with exit code {decompress.returncode}")
    return result


def build_blastdb(fna_gz_path: Path, outdir: Path) -> Path:
    """Build a BLAST nucleotide DB from a gzipped FASTA.

    Returns the DB path stem (pass this to blastn -db).
    """
    if not fna_gz_path.exists():
        raise FileNotFoundError(f"Input file not found: {fna_gz_path}")

    outdir.mkdir(parents=True, exist_ok=True)
    db_path = outdir / _BLASTDB_NAME

    logger.info("Building BLAST DB from %s -> %s", fna_gz_path.name, db_path)
    result = _run_makeblastdb(fna_gz_path, db_path)

    if result.stderr:
        logger.debug("makeblastdb stderr:\n%s", result.stderr.rstrip())

    generated = [f for ext in _BLASTDB_NUCL_EXTENSIONS if (f := db_path.with_suffix(ext)).exists()]
    if not generated:
        raise RuntimeError(
            f"makeblastdb produced no output files at {db_path}. stdout: {result.stdout.rstrip()}",
        )

    logger.info("BLAST DB ready: %s (%d files)", db_path, len(generated))
    return db_path
