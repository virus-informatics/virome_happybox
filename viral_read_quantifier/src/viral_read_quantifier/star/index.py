# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
import math
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

_STAR_GENOME_CHR_BIN_NBITS_DEFAULT = 8

_STAR_GENOME_SA_INDEX_NBASES_MAX = 14
_STAR_GENOME_SA_INDEX_NBASES_MIN = 1


def compute_genome_sa_index_nbases(genome_length_bases: int) -> int:
    """Compute STAR --genomeSAindexNbases from total genome length.

    Formula: min(14, floor(log2(GenomeLength) / 2 - 1))
    See https://github.com/alexdobin/STAR/blob/master/source/parametersDefault
    """
    if genome_length_bases < 1:
        raise ValueError(f"genome_length_bases must be >= 1, got {genome_length_bases}")
    raw = math.log2(genome_length_bases) / 2 - 1
    return max(
        _STAR_GENOME_SA_INDEX_NBASES_MIN,
        min(_STAR_GENOME_SA_INDEX_NBASES_MAX, int(raw)),
    )


def read_fasta_total_length_bases(fasta_path: Path) -> int:
    """Sum all sequence bases in a FASTA file (excluding header lines)."""
    total = 0
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def build_star_index(
    *,
    star_bin: Path,
    genome_fasta: Path,
    index_dir: Path,
    threads: int,
    genome_sa_index_nbases: int,
    genome_chr_bin_nbits: int = _STAR_GENOME_CHR_BIN_NBITS_DEFAULT,
) -> None:
    index_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(star_bin),
        "--runThreadN",
        str(threads),
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        str(index_dir),
        "--genomeFastaFiles",
        str(genome_fasta),
        "--genomeSAindexNbases",
        str(genome_sa_index_nbases),
        "--genomeChrBinNbits",
        str(genome_chr_bin_nbits),
    ]

    logger.info("Running STAR genomeGenerate: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("STAR index created: %s", index_dir)
