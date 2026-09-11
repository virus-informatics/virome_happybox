# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
from pathlib import Path

from rnaseq_mapper.star.binary import check_star_version, resolve_star_binary
from rnaseq_mapper.star.index import (
    build_star_index,
    compute_genome_sa_index_nbases,
    read_fasta_total_length_bases,
)

logger = logging.getLogger(__name__)

_INDEX_COMPLETE_FILES = ("Genome", "SA")


def _genome_name(genome_fasta: Path) -> str:
    return genome_fasta.stem


def _index_dir_for(output_dir: Path, genome_fasta: Path) -> Path:
    return output_dir / "star_index" / _genome_name(genome_fasta)


def _index_exists(index_dir: Path) -> bool:
    return all(
        (p := index_dir / f).is_file() and p.stat().st_size > 0 for f in _INDEX_COMPLETE_FILES
    )


def run_index(
    *,
    genome_fasta: Path,
    output_dir: Path,
    star_bin: Path | None,
    threads: int,
    force: bool,
    genome_sa_index_nbases: int | None = None,
    genome_chr_bin_nbits: int | None = None,
) -> Path:
    """Build STAR genome index. Returns the index directory path."""
    if not genome_fasta.is_file():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fasta}")

    resolved_bin = resolve_star_binary(star_bin)
    check_star_version(resolved_bin)

    index_dir = _index_dir_for(output_dir, genome_fasta)

    if _index_exists(index_dir) and not force:
        logger.info("STAR index already exists, skipping: %s", index_dir)
        return index_dir

    if genome_sa_index_nbases is None:
        genome_length = read_fasta_total_length_bases(genome_fasta)
        genome_sa_index_nbases = compute_genome_sa_index_nbases(genome_length)
        logger.info(
            "Auto-computed genomeSAindexNbases=%d from genome length %d bp",
            genome_sa_index_nbases,
            genome_length,
        )

    extra_kwargs: dict[str, int] = {}
    if genome_chr_bin_nbits is not None:
        extra_kwargs["genome_chr_bin_nbits"] = genome_chr_bin_nbits

    build_star_index(
        star_bin=resolved_bin,
        genome_fasta=genome_fasta,
        index_dir=index_dir,
        threads=threads,
        genome_sa_index_nbases=genome_sa_index_nbases,
        **extra_kwargs,
    )

    return index_dir
