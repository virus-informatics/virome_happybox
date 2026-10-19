# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""update_refseq gen_db: build BLAST DB from downloaded FASTA."""

from __future__ import annotations

import logging
from pathlib import Path

from viral_genome_annotator.update_refseq.blastdb import build_blastdb
from viral_genome_annotator.update_refseq.download import NCBI_DOWNLOAD_FILES

logger = logging.getLogger(__name__)

_FNA_GZ_FILENAME = next(k for k in NCBI_DOWNLOAD_FILES if k.endswith(".fna.gz"))
_BLASTDB_SUBDIR = "blast_db"


def run_gen_db(output_dir: Path) -> None:
    from viral_genome_annotator.external import require_command

    require_command("makeblastdb")

    downloads_dir = output_dir / "downloads"
    fna_gz_path = downloads_dir / _FNA_GZ_FILENAME
    if not fna_gz_path.exists():
        raise FileNotFoundError(
            f"FASTA file not found: {fna_gz_path}. Run 'update_refseq download' first.",
        )

    blastdb_dir = output_dir / _BLASTDB_SUBDIR
    db_path = build_blastdb(fna_gz_path, blastdb_dir)
    logger.info("BLAST DB ready: %s", db_path)
