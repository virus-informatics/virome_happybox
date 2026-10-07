# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""update_refseq download-refseq: fetch reference files from NCBI FTP."""

from __future__ import annotations

import logging
from pathlib import Path

from viral_genome_annotator.update_refseq.download import (
    download_refseq,
    download_taxdump,
    has_taxdump,
)

logger = logging.getLogger(__name__)


def run_download(output_dir: Path, *, update_taxdump: bool = False) -> None:
    downloads_dir = download_refseq(output_dir)
    logger.info("Downloads complete: %s", downloads_dir)

    if update_taxdump or not has_taxdump(output_dir):
        taxdump_dir = download_taxdump(output_dir)
        logger.info("Taxdump ready: %s", taxdump_dir)
    else:
        logger.info("Taxdump update skipped (use --update-taxdump to update)")
