# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""update_refseq download-assembly: fetch GCF→NC mapping from NCBI Datasets API."""

from __future__ import annotations

import logging
from pathlib import Path

from viral_genome_annotator.update_refseq.assembly import (
    _DEFAULT_ASSEMBLY_BATCH_SIZE,
    download_assembly_map,
)

logger = logging.getLogger(__name__)


def run_download_assembly(
    output_dir: Path,
    batch_size: int = _DEFAULT_ASSEMBLY_BATCH_SIZE,
) -> None:
    map_path = download_assembly_map(output_dir, batch_size=batch_size)
    logger.info("Assembly map ready: %s", map_path)
