# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""update_refseq run: download → gen_meta + gen_db pipeline."""

from __future__ import annotations

import logging
from pathlib import Path

from viral_genome_annotator.update_refseq.cmd_download_assembly import run_download_assembly
from viral_genome_annotator.update_refseq.cmd_download_refseq import run_download
from viral_genome_annotator.update_refseq.cmd_gen_db import run_gen_db
from viral_genome_annotator.update_refseq.cmd_gen_meta import run_gen_meta

logger = logging.getLogger(__name__)

_TAXDUMP_SUBDIR = "taxdump"


def run_all(
    output_dir: Path,
    taxdump_dir: Path | None = None,
    *,
    update_taxdump: bool = False,
) -> None:
    from viral_genome_annotator.external import require_command

    require_command("taxonkit")
    require_command("makeblastdb")

    logger.info("=== update_refseq run: starting full pipeline ===")

    run_download(output_dir, update_taxdump=update_taxdump)
    run_download_assembly(output_dir)

    resolved_taxdump_dir = taxdump_dir if taxdump_dir is not None else output_dir / _TAXDUMP_SUBDIR
    run_gen_meta(output_dir, resolved_taxdump_dir)
    run_gen_db(output_dir)

    logger.info("=== update_refseq run: pipeline complete ===")
