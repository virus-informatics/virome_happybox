# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""update_refseq gen_meta: parse GBFF → GFF3 + taxonomy.tsv."""

from __future__ import annotations

import json
import logging
from collections.abc import Iterator
from pathlib import Path

from viral_genome_annotator.update_refseq.assembly import ASSEMBLY_MAP_FILENAME, ASSEMBLY_SUBDIR
from viral_genome_annotator.update_refseq.download import NCBI_DOWNLOAD_FILES
from viral_genome_annotator.update_refseq.gbff import GbffRecord, parse_gbff_gzipped
from viral_genome_annotator.update_refseq.gff3 import write_gff3_file
from viral_genome_annotator.update_refseq.taxonomy import (
    resolve_taxonomy,
    write_taxonomy_tsv,
)

logger = logging.getLogger(__name__)

_GBFF_FILENAME = next(k for k in NCBI_DOWNLOAD_FILES if k.endswith(".gbff.gz"))
_GFF3_OUTPUT_FILENAME = "annotations.gff3"
_TAXONOMY_OUTPUT_FILENAME = "taxonomy.tsv"


def _tap_taxids(
    records: Iterator[GbffRecord],
    taxids: set[int],
) -> Iterator[GbffRecord]:
    for record in records:
        if record.taxid is not None:
            taxids.add(record.taxid)
        yield record


def _read_assembly_taxids(assembly_map_path: Path) -> set[int]:
    """Read the Assembly-level taxids recorded in assembly_map.json.

    NCBI assigns these per assembly, and they do not always appear among the
    per-accession GBFF taxids, so resolving from the GBFF alone would leave them out
    of taxonomy.tsv.
    """
    mappings = json.loads(assembly_map_path.read_text()).get("mappings", {})
    taxids = {entry["taxid"] for entry in mappings.values() if entry.get("taxid") is not None}
    logger.info("Collected %d assembly-level taxids from %s", len(taxids), assembly_map_path)
    return taxids


def run_gen_meta(output_dir: Path, taxdump_dir: Path) -> None:
    from viral_genome_annotator.external import require_command

    require_command("taxonkit")

    downloads_dir = output_dir / "downloads"
    gbff_path = downloads_dir / _GBFF_FILENAME
    if not gbff_path.exists():
        raise FileNotFoundError(
            f"GBFF file not found: {gbff_path}. Run 'update_refseq download' first.",
        )

    # Checked before parsing the GBFF so the run does not fail after the expensive pass
    assembly_map_path = output_dir / ASSEMBLY_SUBDIR / ASSEMBLY_MAP_FILENAME
    if not assembly_map_path.exists():
        raise FileNotFoundError(
            f"Assembly map not found: {assembly_map_path}."
            " Run 'update_refseq download-assembly' first.",
        )

    logger.info("Parsing GBFF: %s", gbff_path)
    taxids: set[int] = set()
    records = _tap_taxids(parse_gbff_gzipped(gbff_path), taxids)

    gff3_path = output_dir / _GFF3_OUTPUT_FILENAME
    cds_count = write_gff3_file(records, gff3_path)
    logger.info("Wrote %d CDS lines to %s", cds_count, gff3_path)

    logger.info("Collected %d unique taxids from GBFF", len(taxids))
    taxids |= _read_assembly_taxids(assembly_map_path)

    logger.info("Resolving %d unique taxids", len(taxids))
    entries = resolve_taxonomy(taxids, taxdump_dir)
    taxonomy_path = output_dir / _TAXONOMY_OUTPUT_FILENAME
    write_taxonomy_tsv(entries, taxonomy_path)
