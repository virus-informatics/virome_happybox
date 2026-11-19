# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Build the lineage metadata table from the fetched NCBI virus genome metadata."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from viral_phylogeny_builder._accession import strip_accession_version
from viral_phylogeny_builder.update_reference.ncbi_datasets import (
    NCBI_ACCESSION_COLUMN,
    NCBI_COMPLETENESS_COLUMN,
    NCBI_HOST_NAME_COLUMN,
    NCBI_HOST_TAXID_COLUMN,
    NCBI_LENGTH_COLUMN,
    NCBI_VIRUS_NAME_COLUMN,
    NCBI_VIRUS_TAXID_COLUMN,
    VIRUS_GENOME_METADATA_TSV_HEADER,
)

logger = logging.getLogger(__name__)

# The metadata CSV holds every record NCBI has for the family, while the rows that survive
# the filter are a small part of it, so the chunks are filtered as they are read.
_METADATA_READ_CHUNK_ROWS = 100_000


def read_filtered_metadata(metadata_csv: Path, *, exclude_virus_taxid: int) -> pd.DataFrame:
    """Read the metadata CSV and return the rows that pass the filter, with derived columns.

    Raises ValueError when the CSV lacks a required column or when no row survives.
    """
    exclude_virus_taxid_text = str(exclude_virus_taxid)
    kept_chunks: list[pd.DataFrame] = []
    read_row_count = 0

    for chunk in pd.read_csv(metadata_csv, dtype=str, chunksize=_METADATA_READ_CHUNK_ROWS):
        # dtype=str leaves an empty cell as NaN, which none of the derivations accept
        chunk = chunk.fillna("")

        if read_row_count == 0:
            _validate_source_columns(chunk, metadata_csv)
        read_row_count += len(chunk)

        # Some of these only repeat a source column under a snake_case name; they are kept
        # because the port leaves the columns of the archive output unchanged.
        chunk["accession_base"] = chunk[NCBI_ACCESSION_COLUMN].map(strip_accession_version)
        chunk["virus_taxid"] = chunk[NCBI_VIRUS_TAXID_COLUMN].str.strip()
        chunk["virus_name_for_output"] = chunk[NCBI_VIRUS_NAME_COLUMN]
        chunk["host_taxid"] = chunk[NCBI_HOST_TAXID_COLUMN].str.strip()
        chunk["host_name_for_output"] = chunk[NCBI_HOST_NAME_COLUMN]
        chunk["length_numeric"] = pd.to_numeric(chunk[NCBI_LENGTH_COLUMN], errors="coerce")
        chunk["completeness_for_output"] = chunk[NCBI_COMPLETENESS_COLUMN]
        chunk["completeness_upper"] = chunk["completeness_for_output"].str.strip().str.upper()

        chunk = chunk.loc[chunk["accession_base"] != ""]
        chunk = chunk.loc[chunk["virus_taxid"] != ""]
        chunk = chunk.loc[chunk["virus_taxid"] != exclude_virus_taxid_text]
        kept_chunks.append(chunk.loc[chunk["length_numeric"].notna()])

    metadata = (
        pd.concat(kept_chunks, ignore_index=True)
        if kept_chunks
        else pd.DataFrame(columns=list(VIRUS_GENOME_METADATA_TSV_HEADER))
    )
    if metadata.empty:
        raise ValueError(
            f"No metadata row survived the row filter of {metadata_csv}"
            f" (read {read_row_count} rows, excluded virus taxid {exclude_virus_taxid})",
        )

    logger.info("Input rows after basic filtering: %d of %d", len(metadata), read_row_count)
    return metadata


def _validate_source_columns(chunk: pd.DataFrame, metadata_csv: Path) -> None:
    missing_columns = [
        column for column in VIRUS_GENOME_METADATA_TSV_HEADER if column not in chunk.columns
    ]
    if missing_columns:
        raise ValueError(
            f"Metadata CSV is missing columns: {missing_columns} (in {metadata_csv})",
        )
