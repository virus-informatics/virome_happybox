# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Concatenate reference and added-sequence taxonomy metadata."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from viral_phylogeny_builder._fasta import iter_fasta

logger = logging.getLogger(__name__)

CONCAT_METADATA_FILENAME_SUFFIX = ".concat_metadata.csv"

FINAL_COLUMNS = [
    "node",
    "Accession",
    "Genus",
    "manual_host",
    "Organism_Name",
    "GenBank_RefSeq",
    "VirusFamily",
    "Host",
    "Lineage",
]


def concat_metadata_filename(virus_family: str) -> str:
    return f"{virus_family}{CONCAT_METADATA_FILENAME_SUFFIX}"


def _label_key(value: str) -> str:
    return value.strip().split()[0].split(".")[0]


def _read_valid_nodes(alignment_fasta: Path) -> set[str]:
    """Read FASTA headers and return versionless IDs as the valid node set."""
    nodes: set[str] = set()
    for header, _seq in iter_fasta(alignment_fasta):
        value = _label_key(header)
        if value:
            nodes.add(value)
    return nodes


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def _normalized_row(
    row: dict[str, str],
    valid_nodes: set[str],
    is_ref: bool,
) -> dict[str, str] | None:
    node = row.get("node") or row.get("Accession") or ""
    node = _label_key(node)
    if valid_nodes and node not in valid_nodes:
        return None

    output = {column: row.get(column, "") for column in FINAL_COLUMNS}
    output["node"] = node
    output["Accession"] = _label_key(output.get("Accession", node)) or node
    if is_ref:
        output["VirusFamily"] = output.get("VirusFamily", "")
    else:
        output["Genus"] = output.get("Genus", "")
    return output


def concat_metadata(
    ref_csv: Path,
    add_csv: Path,
    alignment_fasta: Path,
    output_csv: Path,
) -> None:
    """Concatenate ref and add metadata, filtering to valid sequences in the alignment."""
    valid_nodes = _read_valid_nodes(alignment_fasta)

    rows: list[dict[str, str]] = []
    for row in _read_rows(ref_csv):
        normalized = _normalized_row(row, valid_nodes, is_ref=True)
        if normalized is not None:
            rows.append(normalized)

    for row in _read_rows(add_csv):
        normalized = _normalized_row(row, valid_nodes, is_ref=False)
        if normalized is not None:
            rows.append(normalized)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FINAL_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Combined metadata: %d rows -> %s", len(rows), output_csv)
