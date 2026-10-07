# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Parse tblastx outfmt-6 TSV (14 columns) and extract stem from filename."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

_TBLASTX_EXPECTED_COLUMN_COUNT = 14


@dataclass(frozen=True, slots=True)
class BlastHit:
    tsv_stem: str
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int


def extract_tsv_stem(tsv_path: Path) -> str:
    """Extract stem from TSV filename.

    '.tblastx.tsv' suffix is stripped if present, otherwise the standard
    Path.stem is used.
    """
    name = tsv_path.name
    tblastx_suffix = ".tblastx.tsv"
    if name.endswith(tblastx_suffix):
        stem = name[: -len(tblastx_suffix)]
    else:
        stem = tsv_path.stem
    if not stem:
        raise ValueError(f"Empty stem in filename: {name}")
    return stem


def _parse_row(fields: list[str], tsv_stem: str) -> BlastHit:
    return BlastHit(
        tsv_stem=tsv_stem,
        qseqid=fields[0],
        sseqid=fields[1],
        pident=float(fields[2]),
        length=int(fields[3]),
        mismatch=int(fields[4]),
        gapopen=int(fields[5]),
        qstart=int(fields[6]),
        qend=int(fields[7]),
        sstart=int(fields[8]),
        send=int(fields[9]),
        evalue=float(fields[10]),
        bitscore=float(fields[11]),
        qlen=int(fields[12]),
        slen=int(fields[13]),
    )


def parse_blast_tsv(tsv_path: Path) -> list[BlastHit]:
    """Parse a tblastx outfmt-6 TSV (14 columns). Stem is extracted from filename."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"BLAST TSV not found: {tsv_path}")
    if tsv_path.stat().st_size == 0:
        logger.warning("BLAST TSV is empty: %s", tsv_path)
        return []

    tsv_stem = extract_tsv_stem(tsv_path)
    hits: list[BlastHit] = []

    with tsv_path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for line_num, row in enumerate(reader, 1):
            ncols = len(row)
            if ncols != _TBLASTX_EXPECTED_COLUMN_COUNT:
                raise ValueError(
                    f"{tsv_path}:{line_num}: "
                    f"expected {_TBLASTX_EXPECTED_COLUMN_COUNT} columns, "
                    f"got {ncols}",
                )
            hits.append(_parse_row(row, tsv_stem))

    logger.info("Parsed %d BLAST hits from %s (tsv_stem=%s)", len(hits), tsv_path, tsv_stem)
    return hits
