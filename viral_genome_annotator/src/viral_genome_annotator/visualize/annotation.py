# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Extract CDS annotations from a concatenated GFF3 file."""

from __future__ import annotations

import logging
import re
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

logger = logging.getLogger(__name__)

_GFF3_COL9_ATTR_PATTERN = re.compile(r"(?:^|;)([^=;]+)=([^;]*)")
GFF3_COLUMN_COUNT = 9
_CDS_FALLBACK_PRODUCT = "CDS"


def parse_gff3_attributes(raw: str) -> dict[str, str]:
    """Parse GFF3 column-9 key=value pairs with RFC 3986 percent-decoding."""
    return {unquote(m.group(1)): unquote(m.group(2)) for m in _GFF3_COL9_ATTR_PATTERN.finditer(raw)}


@dataclass(frozen=True, slots=True)
class CdsAnnotation:
    sseqid: str
    parts: tuple[tuple[int, int], ...]
    gene: str | None
    product: str
    cds_id: str

    @property
    def min_start(self) -> int:
        return min(s for s, _ in self.parts)


def _resolve_cds_id(
    attrs: dict[str, str],
    seqid: str,
    start: int,
    end: int,
) -> str:
    """Resolve a CDS identifier: ID → locus_tag → positional fallback."""
    gff3_id = attrs.get("ID")
    if gff3_id:
        return gff3_id
    locus_tag = attrs.get("locus_tag")
    if locus_tag:
        return locus_tag
    return f"{seqid}_{start}_{end}"


def extract_cds_annotations(
    gff3_path: Path,
    target_sseqids: set[str],
) -> list[CdsAnnotation]:
    """Extract CDS rows from a concatenated GFF3 for the given accessions.

    Multiple GFF3 lines sharing the same ID are merged into a single
    CdsAnnotation with multiple parts (e.g. origin-crossing CDS in
    circular genomes).
    """
    if not gff3_path.exists():
        raise FileNotFoundError(f"GFF3 not found: {gff3_path}")

    @dataclass
    class _CdsPending:
        sseqid: str
        parts: list[tuple[int, int]]
        gene: str | None
        product: str
        cds_id: str

    pending: OrderedDict[str, _CdsPending] = OrderedDict()

    with gff3_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != GFF3_COLUMN_COUNT:
                continue

            seqid = fields[0]
            if seqid not in target_sseqids:
                continue
            if fields[2] != "CDS":
                continue

            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_gff3_attributes(fields[8])
            cds_id = _resolve_cds_id(attrs, seqid, start, end)

            if cds_id in pending:
                pending[cds_id].parts.append((start, end))
            else:
                gene = attrs.get("gene")
                product = attrs.get("product") or gene or _CDS_FALLBACK_PRODUCT
                pending[cds_id] = _CdsPending(
                    sseqid=seqid,
                    parts=[(start, end)],
                    gene=gene,
                    product=product,
                    cds_id=cds_id,
                )

    results = [
        CdsAnnotation(
            sseqid=p.sseqid,
            parts=tuple(p.parts),
            gene=p.gene,
            product=p.product,
            cds_id=p.cds_id,
        )
        for p in pending.values()
    ]

    logger.info(
        "Extracted %d CDS annotations for %d accessions from %s",
        len(results),
        len(target_sseqids),
        gff3_path,
    )
    return results
