# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Parse NCBI GenBank Flat File (GBFF) to extract CDS features and metadata.

Delegates format parsing to BioPython's SeqIO (GenBank parser), then maps
the result into lightweight domain dataclasses.
"""

from __future__ import annotations

import gzip
import logging
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, TextIO

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

logger = logging.getLogger(__name__)


@dataclass
class CdsFeature:
    """CDS feature extracted from a GBFF record."""

    location_spans: list[tuple[int, int]]
    strand: int | None
    gene: Optional[str] = None
    product: Optional[str] = None
    protein_id: Optional[str] = None
    locus_tag: Optional[str] = None
    codon_start: Optional[int] = None
    transl_table: Optional[int] = None


@dataclass
class GbffRecord:
    """Metadata and CDS features for a single GBFF record."""

    accession: str
    organism: Optional[str] = None
    taxid: Optional[int] = None
    length: Optional[int] = None
    segment: Optional[str] = None
    host: Optional[str] = None
    cds_features: list[CdsFeature] = field(default_factory=list)


def _extract_taxid(feature: SeqFeature) -> int | None:
    for xref in feature.qualifiers.get("db_xref", []):
        if xref.startswith("taxon:"):
            return int(xref.split("taxon:", 1)[1])
    return None


def _first_qualifier(feature: SeqFeature, key: str) -> str | None:
    values = feature.qualifiers.get(key)
    if values:
        return values[0]
    return None


def _int_qualifier(feature: SeqFeature, key: str) -> int | None:
    raw = _first_qualifier(feature, key)
    if raw is not None:
        return int(raw)
    return None


def _has_remote_reference(feature: SeqFeature) -> bool:
    if feature.location is None:
        return False
    return any(getattr(p, "ref", None) is not None for p in feature.location.parts)


def _build_cds_feature(feature: SeqFeature) -> CdsFeature | None:
    if feature.location is None or _has_remote_reference(feature):
        return None
    spans = [(int(part.start) + 1, int(part.end)) for part in feature.location.parts]
    if not spans:
        return None
    return CdsFeature(
        location_spans=spans,
        strand=feature.location.strand,
        gene=_first_qualifier(feature, "gene"),
        product=_first_qualifier(feature, "product"),
        protein_id=_first_qualifier(feature, "protein_id"),
        locus_tag=_first_qualifier(feature, "locus_tag"),
        codon_start=_int_qualifier(feature, "codon_start"),
        transl_table=_int_qualifier(feature, "transl_table"),
    )


def _convert_record(seq_record) -> GbffRecord | None:
    accession = seq_record.id
    if not accession or accession == "<unknown id>":
        return None

    organism: str | None = seq_record.annotations.get("organism")
    length = len(seq_record.seq)
    taxid: int | None = None
    segment: str | None = None
    host: str | None = None
    cds_features: list[CdsFeature] = []

    for feature in seq_record.features:
        if feature.type == "source" and taxid is None:
            taxid = _extract_taxid(feature)
            segment = _first_qualifier(feature, "segment")
            host = _first_qualifier(feature, "host")
        elif feature.type == "CDS":
            cds = _build_cds_feature(feature)
            if cds is not None:
                cds_features.append(cds)

    return GbffRecord(
        accession=accession,
        organism=organism,
        taxid=taxid,
        length=length,
        segment=segment,
        host=host,
        cds_features=cds_features,
    )


def parse_gbff(stream: TextIO) -> Iterator[GbffRecord]:
    """Parse a GBFF stream, yielding one GbffRecord per record."""
    for seq_record in SeqIO.parse(stream, "genbank"):
        record = _convert_record(seq_record)
        if record is not None:
            yield record


def parse_gbff_gzipped(path: Path) -> Iterator[GbffRecord]:
    """Parse a gzip-compressed GBFF file."""
    with gzip.open(path, "rt") as f:
        yield from parse_gbff(f)
