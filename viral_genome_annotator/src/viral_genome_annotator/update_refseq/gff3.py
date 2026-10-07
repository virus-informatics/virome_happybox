# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Write GbffRecord features to GFF3 format.

Outputs a region feature and ##sequence-region directive per record,
followed by CDS lines.
Multi-span CDS (join locations) produce one line per span.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TextIO

from viral_genome_annotator.update_refseq.gbff import CdsFeature, GbffRecord

_GFF3_VERSION_HEADER = "##gff-version 3"
_GFF3_CDS_SOURCE = "RefSeq"

# GFF3 spec "Encoding" section: tab/newline/CR/% must be RFC 3986 percent-encoded.
# "Column 9 reserved characters": ; = & must be escaped in attribute values.
# % must be first to avoid double-encoding.
_GFF3_VALUE_ESCAPE_PAIRS = [
    ("%", "%25"),
    ("\t", "%09"),
    ("\n", "%0A"),
    ("\r", "%0D"),
    (";", "%3B"),
    ("=", "%3D"),
    ("&", "%26"),
]


def _escape_gff3_value(raw: str) -> str:
    for char, encoded in _GFF3_VALUE_ESCAPE_PAIRS:
        raw = raw.replace(char, encoded)
    return raw


# GFF3 spec Column 7: +/- for stranded, . for unstranded, ? for unknown.
# BioPython strand: 1=forward, -1=reverse, 0=unstranded, None=unknown.
_STRAND_TO_GFF3_SYMBOL: dict[int | None, str] = {1: "+", -1: "-", 0: ".", None: "?"}


def _strand_to_gff3(strand: int | None) -> str:
    try:
        return _STRAND_TO_GFF3_SYMBOL[strand]
    except KeyError:
        raise ValueError(f"unexpected strand value: {strand!r}") from None


def _compute_cds_phases(cds: CdsFeature) -> list[str]:
    """Compute GFF3 phase for each location span.

    For multi-span (join) CDS the phase of each subsequent span depends on
    cumulative coding length.  Reading order follows genomic order for
    +strand and reversed genomic order for -strand.
    """
    if cds.codon_start is None:
        return ["." for _ in cds.location_spans]

    first_phase = cds.codon_start - 1
    n = len(cds.location_spans)
    # TODO(VII-23): -strand join phase is spec-based only; no real-data verification yet.
    reading_indices = list(reversed(range(n))) if cds.strand == -1 else list(range(n))

    phases: list[str] = [""] * n
    accumulated_span_len = 0

    for reading_rank, genomic_idx in enumerate(reading_indices):
        start, end = cds.location_spans[genomic_idx]
        span_len = end - start + 1
        if reading_rank == 0:
            phases[genomic_idx] = str(first_phase)
        else:
            coding_bases = accumulated_span_len - first_phase
            phase = (3 - (coding_bases % 3)) % 3
            phases[genomic_idx] = str(phase)
        accumulated_span_len += span_len

    return phases


def _build_cds_attributes(cds: CdsFeature) -> str:
    pairs: list[str] = []
    if cds.protein_id is not None:
        pairs.append(f"ID=cds-{cds.protein_id}")
        pairs.append(f"Name={cds.protein_id}")
    pairs.append("gbkey=CDS")
    if cds.gene is not None:
        pairs.append(f"gene={_escape_gff3_value(cds.gene)}")
    if cds.locus_tag is not None:
        pairs.append(f"locus_tag={_escape_gff3_value(cds.locus_tag)}")
    if cds.product is not None:
        pairs.append(f"product={_escape_gff3_value(cds.product)}")
    if cds.protein_id is not None:
        pairs.append(f"protein_id={cds.protein_id}")
    if cds.transl_table is not None:
        pairs.append(f"transl_table={cds.transl_table}")
    return ";".join(pairs)


def _format_cds_lines(accession: str, cds: CdsFeature) -> list[str]:
    strand = _strand_to_gff3(cds.strand)
    phases = _compute_cds_phases(cds)
    attributes = _build_cds_attributes(cds)

    lines: list[str] = []
    for (start, end), phase in zip(cds.location_spans, phases):
        fields = (
            accession,
            _GFF3_CDS_SOURCE,
            "CDS",
            str(start),
            str(end),
            ".",
            strand,
            phase,
            attributes,
        )
        lines.append("\t".join(fields))
    return lines


def _build_region_attributes(record: GbffRecord) -> str:
    # ID is generated on-the-fly following NCBI annotwriter convention.
    # Per NCBI GFF3 docs, IDs are not stable identifiers:
    # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/
    pairs: list[str] = [f"ID={record.accession}:1..{record.length}"]
    if record.taxid is not None:
        pairs.append(f"Dbxref=taxon:{record.taxid}")
    # "gbkey=Src" identifies the region feature as originating from
    # the GenBank source feature, per NCBI GFF3 format specification:
    # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/#the-source-region-feature
    pairs.append("gbkey=Src")
    if record.host is not None:
        # NCBI annotwriter renders the INSDC /host qualifier as "nat-host"
        # in GFF3 output. No official spec documents this renaming, but
        # BioSample attributes list "nat host" as a synonym for "host",
        # and all NCBI-produced GFF3 files consistently use "nat-host".
        pairs.append(f"nat-host={_escape_gff3_value(record.host)}")
    if record.segment is not None:
        pairs.append(f"segment={_escape_gff3_value(record.segment)}")
    return ";".join(pairs)


def _format_region_line(record: GbffRecord) -> str:
    attrs = _build_region_attributes(record)
    fields = (
        record.accession,
        _GFF3_CDS_SOURCE,
        "region",
        "1",
        str(record.length),
        ".",
        "+",
        ".",
        attrs,
    )
    return "\t".join(fields)


def write_gff3(records: Iterable[GbffRecord], dest: TextIO) -> int:
    """Write region + CDS features in GFF3 format.

    Returns the number of CDS lines written (excluding region and directives).
    """
    dest.write(_GFF3_VERSION_HEADER + "\n")
    feature_line_count = 0
    for record in records:
        if record.length is not None:
            dest.write(f"##sequence-region {record.accession} 1 {record.length}\n")
            dest.write(_format_region_line(record) + "\n")
        for cds in record.cds_features:
            for line in _format_cds_lines(record.accession, cds):
                dest.write(line + "\n")
                feature_line_count += 1
    return feature_line_count


def write_gff3_file(records: Iterable[GbffRecord], path: Path) -> int:
    """Write features to a GFF3 file. Returns the number of feature lines written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        return write_gff3(records, f)
