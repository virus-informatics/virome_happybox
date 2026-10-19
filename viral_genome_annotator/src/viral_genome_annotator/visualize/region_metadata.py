# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Extract region metadata (taxid, host, segment) from GFF3 and resolve taxonomy."""

from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from viral_genome_annotator.accession import strip_accession_version
from viral_genome_annotator.visualize.annotation import (
    GFF3_COLUMN_COUNT,
    parse_gff3_attributes,
)

logger = logging.getLogger(__name__)

_TAXON_DBXREF_PATTERN = re.compile(r"taxon:(\d+)")
_UNRECORDED_TAXID_LOG_SAMPLE_COUNT = 3


@dataclass(frozen=True, slots=True)
class RegionMeta:
    sseqid: str
    slen: int
    accession_taxid: int
    host: str | None
    segment: str | None


@dataclass(frozen=True, slots=True)
class ResolvedRegionMeta:
    sseqid: str
    slen: int
    accession_taxid: int
    family: str
    genus: str
    host: str | None
    segment: str | None
    assembly: str | None
    genome_len: int

    @property
    def group_key(self) -> str:
        """Key that groups tracks into one figure.

        The Assembly when NCBI Assembly covers this accession. Otherwise the accession
        itself, forming a single-sequence group.
        """
        return self.assembly if self.assembly is not None else self.sseqid


@dataclass(frozen=True, slots=True)
class GroupTaxonomy:
    """Taxonomy shown for one figure, and the Assembly-level taxid it came from.

    assembly_taxid is None when the group has no Assembly-level taxid and fell back to
    the taxid its own accession declares.
    """

    family: str
    genus: str
    assembly_taxid: int | None


def _load_assembly_mappings(assembly_map_path: Path) -> dict:
    if not assembly_map_path.exists():
        raise FileNotFoundError(f"assembly_map.json not found: {assembly_map_path}")
    return json.loads(assembly_map_path.read_text()).get("mappings", {})


def build_nc_to_assembly_map(
    assembly_map_path: Path,
) -> dict[str, str]:
    """Build {NC_versioned: GCF_stripped} reverse lookup from assembly_map.json.

    Also builds an internal version-stripped index for fallback matching.
    """
    mappings = _load_assembly_mappings(assembly_map_path)

    nc_to_assembly: dict[str, str] = {}
    for gcf_key, entry in mappings.items():
        for seq in entry.get("sequences", []):
            refseq_acc = seq.get("refseq_accession")
            if refseq_acc:
                nc_to_assembly[refseq_acc] = gcf_key

    logger.info(
        "Built NC→Assembly reverse map: %d accessions from %s",
        len(nc_to_assembly),
        assembly_map_path,
    )
    return nc_to_assembly


def build_assembly_to_nc_map(
    assembly_map_path: Path,
) -> dict[str, list[str]]:
    """Build {GCF_stripped: [NC_versioned, ...]} forward lookup."""
    mappings = _load_assembly_mappings(assembly_map_path)

    assembly_to_ncs: dict[str, list[str]] = {}
    for gcf_key, entry in mappings.items():
        ncs: list[str] = []
        for seq in entry.get("sequences", []):
            refseq_acc = seq.get("refseq_accession")
            if refseq_acc:
                ncs.append(refseq_acc)
        if ncs:
            assembly_to_ncs[gcf_key] = ncs

    logger.info(
        "Built Assembly→NC forward map: %d assemblies from %s",
        len(assembly_to_ncs),
        assembly_map_path,
    )
    return assembly_to_ncs


def build_assembly_to_taxid_map(
    assembly_map_path: Path,
) -> dict[str, int]:
    """Build {GCF_stripped: assembly-level taxid} from assembly_map.json.

    This taxid is the assignment NCBI made for the Assembly as a whole, which is
    distinct from the per-accession taxid each GBFF record carries.
    """
    mappings = _load_assembly_mappings(assembly_map_path)

    assembly_to_taxid: dict[str, int] = {}
    unrecorded: list[str] = []
    for gcf_key, entry in mappings.items():
        taxid = entry.get("taxid")
        if taxid is None:
            unrecorded.append(gcf_key)
            continue
        assembly_to_taxid[gcf_key] = taxid

    if unrecorded:
        # An assembly_map.json written before assembly-level taxids were recorded
        logger.warning(
            "No assembly-level taxid for %d assemblies, e.g. %s",
            len(unrecorded),
            ", ".join(unrecorded[:_UNRECORDED_TAXID_LOG_SAMPLE_COUNT]),
        )

    logger.info(
        "Built Assembly→taxid map: %d assemblies from %s",
        len(assembly_to_taxid),
        assembly_map_path,
    )
    return assembly_to_taxid


def _resolve_one_group_taxonomy(
    metas: list[ResolvedRegionMeta],
    assembly_to_taxid: dict[str, int],
    taxonomy: dict[int, tuple[str, str]],
) -> GroupTaxonomy:
    """Resolve the taxonomy of one figure group.

    group_key is derived from assembly, so every member of a group shares its Assembly.
    """
    assembly = metas[0].assembly
    assembly_taxid = assembly_to_taxid.get(assembly) if assembly is not None else None

    if assembly_taxid is not None:
        tax = taxonomy.get(assembly_taxid)
        if tax is None:
            logger.warning(
                "Assembly-level taxid %d of %s not found in taxonomy",
                assembly_taxid,
                assembly,
            )
        family, genus = tax if tax else ("", "")
        return GroupTaxonomy(family=family, genus=genus, assembly_taxid=assembly_taxid)

    if assembly is None:
        # RefSeq release and NCBI Assembly differ in coverage, so this is a normal path
        logger.debug("%s has no Assembly, using the taxonomy of its own taxid", metas[0].sseqid)
    else:
        logger.warning(
            "No Assembly-level taxid for %s, using the taxonomy of its accessions",
            assembly,
        )

    first = min(metas, key=lambda m: m.sseqid)
    return GroupTaxonomy(family=first.family, genus=first.genus, assembly_taxid=None)


def resolve_group_taxonomy(
    region_meta: list[ResolvedRegionMeta],
    assembly_to_taxid: dict[str, int],
    taxonomy: dict[int, tuple[str, str]],
) -> dict[str, GroupTaxonomy]:
    """Resolve the taxonomy shown for each figure group.

    The Assembly-level taxid decides it when NCBI assigned one. Otherwise the group
    falls back to the taxid its own accession declares, and assembly_taxid stays None.
    """
    metas_by_group_key: dict[str, list[ResolvedRegionMeta]] = {}
    for meta in region_meta:
        metas_by_group_key.setdefault(meta.group_key, []).append(meta)

    return {
        group_key: _resolve_one_group_taxonomy(metas, assembly_to_taxid, taxonomy)
        for group_key, metas in metas_by_group_key.items()
    }


def _build_version_stripped_index(
    nc_to_assembly: dict[str, str],
) -> dict[str, tuple[str, str]]:
    """Build {accession_base: (versioned_accession, assembly)} for fallback."""
    index: dict[str, tuple[str, str]] = {}
    for nc_versioned, assembly in nc_to_assembly.items():
        base = strip_accession_version(nc_versioned)
        index[base] = (nc_versioned, assembly)
    return index


def resolve_assembly_for_accession(
    sseqid: str,
    nc_to_assembly: dict[str, str],
    version_stripped_index: dict[str, tuple[str, str]],
) -> str | None:
    """Resolve the Assembly of a single accession, with a version-stripped fallback.

    Returns None when assembly_map.json holds no Assembly for the accession. Callers
    decide what to use in its place; "not found" is not an Assembly.
    """
    exact = nc_to_assembly.get(sseqid)
    if exact is not None:
        return exact

    base = strip_accession_version(sseqid)
    fallback = version_stripped_index.get(base)
    if fallback is not None:
        matched_versioned, assembly = fallback
        logger.warning(
            "%s not in assembly_map, but %s matched %s (possible version mismatch)",
            sseqid,
            base,
            matched_versioned,
        )
        return assembly

    logger.debug("%s has no Assembly in assembly_map", sseqid)
    return None


def _extract_taxid(attrs: dict[str, str]) -> int | None:
    dbxref = attrs.get("Dbxref", "")
    m = _TAXON_DBXREF_PATTERN.search(dbxref)
    return int(m.group(1)) if m else None


def extract_region_metadata(
    gff3_path: Path,
    target_sseqids: set[str],
) -> list[RegionMeta]:
    """Extract region rows (gbkey=Src) from GFF3 for given accessions."""
    if not gff3_path.exists():
        raise FileNotFoundError(f"GFF3 not found: {gff3_path}")

    results: list[RegionMeta] = []

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
            if fields[2] != "region":
                continue

            attrs = parse_gff3_attributes(fields[8])
            if attrs.get("gbkey") != "Src":
                continue

            accession_taxid = _extract_taxid(attrs)
            if accession_taxid is None:
                logger.warning("No taxon Dbxref in region for %s, skipping", seqid)
                continue

            slen = int(fields[4]) - int(fields[3]) + 1

            results.append(
                RegionMeta(
                    sseqid=seqid,
                    slen=slen,
                    accession_taxid=accession_taxid,
                    host=attrs.get("nat-host"),
                    segment=attrs.get("segment"),
                ),
            )

    logger.info(
        "Extracted %d region metadata for %d target accessions from %s",
        len(results),
        len(target_sseqids),
        gff3_path,
    )
    return results


def load_taxonomy_lookup(
    taxonomy_tsv_path: Path,
) -> dict[int, tuple[str, str]]:
    """Load taxonomy.tsv → {taxid: (family, genus)}."""
    if not taxonomy_tsv_path.exists():
        raise FileNotFoundError(f"Taxonomy TSV not found: {taxonomy_tsv_path}")

    lookup: dict[int, tuple[str, str]] = {}
    with taxonomy_tsv_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxid = int(row["taxid"])
            lookup[taxid] = (row["Family"], row["Genus"])

    logger.info("Loaded %d taxonomy entries from %s", len(lookup), taxonomy_tsv_path)
    return lookup


def resolve_region_metadata(
    regions: list[RegionMeta],
    taxonomy: dict[int, tuple[str, str]],
    nc_to_assembly: dict[str, str],
) -> list[ResolvedRegionMeta]:
    """Join region metadata with taxonomy and assembly map."""
    missing_taxids: set[int] = set()
    stripped_index = _build_version_stripped_index(nc_to_assembly)

    # Resolve taxonomy + assembly per accession, accumulate genome_len by group_key
    resolved_fields: list[tuple[RegionMeta, str, str, str | None]] = []
    slen_by_group_key: dict[str, int] = {}
    for r in regions:
        tax = taxonomy.get(r.accession_taxid)
        if tax is None and r.accession_taxid not in missing_taxids:
            missing_taxids.add(r.accession_taxid)
            logger.warning(
                "taxid %d not found in taxonomy for %s",
                r.accession_taxid,
                r.sseqid,
            )

        family, genus = tax if tax else ("", "")
        assembly = resolve_assembly_for_accession(
            r.sseqid,
            nc_to_assembly,
            stripped_index,
        )
        resolved_fields.append((r, family, genus, assembly))

        group_key = assembly if assembly is not None else r.sseqid
        slen_by_group_key[group_key] = slen_by_group_key.get(group_key, 0) + r.slen

    # Construct ResolvedRegionMeta with genome_len resolved
    results: list[ResolvedRegionMeta] = []
    for r, family, genus, assembly in resolved_fields:
        group_key = assembly if assembly is not None else r.sseqid
        results.append(
            ResolvedRegionMeta(
                sseqid=r.sseqid,
                slen=r.slen,
                accession_taxid=r.accession_taxid,
                family=family,
                genus=genus,
                host=r.host,
                segment=r.segment,
                assembly=assembly,
                genome_len=slen_by_group_key[group_key],
            ),
        )

    return results
