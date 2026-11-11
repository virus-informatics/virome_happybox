# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Build added-sequence taxonomy metadata from name_map, host CSV, and taxonkit."""

from __future__ import annotations

import csv
import logging
import subprocess
from pathlib import Path

from viral_phylogeny_builder._accession import strip_accession_version
from viral_phylogeny_builder._fasta import iter_fasta
from viral_phylogeny_builder.update_reference.host_classification import (
    HOST_MANUAL_KEY_PATH,
    HostManualKeyRule,
    classify_host_lineage,
    read_host_manual_key,
)

logger = logging.getLogger(__name__)

QUERY_METADATA_COLUMNS = (
    "Accession",
    "Genus",
    "manual_host",
    "Organism_Name",
    "GenBank_RefSeq",
    "Host",
    "Lineage",
)

QUERY_METADATA_FILENAME_SUFFIX = ".query_metadata.csv"

UNKNOWN = "unknown"

_NAME_MAP_NODE_COLUMN = "node"
_NAME_MAP_CONTIG_COLUMN = "contig"
_NAME_MAP_ORF_ID_COLUMN = "orf_id"

_HOST_CSV_CONTIG_CANDIDATES = (
    "contig",
    "query_name",
    "sequence",
    "sequence_id",
    "id",
    "accession",
)
_HOST_CSV_HOST_CANDIDATES = (
    "host",
    "Host",
    "ScientificName",
    "OrganismName",
    "Organism Name",
)


def query_metadata_filename(virus_family: str) -> str:
    return f"{virus_family}{QUERY_METADATA_FILENAME_SUFFIX}"


def _normalize_column_name(name: str) -> str:
    return name.strip().lower().replace("_", "").replace(" ", "").replace("-", "")


def _pick_column(
    fieldnames: list[str],
    candidates: tuple[str, ...],
    label: str,
) -> str:
    norm = {_normalize_column_name(name): name for name in fieldnames}
    for candidate in candidates:
        key = _normalize_column_name(candidate)
        if key in norm:
            return norm[key]
    raise ValueError(
        f"Could not find {label} column in CSV. "
        f"Expected one of: {', '.join(candidates)}. "
        f"Found: {', '.join(fieldnames)}",
    )


def _read_name_map(name_map_tsv: Path) -> dict[str, tuple[str, str]]:
    """Read name_map TSV and return node -> (contig, orf_id) mapping."""
    result: dict[str, tuple[str, str]] = {}
    with name_map_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"name_map TSV has no header: {name_map_tsv}")
        for row in reader:
            node = row.get(_NAME_MAP_NODE_COLUMN, "").strip()
            contig = row.get(_NAME_MAP_CONTIG_COLUMN, "").strip()
            orf_id = row.get(_NAME_MAP_ORF_ID_COLUMN, "").strip()
            if node and contig:
                result[node] = (contig, orf_id)
    return result


def _extract_valid_ids(alignment_fasta: Path) -> list[str]:
    """Read FASTA headers and return unique record IDs (order-preserved)."""
    seen: set[str] = set()
    labels: list[str] = []
    for header, _seq in iter_fasta(alignment_fasta):
        record_id = header.split()[0]
        if record_id and record_id not in seen:
            labels.append(record_id)
            seen.add(record_id)
    return labels


def _read_host_csv(host_csv: Path) -> dict[str, str]:
    """Read host CSV and return normalized_contig -> host mapping."""
    with host_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"host CSV has no header: {host_csv}")

        contig_col = _pick_column(
            list(reader.fieldnames),
            _HOST_CSV_CONTIG_CANDIDATES,
            "contig",
        )
        host_col = _pick_column(
            list(reader.fieldnames),
            _HOST_CSV_HOST_CANDIDATES,
            "host",
        )

        host_map: dict[str, str] = {}
        for row in reader:
            raw_contig = row.get(contig_col, "").strip()
            contig = strip_accession_version(raw_contig)
            host = row.get(host_col, "").strip() or UNKNOWN

            if not contig:
                continue
            if contig in host_map and host_map[contig] != host:
                logger.warning(
                    "Conflicting host for %s: %s vs %s; keeping first",
                    contig,
                    host_map[contig],
                    host,
                )
                continue
            host_map[contig] = host

    if not host_map:
        raise ValueError(f"No contig-host rows read from {host_csv}")

    return host_map


def _run_taxonkit_name2taxid(
    names: list[str],
    data_dir_args: list[str],
    *,
    taxonkit_cmd: str = "taxonkit",
) -> dict[str, str]:
    """Run taxonkit name2taxid and return name -> taxid mapping."""
    result = subprocess.run(
        [taxonkit_cmd, "name2taxid", "--sci-name", *data_dir_args],
        input="\n".join(names) + "\n",
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"taxonkit name2taxid failed (exit {result.returncode}): {result.stderr}",
        )
    name_to_taxid: dict[str, str] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].strip():
            name_to_taxid[parts[0].strip()] = parts[1].strip()
    return name_to_taxid


def _run_taxonkit_lineage(
    taxids: list[str],
    data_dir_args: list[str],
    *,
    taxonkit_cmd: str = "taxonkit",
) -> dict[str, str]:
    """Run taxonkit lineage and return taxid -> lineage mapping."""
    result = subprocess.run(
        [taxonkit_cmd, "lineage", *data_dir_args],
        input="\n".join(taxids) + "\n",
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"taxonkit lineage failed (exit {result.returncode}): {result.stderr}",
        )
    taxid_to_lineage: dict[str, str] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            taxid_to_lineage[parts[0].strip()] = parts[1].strip()
    return taxid_to_lineage


def _resolve_hosts(
    hosts: set[str],
    rules: tuple[HostManualKeyRule, ...],
    taxdump_dir: Path,
    *,
    taxonkit_cmd: str = "taxonkit",
) -> dict[str, tuple[str, str]]:
    """Resolve host names to (lineage, manual_host) via taxonkit + classify.

    Mirrors archive resolve_hosts: returns host -> (lineage, manual_host).
    """
    result: dict[str, tuple[str, str]] = {}
    unique_hosts = sorted(hosts)

    resolvable = [h for h in unique_hosts if h != UNKNOWN]
    if not resolvable:
        for h in unique_hosts:
            result[h] = (UNKNOWN, UNKNOWN)
        return result

    data_dir_args = ["--data-dir", str(taxdump_dir)]

    name_to_taxid = _run_taxonkit_name2taxid(resolvable, data_dir_args, taxonkit_cmd=taxonkit_cmd)

    unresolved = set(resolvable) - set(name_to_taxid)
    if unresolved:
        logger.warning(
            "%d host names not resolved by taxonkit name2taxid (e.g. %s)",
            len(unresolved),
            ", ".join(sorted(unresolved)[:5]),
        )

    taxid_to_lineage: dict[str, str] = {}
    if name_to_taxid:
        taxid_to_lineage = _run_taxonkit_lineage(
            sorted(set(name_to_taxid.values())),
            data_dir_args,
            taxonkit_cmd=taxonkit_cmd,
        )

    for host in unique_hosts:
        if host == UNKNOWN:
            result[host] = (UNKNOWN, UNKNOWN)
            continue

        taxid = name_to_taxid.get(host, "")
        lineage = taxid_to_lineage.get(taxid, "") if taxid else ""

        if not lineage:
            logger.warning("Could not resolve host taxonomy with taxonkit: %s", host)
            result[host] = (UNKNOWN, UNKNOWN)
        else:
            classification = classify_host_lineage(lineage, rules)
            result[host] = (lineage, classification.manual_host_class)

    return result


def prepare_query_metadata(
    name_map_tsv: Path,
    alignment_fasta: Path,
    host_csv: Path,
    virus_family: str,
    output_csv: Path,
    *,
    taxdump_dir: Path,
    taxonkit_cmd: str = "taxonkit",
) -> None:
    """Build added-sequence taxonomy metadata CSV.

    Produces one row per added sequence (those present in both alignment_fasta
    and name_map_tsv). Reference sequences in the alignment are ignored.
    """
    valid_ids = _extract_valid_ids(alignment_fasta)
    node_to_contig_orf = _read_name_map(name_map_tsv)
    host_map = _read_host_csv(host_csv)
    rules = read_host_manual_key(HOST_MANUAL_KEY_PATH)

    # ── 1. Build rows with placeholder manual_host (archive L485-507) ──
    rows: list[dict[str, str]] = []
    missing_host_count = 0

    for valid_id in valid_ids:
        if valid_id not in node_to_contig_orf:
            continue
        contig, orf_id = node_to_contig_orf[valid_id]
        host = host_map.get(contig) or host_map.get(strip_accession_version(contig))
        if host is None:
            missing_host_count += 1
            host = UNKNOWN

        rows.append(
            {
                "Accession": valid_id,
                "Genus": virus_family,
                "manual_host": UNKNOWN,
                "Organism_Name": orf_id,
                "GenBank_RefSeq": "Input",
                "Host": host,
                "Lineage": UNKNOWN,
            },
        )

    if not rows:
        raise ValueError(
            "No added sequences found in both alignment and name_map. "
            "Check that name_map node IDs match alignment FASTA headers.",
        )

    if missing_host_count:
        logger.warning(
            "%d/%d added sequences had no host in host_csv",
            missing_host_count,
            len(rows),
        )

    # ── 2. Resolve host → lineage → manual_host, then overwrite (archive L508-516) ──
    unique_hosts = {row["Host"] for row in rows}
    host_taxonomy = _resolve_hosts(unique_hosts, rules, taxdump_dir, taxonkit_cmd=taxonkit_cmd)

    for row in rows:
        lineage, manual_host = host_taxonomy.get(row["Host"], (UNKNOWN, UNKNOWN))
        row["Lineage"] = lineage or UNKNOWN
        row["manual_host"] = manual_host or UNKNOWN

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(QUERY_METADATA_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)

    logger.info(
        "Query metadata: %d added sequences (%d with host match) -> %s",
        len(rows),
        len(rows) - missing_host_count,
        output_csv,
    )
