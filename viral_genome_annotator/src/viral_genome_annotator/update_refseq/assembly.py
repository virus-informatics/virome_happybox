# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Assembly summary parser and assembly map builder for GCF → NC mapping."""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path

from viral_genome_annotator.accession import strip_accession_version
from viral_genome_annotator.update_refseq.datasets_api import fetch_and_verify_batch_zip
from viral_genome_annotator.update_refseq.download import (
    DOWNLOAD_META_FILENAME,
    _download_with_retry,
    _head_last_modified,
    _load_meta,
    _save_meta,
)

logger = logging.getLogger(__name__)

_NCBI_VIRAL_ASSEMBLY_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt"
)

ASSEMBLY_SUMMARY_FILENAME = "assembly_summary.txt"
ASSEMBLY_MAP_FILENAME = "assembly_map.json"
ASSEMBLY_SUBDIR = "assembly"
_BATCH_MANIFEST_FILENAME = "batch_manifest.jsonl"

_DEFAULT_ASSEMBLY_BATCH_SIZE = 500
# Rate limit is 5 rps (200ms); 3s is intentionally conservative
_ASSEMBLY_BATCH_SLEEP_SEC = 3

# Column layout of assembly_summary.txt, described in its header comment and in
# ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
_ASSEMBLY_SUMMARY_ACCESSION_COLUMN = 0
_ASSEMBLY_SUMMARY_TAXID_COLUMN = 5


@dataclass(frozen=True, slots=True)
class AssemblySummaryRow:
    """One data row of assembly_summary.txt.

    taxid is the assignment NCBI made for the assembly as a whole, which is distinct
    from the per-accession taxid carried by each GBFF record.
    """

    accession: str
    taxid: int | None


def _parse_summary_taxid(columns: list[str], accession: str) -> int | None:
    if len(columns) <= _ASSEMBLY_SUMMARY_TAXID_COLUMN:
        logger.warning("%s: %s has no taxid column", ASSEMBLY_SUMMARY_FILENAME, accession)
        return None
    raw_taxid = columns[_ASSEMBLY_SUMMARY_TAXID_COLUMN].strip()
    try:
        return int(raw_taxid)
    except ValueError:
        logger.warning(
            "%s: %s has non-numeric taxid %r",
            ASSEMBLY_SUMMARY_FILENAME,
            accession,
            raw_taxid,
        )
        return None


def parse_assembly_summary(path: Path) -> list[AssemblySummaryRow]:
    """Parse assembly_summary.txt, return GCF_ accessions with their assembly-level taxid."""
    rows: list[AssemblySummaryRow] = []
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            accession = columns[_ASSEMBLY_SUMMARY_ACCESSION_COLUMN].strip()
            if not accession:
                continue
            rows.append(
                AssemblySummaryRow(
                    accession=accession,
                    taxid=_parse_summary_taxid(columns, accession),
                ),
            )
    return rows


def build_assembly_taxid_map(rows: list[AssemblySummaryRow]) -> dict[str, int]:
    """Build {GCF_stripped: assembly-level taxid} from parsed summary rows."""
    return {
        strip_accession_version(row.accession): row.taxid for row in rows if row.taxid is not None
    }


def determine_fetch_targets(
    all_accessions: list[str],
    existing_mappings: dict[str, dict],
) -> list[str]:
    """Return accessions needing fetch: new or version-updated."""
    targets: list[str] = []
    for acc in all_accessions:
        key = strip_accession_version(acc)
        existing = existing_mappings.get(key)
        if existing is None or existing.get("accession") != acc:
            targets.append(acc)
    return targets


def _load_assembly_map(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {
        "source": ASSEMBLY_SUMMARY_FILENAME,
        "source_record_count": 0,
        "mappings": {},
    }


def _save_assembly_map(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def _ensure_assembly_summary(downloads_dir: Path) -> Path:
    """Fetch assembly_summary.txt if missing or stale."""
    dest = downloads_dir / ASSEMBLY_SUMMARY_FILENAME
    meta_path = downloads_dir / DOWNLOAD_META_FILENAME
    meta = _load_meta(meta_path)

    logger.info("%s: checking Last-Modified...", ASSEMBLY_SUMMARY_FILENAME)
    remote_modified = _head_last_modified(_NCBI_VIRAL_ASSEMBLY_SUMMARY_URL)

    cached = meta.get(ASSEMBLY_SUMMARY_FILENAME, {})
    if dest.exists() and cached.get("last_modified") and cached["last_modified"] == remote_modified:
        logger.info("%s: cache valid, skipping", ASSEMBLY_SUMMARY_FILENAME)
        return dest

    logger.info("%s: downloading...", ASSEMBLY_SUMMARY_FILENAME)
    size = _download_with_retry(_NCBI_VIRAL_ASSEMBLY_SUMMARY_URL, dest)

    meta[ASSEMBLY_SUMMARY_FILENAME] = {
        "url": _NCBI_VIRAL_ASSEMBLY_SUMMARY_URL,
        "last_modified": remote_modified,
        "size": size,
    }
    _save_meta(meta_path, meta)
    logger.info("%s: done (%d KB)", ASSEMBLY_SUMMARY_FILENAME, size // 1024)

    return dest


def download_assembly_map(
    output_dir: Path,
    batch_size: int = _DEFAULT_ASSEMBLY_BATCH_SIZE,
) -> Path:
    """Fetch NCBI Datasets SEQUENCE_REPORT and build assembly map."""
    downloads_dir = output_dir / "downloads"
    assembly_dir = output_dir / ASSEMBLY_SUBDIR
    zips_dir = downloads_dir / "assembly_zips"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    assembly_dir.mkdir(parents=True, exist_ok=True)
    zips_dir.mkdir(parents=True, exist_ok=True)

    summary_path = _ensure_assembly_summary(downloads_dir)

    summary_rows = parse_assembly_summary(summary_path)
    assembly_taxid_map = build_assembly_taxid_map(summary_rows)
    logger.info("Total assemblies in summary: %d", len(summary_rows))

    map_path = assembly_dir / ASSEMBLY_MAP_FILENAME
    assembly_map = _load_assembly_map(map_path)
    targets = determine_fetch_targets(
        [row.accession for row in summary_rows],
        assembly_map.get("mappings", {}),
    )
    logger.info(
        "Fetch targets: %d new/updated, %d already up-to-date",
        len(targets),
        len(summary_rows) - len(targets),
    )

    if not targets:
        logger.info("All assemblies up-to-date, nothing to fetch")
        return map_path

    manifest_path = zips_dir / _BATCH_MANIFEST_FILENAME
    total_batches = (len(targets) + batch_size - 1) // batch_size

    for batch_idx in range(total_batches):
        start = batch_idx * batch_size
        batch = targets[start : start + batch_size]

        if batch_idx > 0:
            time.sleep(_ASSEMBLY_BATCH_SLEEP_SEC)

        logger.info(
            "Batch %d/%d: fetching %d accessions...",
            batch_idx + 1,
            total_batches,
            len(batch),
        )

        fetched = fetch_and_verify_batch_zip(
            batch,
            zips_dir,
            assembly_taxid_map,
            ASSEMBLY_SUMMARY_FILENAME,
        )

        # The manifest records only batches whose zip was extracted successfully
        manifest_entry = {
            "batch": batch_idx,
            "zip": fetched.zip_path.name,
            "fetched_at": fetched.fetched_at,
            "accessions": batch,
        }
        with manifest_path.open("a") as f:
            f.write(json.dumps(manifest_entry, ensure_ascii=False) + "\n")

        assembly_map.setdefault("mappings", {}).update(fetched.mappings)
        assembly_map["source"] = ASSEMBLY_SUMMARY_FILENAME
        assembly_map["source_record_count"] = len(summary_rows)
        _save_assembly_map(map_path, assembly_map)

        logger.info(
            "Batch %d/%d: extracted %d mappings",
            batch_idx + 1,
            total_batches,
            len(fetched.mappings),
        )

    logger.info(
        "Assembly map complete: %d total mappings in %s",
        len(assembly_map.get("mappings", {})),
        map_path,
    )
    return map_path
