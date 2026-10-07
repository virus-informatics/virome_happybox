# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Assembly summary parser and NCBI Datasets API client for GCF → NC mapping."""

from __future__ import annotations

import json
import logging
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

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
_NCBI_DATASETS_GENOME_DOWNLOAD_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"

ASSEMBLY_SUMMARY_FILENAME = "assembly_summary.txt"
ASSEMBLY_MAP_FILENAME = "assembly_map.json"
ASSEMBLY_SUBDIR = "assembly"
_BATCH_MANIFEST_FILENAME = "batch_manifest.jsonl"

_DEFAULT_ASSEMBLY_BATCH_SIZE = 500
# Rate limit is 5 rps (200ms); 3s is intentionally conservative
_ASSEMBLY_BATCH_SLEEP_SEC = 3
_ASSEMBLY_API_MAX_ATTEMPTS = 3
_ASSEMBLY_API_BACKOFF_BASE_SEC = 5

_SEQUENCE_REPORT_FIELD_MAP: dict[str, str] = {
    "refseqAccession": "refseq_accession",
    "chrName": "chr_name",
    "length": "length",
    "genbankAccession": "genbank_accession",
}

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


def strip_gcf_version(accession: str) -> str:
    """GCF_000839185.1 -> GCF_000839185"""
    dot_idx = accession.rfind(".")
    return accession[:dot_idx] if dot_idx != -1 else accession


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
    return {strip_gcf_version(row.accession): row.taxid for row in rows if row.taxid is not None}


def determine_fetch_targets(
    all_accessions: list[str],
    existing_mappings: dict[str, dict],
) -> list[str]:
    """Return accessions needing fetch: new or version-updated."""
    targets: list[str] = []
    for acc in all_accessions:
        key = strip_gcf_version(acc)
        existing = existing_mappings.get(key)
        if existing is None or existing.get("accession") != acc:
            targets.append(acc)
    return targets


def _fetch_batch_zip(accessions: list[str], dest: Path) -> None:
    body = json.dumps(
        {
            "accessions": accessions,
            "include_annotation_type": ["SEQUENCE_REPORT"],
        },
    ).encode()

    req = urllib.request.Request(
        _NCBI_DATASETS_GENOME_DOWNLOAD_URL,
        data=body,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/zip",
        },
        method="POST",
    )

    tmp = dest.with_suffix(".zip.tmp")
    try:
        with urllib.request.urlopen(req) as resp:  # noqa: S310
            with tmp.open("wb") as f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
        tmp.rename(dest)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _fetch_batch_zip_with_retry(accessions: list[str], dest: Path) -> None:
    for attempt in range(_ASSEMBLY_API_MAX_ATTEMPTS):
        try:
            _fetch_batch_zip(accessions, dest)
            return
        except (urllib.error.URLError, OSError) as e:
            if attempt == _ASSEMBLY_API_MAX_ATTEMPTS - 1:
                raise
            wait = _ASSEMBLY_API_BACKOFF_BASE_SEC * (attempt + 1)
            logger.warning(
                "API batch request: attempt %d/%d failed (%s), retrying in %ds",
                attempt + 1,
                _ASSEMBLY_API_MAX_ATTEMPTS,
                e,
                wait,
            )
            time.sleep(wait)


def extract_mappings_from_zip(
    zip_path: Path,
    fetched_at: str,
    assembly_taxid_map: dict[str, int],
) -> dict[str, dict]:
    """Read sequence_report.jsonl from ZIP, return {stripped_gcf: mapping}."""
    mappings: dict[str, dict] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if not name.endswith("/sequence_report.jsonl"):
                continue
            parts = name.split("/")
            if len(parts) < 4:
                continue
            gcf_versioned = parts[-2]

            data = zf.read(name).decode()
            sequences: list[dict] = []
            for line in data.strip().split("\n"):
                if not line:
                    continue
                rec = json.loads(line)
                seq = {
                    snake: rec[camel]
                    for camel, snake in _SEQUENCE_REPORT_FIELD_MAP.items()
                    if camel in rec
                }
                sequences.append(seq)

            key = strip_gcf_version(gcf_versioned)
            taxid = assembly_taxid_map.get(key)
            if taxid is None:
                logger.warning("%s: no assembly-level taxid in %s", key, ASSEMBLY_SUMMARY_FILENAME)
            mappings[key] = {
                "accession": gcf_versioned,
                "fetched_at": fetched_at,
                "taxid": taxid,
                "sequences": sequences,
            }
    return mappings


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

        now = datetime.now(timezone.utc)
        timestamp_compact = now.strftime("%Y%m%dT%H%M%SZ")
        fetched_at = now.isoformat()
        zip_name = f"batch_{timestamp_compact}.zip"
        zip_path = zips_dir / zip_name

        logger.info(
            "Batch %d/%d: fetching %d accessions...",
            batch_idx + 1,
            total_batches,
            len(batch),
        )

        _fetch_batch_zip_with_retry(batch, zip_path)

        manifest_entry = {
            "batch": batch_idx,
            "zip": zip_name,
            "fetched_at": fetched_at,
            "accessions": batch,
        }
        with manifest_path.open("a") as f:
            f.write(json.dumps(manifest_entry, ensure_ascii=False) + "\n")

        batch_mappings = extract_mappings_from_zip(zip_path, fetched_at, assembly_taxid_map)
        assembly_map.setdefault("mappings", {}).update(batch_mappings)
        assembly_map["source"] = ASSEMBLY_SUMMARY_FILENAME
        assembly_map["source_record_count"] = len(summary_rows)
        _save_assembly_map(map_path, assembly_map)

        logger.info(
            "Batch %d/%d: extracted %d mappings",
            batch_idx + 1,
            total_batches,
            len(batch_mappings),
        )

    logger.info(
        "Assembly map complete: %d total mappings in %s",
        len(assembly_map.get("mappings", {})),
        map_path,
    )
    return map_path
