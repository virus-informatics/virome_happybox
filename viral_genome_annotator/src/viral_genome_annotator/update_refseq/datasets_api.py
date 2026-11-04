# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""NCBI Datasets API client: fetch SEQUENCE_REPORT zips and extract GCF → NC mappings."""

from __future__ import annotations

import http.client
import json
import logging
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from dataclasses import dataclass
from datetime import datetime, timezone
from http import HTTPStatus
from pathlib import Path

from viral_genome_annotator.accession import strip_accession_version

logger = logging.getLogger(__name__)

_NCBI_DATASETS_GENOME_DOWNLOAD_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"

_ASSEMBLY_API_MAX_ATTEMPTS = 5
_ASSEMBLY_API_BACKOFF_BASE_SEC = 5

_SEQUENCE_REPORT_FIELD_MAP: dict[str, str] = {
    "refseqAccession": "refseq_accession",
    "chrName": "chr_name",
    "length": "length",
    "genbankAccession": "genbank_accession",
}


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


@dataclass(frozen=True, slots=True)
class VerifiedBatchZip:
    """One batch zip verified by extracting its mappings, named after the attempt that succeeded."""

    zip_path: Path
    fetched_at: str
    mappings: dict[str, dict]


def fetch_and_verify_batch_zip(
    accessions: list[str],
    zips_dir: Path,
    assembly_taxid_map: dict[str, int],
    taxid_source_filename: str,
) -> VerifiedBatchZip:
    """Fetch one batch zip and verify it by extracting its mappings.

    The zip name and fetched_at come from the attempt that succeeded.
    """
    for attempt in range(_ASSEMBLY_API_MAX_ATTEMPTS):
        fetched_at = datetime.now(timezone.utc)
        zip_path = zips_dir / f"batch_{fetched_at.strftime('%Y%m%dT%H%M%SZ')}.zip"
        try:
            _fetch_batch_zip(accessions, zip_path)
            mappings = extract_mappings_from_zip(
                zip_path,
                fetched_at.isoformat(),
                assembly_taxid_map,
                taxid_source_filename,
            )
            return VerifiedBatchZip(
                zip_path=zip_path,
                fetched_at=fetched_at.isoformat(),
                mappings=mappings,
            )
        # 0-byte and truncated bodies surface as BadZipFile on extraction; a corrupted deflate
        # stream surfaces as either BadZipFile or zlib.error
        except (
            urllib.error.URLError,
            OSError,
            http.client.HTTPException,
            zipfile.BadZipFile,
            zlib.error,
        ) as e:
            if isinstance(e, urllib.error.HTTPError) and e.code == HTTPStatus.TOO_MANY_REQUESTS:
                logger.error(
                    "API batch request: rate limited (HTTP 429), aborting; batches fetched so far"
                    " are saved, rerun later to resume",
                )
                raise
            # Resending the same request cannot change a client error
            if isinstance(e, urllib.error.HTTPError) and e.code < HTTPStatus.INTERNAL_SERVER_ERROR:
                raise
            if isinstance(e, (zipfile.BadZipFile, zlib.error)):
                # Unlike a failed download, the zip was already renamed into place
                zip_path.unlink(missing_ok=True)
            if attempt == _ASSEMBLY_API_MAX_ATTEMPTS - 1:
                raise
            wait = _ASSEMBLY_API_BACKOFF_BASE_SEC * 2**attempt
            logger.warning(
                "API batch request: attempt %d/%d failed (%s), retrying in %ds",
                attempt + 1,
                _ASSEMBLY_API_MAX_ATTEMPTS,
                e,
                wait,
            )
            time.sleep(wait)
    raise RuntimeError("unreachable")


def extract_mappings_from_zip(
    zip_path: Path,
    fetched_at: str,
    assembly_taxid_map: dict[str, int],
    taxid_source_filename: str,
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

            key = strip_accession_version(gcf_versioned)
            taxid = assembly_taxid_map.get(key)
            if taxid is None:
                logger.warning("%s: no assembly-level taxid in %s", key, taxid_source_filename)
            mappings[key] = {
                "accession": gcf_versioned,
                "fetched_at": fetched_at,
                "taxid": taxid,
                "sequences": sequences,
            }
    return mappings
