# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Download and cache the NCBI taxdump that taxonkit reads."""

from __future__ import annotations

import json
import logging
import tarfile
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from http import HTTPStatus
from io import BytesIO
from pathlib import Path

logger = logging.getLogger(__name__)

# Directory the taxdump is downloaded into, under the run's output directory. It is shared
# by every virus family, so it sits beside the per-family directories rather than inside one.
TAXDUMP_DIR_NAME = "taxdump"

_NCBI_TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
_NCBI_TAXDUMP_ARCHIVE_FILENAME = "taxdump.tar.gz"
# taxonkit reads names.dmp and nodes.dmp; delnodes.dmp and merged.dmp let it resolve
# deleted and merged taxids instead of reporting them as unknown.
_TAXDUMP_REQUIRED_FILENAMES = ("names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp")
_TAXDUMP_DOWNLOAD_META_FILENAME = "download_meta.json"
_TAXDUMP_DOWNLOAD_CHUNK_BYTES = 1024 * 1024
_TAXDUMP_DOWNLOAD_MAX_ATTEMPTS = 3
_TAXDUMP_DOWNLOAD_BACKOFF_BASE_SEC = 2
_TAXDUMP_PROGRESS_LOG_STEP_PCT = 10
_TAXDUMP_HEAD_MAX_ATTEMPTS = 2
_TAXDUMP_HEAD_BACKOFF_BASE_SEC = 5


def _head_last_modified(url: str) -> str | None:
    request = urllib.request.Request(url, method="HEAD")

    def send_head() -> str | None:
        with urllib.request.urlopen(request) as response:  # noqa: S310
            return response.headers.get("Last-Modified")

    try:
        return _call_http_with_retry(
            send_head,
            label=f"{_NCBI_TAXDUMP_ARCHIVE_FILENAME} HEAD",
            max_attempts=_TAXDUMP_HEAD_MAX_ATTEMPTS,
            backoff_base_sec=_TAXDUMP_HEAD_BACKOFF_BASE_SEC,
        )
    except urllib.error.HTTPError as e:
        raise FileNotFoundError(f"Remote file not found: {url} (HTTP {e.code})") from e


def _load_download_meta(meta_path: Path) -> dict:
    if meta_path.exists():
        return json.loads(meta_path.read_text())
    return {}


def _save_download_meta(meta_path: Path, meta: dict) -> None:
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")


def _download_file(url: str, dest: Path) -> int:
    partial_path = dest.with_name(dest.name + ".partial")
    request = urllib.request.Request(url)
    try:
        with urllib.request.urlopen(request) as response:  # noqa: S310
            total_bytes = int(response.headers.get("Content-Length", 0))
            downloaded_bytes = 0
            last_logged_pct = -1
            with partial_path.open("wb") as handle:
                while chunk := response.read(_TAXDUMP_DOWNLOAD_CHUNK_BYTES):
                    handle.write(chunk)
                    downloaded_bytes += len(chunk)
                    if total_bytes <= 0:
                        continue
                    pct = downloaded_bytes * 100 // total_bytes
                    if pct >= last_logged_pct + _TAXDUMP_PROGRESS_LOG_STEP_PCT:
                        last_logged_pct = pct
                        logger.info(
                            "%s: %d / %d MB (%d%%)",
                            dest.name,
                            downloaded_bytes // (1024 * 1024),
                            total_bytes // (1024 * 1024),
                            pct,
                        )
        partial_path.replace(dest)
    except BaseException:
        partial_path.unlink(missing_ok=True)
        raise
    return downloaded_bytes


def _call_http_with_retry[T](
    func: Callable[[], T],
    *,
    label: str,
    max_attempts: int,
    backoff_base_sec: int,
) -> T:
    for attempt in range(max_attempts):
        try:
            return func()
        except (urllib.error.URLError, OSError) as e:
            if isinstance(e, urllib.error.HTTPError) and e.code == HTTPStatus.NOT_FOUND:
                raise
            if attempt == max_attempts - 1:
                raise
            wait_sec = backoff_base_sec ** (attempt + 1)
            logger.warning(
                "%s: attempt %d/%d failed (%s), retrying in %ds",
                label,
                attempt + 1,
                max_attempts,
                e,
                wait_sec,
            )
            time.sleep(wait_sec)
    raise RuntimeError("unreachable")


def _download_with_retry(url: str, dest: Path) -> int:
    return _call_http_with_retry(
        lambda: _download_file(url, dest),
        label=dest.name,
        max_attempts=_TAXDUMP_DOWNLOAD_MAX_ATTEMPTS,
        backoff_base_sec=_TAXDUMP_DOWNLOAD_BACKOFF_BASE_SEC,
    )


def validate_taxdump_dir(taxdump_dir: Path) -> None:
    """Reject a taxdump directory that lacks a file taxonkit reads.

    For a directory the caller supplied: it is never downloaded into, so a missing file has
    to be reported rather than fetched.
    """
    missing = [name for name in _TAXDUMP_REQUIRED_FILENAMES if not (taxdump_dir / name).exists()]
    if missing:
        raise FileNotFoundError(
            f"taxdump directory is missing required files: {missing} (in {taxdump_dir})",
        )


def _extract_taxdump(archive_bytes: bytes, taxdump_dir: Path) -> None:
    with tarfile.open(fileobj=BytesIO(archive_bytes), mode="r:gz") as archive:
        for member in archive.getmembers():
            if member.name in _TAXDUMP_REQUIRED_FILENAMES:
                archive.extract(member, taxdump_dir, filter="data")
    missing = [name for name in _TAXDUMP_REQUIRED_FILENAMES if not (taxdump_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"taxdump archive is missing required files: {missing}")


def has_taxdump(taxdump_dir: Path) -> bool:
    """Check whether all required taxdump files exist in the directory."""
    return all((taxdump_dir / name).exists() for name in _TAXDUMP_REQUIRED_FILENAMES)


def download_taxdump(taxdump_dir: Path) -> Path:
    """Download and extract the NCBI taxdump into taxdump_dir, reusing an up-to-date copy.

    taxdump_dir is the directory taxonkit is pointed at, not its parent, so that the
    caller decides where the extracted files land.
    """
    meta_path = taxdump_dir / _TAXDUMP_DOWNLOAD_META_FILENAME
    meta = _load_download_meta(meta_path)

    logger.info("taxdump: checking Last-Modified...")
    remote_last_modified = _head_last_modified(_NCBI_TAXDUMP_URL)

    cached = meta.get(_NCBI_TAXDUMP_ARCHIVE_FILENAME, {})
    # The metadata alone is not enough: it would stand in for files that are no longer there
    has_all_files = all((taxdump_dir / name).exists() for name in _TAXDUMP_REQUIRED_FILENAMES)
    if (
        has_all_files
        and cached.get("last_modified")
        and cached["last_modified"] == remote_last_modified
    ):
        logger.info("taxdump: cache valid, skipping")
        return taxdump_dir

    logger.info("taxdump: starting download")
    taxdump_dir.mkdir(parents=True, exist_ok=True)
    archive_path = taxdump_dir / _NCBI_TAXDUMP_ARCHIVE_FILENAME
    downloaded_bytes = _download_with_retry(_NCBI_TAXDUMP_URL, archive_path)

    logger.info("taxdump: extracting...")
    _extract_taxdump(archive_path.read_bytes(), taxdump_dir)
    archive_path.unlink()

    meta[_NCBI_TAXDUMP_ARCHIVE_FILENAME] = {
        "last_modified": remote_last_modified,
        "size": downloaded_bytes,
    }
    _save_download_meta(meta_path, meta)
    logger.info("taxdump: done (%d MB)", downloaded_bytes // (1024 * 1024))

    return taxdump_dir
