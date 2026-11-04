# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

from __future__ import annotations

import json
import logging
import tarfile
import time
import urllib.error
import urllib.request
from io import BytesIO
from pathlib import Path

logger = logging.getLogger(__name__)

_NCBI_REFSEQ_VIRAL_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral"
_NCBI_TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
_TAXDUMP_REQUIRED_FILES = ("names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp")

NCBI_DOWNLOAD_FILES: dict[str, str] = {
    "viral.1.1.genomic.fna.gz": f"{_NCBI_REFSEQ_VIRAL_URL}/viral.1.1.genomic.fna.gz",
    "viral.1.genomic.gbff.gz": f"{_NCBI_REFSEQ_VIRAL_URL}/viral.1.genomic.gbff.gz",
}

DOWNLOAD_META_FILENAME = "download_meta.json"
_DOWNLOAD_CHUNK_BYTES = 1024 * 1024
_DOWNLOAD_MAX_ATTEMPTS = 3
_DOWNLOAD_BACKOFF_BASE_SEC = 2  # wait = 2^attempt 秒 (2s, 4s)
_NCBI_REQUEST_INTERVAL_SEC = 1
_PROGRESS_LOG_STEP_PCT = 10


def _head_last_modified(url: str) -> str | None:
    req = urllib.request.Request(url, method="HEAD")
    for attempt in range(_DOWNLOAD_MAX_ATTEMPTS):
        try:
            with urllib.request.urlopen(req) as resp:  # noqa: S310
                return resp.headers.get("Last-Modified")
        except urllib.error.HTTPError as e:
            raise FileNotFoundError(f"Remote file not found: {url} (HTTP {e.code})") from e
        # urllib wraps socket errors (e.g. socket.gaierror) in URLError
        except (urllib.error.URLError, OSError) as e:
            if attempt == _DOWNLOAD_MAX_ATTEMPTS - 1:
                raise
            wait = _DOWNLOAD_BACKOFF_BASE_SEC ** (attempt + 1)
            logger.warning(
                "%s: HEAD attempt %d/%d failed (%s), retrying in %ds",
                url,
                attempt + 1,
                _DOWNLOAD_MAX_ATTEMPTS,
                e,
                wait,
            )
            time.sleep(wait)
    raise RuntimeError("unreachable")


def _load_meta(meta_path: Path) -> dict:
    if meta_path.exists():
        return json.loads(meta_path.read_text())
    return {}


def _save_meta(meta_path: Path, meta: dict) -> None:
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")


def _download_file(url: str, dest: Path) -> int:
    tmp_dest = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url)
    try:
        with urllib.request.urlopen(req) as resp:  # noqa: S310
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            last_logged_pct = -1
            with tmp_dest.open("wb") as f:
                while True:
                    chunk = resp.read(_DOWNLOAD_CHUNK_BYTES)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0:
                        pct = downloaded * 100 // total
                        if pct >= last_logged_pct + _PROGRESS_LOG_STEP_PCT:
                            last_logged_pct = pct
                            logger.info(
                                "%s: %d / %d MB (%d%%)",
                                dest.name,
                                downloaded // (1024 * 1024),
                                total // (1024 * 1024),
                                pct,
                            )
        tmp_dest.rename(dest)
    except BaseException:
        tmp_dest.unlink(missing_ok=True)
        raise
    return downloaded


def _download_with_retry(url: str, dest: Path) -> int:
    for attempt in range(_DOWNLOAD_MAX_ATTEMPTS):
        try:
            return _download_file(url, dest)
        except (urllib.error.URLError, OSError) as e:
            if attempt == _DOWNLOAD_MAX_ATTEMPTS - 1:
                raise
            wait = _DOWNLOAD_BACKOFF_BASE_SEC ** (attempt + 1)
            logger.warning(
                "%s: attempt %d/%d failed (%s), retrying in %ds",
                dest.name,
                attempt + 1,
                _DOWNLOAD_MAX_ATTEMPTS,
                e,
                wait,
            )
            time.sleep(wait)
    raise RuntimeError("unreachable")


def download_refseq(outdir: Path) -> Path:
    """Download reference data from NCBI FTP.

    Returns:
        Path to the downloads directory
    """
    downloads_dir = outdir / "downloads"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    meta_path = downloads_dir / DOWNLOAD_META_FILENAME
    meta = _load_meta(meta_path)

    for i, (filename, url) in enumerate(NCBI_DOWNLOAD_FILES.items()):
        if i > 0:
            time.sleep(_NCBI_REQUEST_INTERVAL_SEC)
        dest = downloads_dir / filename

        logger.info("%s: checking Last-Modified...", filename)
        remote_modified = _head_last_modified(url)

        cached = meta.get(filename, {})
        if (
            dest.exists()
            and cached.get("last_modified")
            and cached["last_modified"] == remote_modified
        ):
            logger.info("%s: cache valid, skipping", filename)
            continue

        logger.info("%s: starting download", filename)
        size = _download_with_retry(url, dest)

        try:
            post_modified = _head_last_modified(url)
        except (urllib.error.URLError, OSError) as e:
            logger.warning("%s: post-download Last-Modified check failed (%s)", filename, e)
            post_modified = None

        if post_modified is not None and remote_modified != post_modified:
            logger.warning(
                "%s: Last-Modified changed during download "
                "(start: %s -> end: %s); file may be incomplete",
                filename,
                remote_modified,
                post_modified,
            )

        meta[filename] = {
            "last_modified": remote_modified,
            "size": size,
        }
        _save_meta(meta_path, meta)
        logger.info("%s: done (%d MB)", filename, size // (1024 * 1024))

    return downloads_dir


def _extract_taxdump(archive_bytes: bytes, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=BytesIO(archive_bytes), mode="r:gz") as tf:
        for member in tf.getmembers():
            if member.name in _TAXDUMP_REQUIRED_FILES:
                tf.extract(member, dest_dir, filter="data")
    missing = [f for f in _TAXDUMP_REQUIRED_FILES if not (dest_dir / f).exists()]
    if missing:
        raise FileNotFoundError(f"taxdump archive is missing required files: {missing}")


def has_taxdump(outdir: Path) -> bool:
    taxdump_dir = outdir / "taxdump"
    return all((taxdump_dir / f).exists() for f in _TAXDUMP_REQUIRED_FILES)


def download_taxdump(outdir: Path) -> Path:
    """Download and extract NCBI taxdump into outdir/taxdump.

    Returns:
        taxdump directory path containing .dmp files
    """
    taxdump_dir = outdir / "taxdump"
    meta_path = taxdump_dir / DOWNLOAD_META_FILENAME
    meta = _load_meta(meta_path)

    logger.info("taxdump: checking Last-Modified...")
    remote_modified = _head_last_modified(_NCBI_TAXDUMP_URL)

    cached = meta.get("taxdump.tar.gz", {})
    has_all_files = all((taxdump_dir / f).exists() for f in _TAXDUMP_REQUIRED_FILES)
    if has_all_files and cached.get("last_modified") and cached["last_modified"] == remote_modified:
        logger.info("taxdump: cache valid, skipping")
        return taxdump_dir

    logger.info("taxdump: starting download")
    archive_dest = taxdump_dir / "taxdump.tar.gz"
    taxdump_dir.mkdir(parents=True, exist_ok=True)
    size = _download_with_retry(_NCBI_TAXDUMP_URL, archive_dest)

    logger.info("taxdump: extracting...")
    archive_bytes = archive_dest.read_bytes()
    _extract_taxdump(archive_bytes, taxdump_dir)
    archive_dest.unlink()

    meta["taxdump.tar.gz"] = {
        "last_modified": remote_modified,
        "size": size,
    }
    _save_meta(meta_path, meta)
    logger.info("taxdump: done (%d MB)", size // (1024 * 1024))

    return taxdump_dir
