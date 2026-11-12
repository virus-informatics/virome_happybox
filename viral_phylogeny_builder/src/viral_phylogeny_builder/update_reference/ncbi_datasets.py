# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""NCBI Datasets CLI wrappers for virus genome metadata packages."""

from __future__ import annotations

import hashlib
import logging
import subprocess
import zipfile
import zlib
from collections.abc import Callable
from pathlib import Path

logger = logging.getLogger(__name__)

# Field order defines the TSV columns emitted by `dataformat tsv virus-genome`
VIRUS_GENOME_METADATA_FIELDS = (
    "accession",
    "virus-tax-id",
    "virus-name",
    "host-name",
    "host-tax-id",
    "completeness",
    "length",
    "sourcedb",
)

# Display names dataformat writes for the fields above, verified against datasets 18.41.0
NCBI_ACCESSION_COLUMN = "Accession"
NCBI_VIRUS_TAXID_COLUMN = "Virus Taxonomic ID"
NCBI_VIRUS_NAME_COLUMN = "Virus Name"
NCBI_HOST_NAME_COLUMN = "Host Name"
NCBI_HOST_TAXID_COLUMN = "Host Taxonomic ID"
NCBI_COMPLETENESS_COLUMN = "Completeness"
NCBI_LENGTH_COLUMN = "Length"
NCBI_SOURCE_DATABASE_COLUMN = "Source database"

VIRUS_GENOME_METADATA_TSV_HEADER = (
    NCBI_ACCESSION_COLUMN,
    NCBI_VIRUS_TAXID_COLUMN,
    NCBI_VIRUS_NAME_COLUMN,
    NCBI_HOST_NAME_COLUMN,
    NCBI_HOST_TAXID_COLUMN,
    NCBI_COMPLETENESS_COLUMN,
    NCBI_LENGTH_COLUMN,
    NCBI_SOURCE_DATABASE_COLUMN,
)

_NCBI_PACKAGE_MD5_FILENAME = "md5sum.txt"
_MD5_READ_CHUNK_BYTES = 1024 * 1024


def _parse_md5sum_file(text: str) -> dict[str, str]:
    """Parse `md5sum.txt` entries of the form `<md5>  <path within package>`."""
    checksums: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            continue
        digest, member = parts
        checksums[member.strip().lstrip("*")] = digest
    return checksums


def _md5_of_member(package: zipfile.ZipFile, member: str) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with package.open(member) as handle:
        while chunk := handle.read(_MD5_READ_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def is_valid_package(zip_path: Path) -> bool:
    """Verify a NCBI Datasets package against the md5sum.txt it ships with."""
    if not zip_path.exists() or zip_path.stat().st_size == 0:
        return False

    try:
        with zipfile.ZipFile(zip_path) as package:
            members = set(package.namelist())
            if _NCBI_PACKAGE_MD5_FILENAME not in members:
                logger.warning("Package has no %s: %s", _NCBI_PACKAGE_MD5_FILENAME, zip_path)
                return False

            checksums = _parse_md5sum_file(package.read(_NCBI_PACKAGE_MD5_FILENAME).decode())
            if not checksums:
                logger.warning(
                    "Package %s has no entries: %s",
                    _NCBI_PACKAGE_MD5_FILENAME,
                    zip_path,
                )
                return False

            for member, expected_md5 in checksums.items():
                if member not in members:
                    logger.warning("Package is missing %s: %s", member, zip_path)
                    return False
                if _md5_of_member(package, member) != expected_md5:
                    logger.warning("MD5 mismatch for %s: %s", member, zip_path)
                    return False
    # A corrupted deflate stream surfaces as zlib.error / EOFError only while reading a member
    except (OSError, zipfile.BadZipFile, zlib.error, EOFError) as e:
        logger.warning("Package is not a readable ZIP (%s): %s", e, zip_path)
        return False

    return True


_PARTIAL_SUFFIX = ".partial"


def download_validated_package(
    download_fn: Callable[[Path], None],
    package_path: Path,
) -> None:
    """Download a NCBI Datasets package, keeping the previous one until validation passes.

    1. Download into a .partial sibling so that package_path is never touched mid-transfer.
    2. Verify the .partial file against the md5sum.txt inside the ZIP.
       e.g. metadata.zip → metadata.zip.partial
    3. Atomically replace package_path only after validation succeeds.
    4. On any failure (network, validation, interrupt), remove .partial and re-raise
       so that the caller's cached package_path survives intact.
    """
    partial_path = package_path.with_name(package_path.name + _PARTIAL_SUFFIX)
    try:
        download_fn(partial_path)
        if not is_valid_package(partial_path):
            raise ValueError(
                f"Downloaded package failed MD5 validation: {package_path}",
            )
        # Validation passed — swap in the new package.
        partial_path.replace(package_path)
    except BaseException:
        # Download or validation failed — discard the partial, keep the old package.
        partial_path.unlink(missing_ok=True)
        raise


def download_virus_genome_metadata(
    virus_family: str,
    zip_path: Path,
    *,
    datasets_cmd: str = "datasets",
) -> None:
    """Download the metadata-only virus genome package of complete sequences for a virus family.

    stdout and stderr are inherited so that the progress display of datasets stays
    visible: a large family exceeds 100 MB and takes minutes to download.
    """
    cmd = [
        datasets_cmd,
        "download",
        "virus",
        "genome",
        "taxon",
        virus_family,
        "--include",
        "none",
        # Reference selection uses only COMPLETE records
        "--complete-only",
        # is_valid_package checks the package with MD5 and CRC-32
        "--fast-zip-validation",
        "--filename",
        str(zip_path),
    ]
    logger.info("Running datasets download virus genome taxon")

    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)


def download_virus_genomes_by_accessions(
    accessions_path: Path,
    zip_path: Path,
    *,
    datasets_cmd: str = "datasets",
) -> None:
    """Download the genome package for the nucleotide accessions listed in a file.

    The accessions are handed over as a file rather than as arguments: a family yields one
    accession per virus taxid, which is thousands of them for the larger families.

    stdout and stderr are inherited so that the progress display of datasets stays visible.
    """
    cmd = [
        datasets_cmd,
        "download",
        "virus",
        "genome",
        "accession",
        "--inputfile",
        str(accessions_path),
        "--include",
        "genome",
        "--filename",
        str(zip_path),
    ]
    logger.info("Running datasets download virus genome accession")

    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)


def convert_package_to_tsv(
    zip_path: Path,
    tsv_path: Path,
    *,
    dataformat_cmd: str = "dataformat",
) -> None:
    """Convert a metadata package to TSV via dataformat.

    stdout is streamed to disk rather than buffered: a large family yields
    over a million records.
    """
    cmd = [
        dataformat_cmd,
        "tsv",
        "virus-genome",
        "--package",
        str(zip_path),
        "--fields",
        ",".join(VIRUS_GENOME_METADATA_FIELDS),
    ]
    logger.info("Running dataformat tsv virus-genome")

    partial_path = tsv_path.with_name(tsv_path.name + ".partial")
    try:
        # binary mode: dataformat writes to the file descriptor, so the text wrapper is unused
        with partial_path.open("wb") as handle:
            result = subprocess.run(
                cmd,
                stdout=handle,
                stderr=subprocess.PIPE,
                encoding="utf-8",
            )
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd, None, result.stderr)
        partial_path.replace(tsv_path)
    except BaseException:
        partial_path.unlink(missing_ok=True)
        raise
