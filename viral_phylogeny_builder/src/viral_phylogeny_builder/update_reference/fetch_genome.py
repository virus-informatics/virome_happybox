# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Download the reference genomes for the accessions selected by fetch-metadata."""

from __future__ import annotations

import hashlib
import logging
import shutil
import zipfile
from pathlib import Path

from viral_phylogeny_builder._commands import (
    UPDATE_REFERENCE_COMMAND,
    UPDATE_REFERENCE_FETCH_GENOME_SUBCOMMAND,
    UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
)
from viral_phylogeny_builder._fasta import iter_fasta
from viral_phylogeny_builder._pipeline import (
    mark_step_done,
    step_done_marker,
    step_is_done,
    step_workdir,
    validate_input_file,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools
from viral_phylogeny_builder.update_reference.ncbi_datasets import (
    download_validated_package,
    download_virus_genomes_by_accessions,
    is_valid_package,
)
from viral_phylogeny_builder.update_reference.reference_selection import (
    REFERENCE_ACCESSIONS_TXT_FILENAME,
)

logger = logging.getLogger(__name__)

REFERENCE_GENOMES_ZIP_FILENAME = "reference_genomes.zip"
REFERENCE_GENOMES_DIR_NAME = "reference_genomes"
# The extracted directory, followed by the path datasets lays the genomic sequences out at
# inside the package. Exported so that the steps reading the FASTA do not rebuild the path.
REFERENCE_GENOME_FASTA_PATH_RELATIVE_TO_WORKDIR = (
    Path(REFERENCE_GENOMES_DIR_NAME) / "ncbi_dataset" / "data" / "genomic.fna"
)

_FETCH_GENOME_REQUIRED_TOOLS = ("datasets",)
_ACCESSIONS_SHA256_FILENAME = "reference_accessions.sha256"
_SHA256_READ_CHUNK_BYTES = 1024 * 1024


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_SHA256_READ_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def _accessions_changed(accessions_path: Path, sha256_path: Path) -> bool:
    """Return True when the accession list differs from the last recorded digest.

    Returns False when no digest has been recorded yet: the first run after
    upgrading to this check should not force a redundant re-download.
    """
    if not sha256_path.exists():
        return False
    return _file_sha256(accessions_path) != sha256_path.read_text().strip()


def _record_accessions_sha256(accessions_path: Path, sha256_path: Path) -> None:
    sha256_path.write_text(_file_sha256(accessions_path) + "\n")


def _read_accessions(accessions_path: Path) -> set[str]:
    with accessions_path.open() as handle:
        return {line.strip() for line in handle if line.strip()}


def _warn_missing_accessions(
    requested: set[str],
    genome_fasta: Path,
) -> None:
    """Warn about accessions requested but absent from the downloaded FASTA."""
    fetched = {header.split()[0] for header, _seq in iter_fasta(genome_fasta)}
    missing = requested - fetched
    if missing:
        logger.warning(
            "%d of %d requested accessions missing from downloaded FASTA: %s",
            len(missing),
            len(requested),
            ", ".join(sorted(missing)),
        )


def _reference_accessions_path(output_dir: Path, virus_family: str) -> Path:
    """Path of the accession list fetch-metadata wrote for this virus family."""
    return (
        step_workdir(
            output_dir,
            UPDATE_REFERENCE_COMMAND,
            UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
            virus_family,
        )
        / REFERENCE_ACCESSIONS_TXT_FILENAME
    )


_EXTRACT_STAGING_SUFFIX = ".tmp"


def _extract_package(zip_path: Path, extract_dir: Path) -> None:
    """Replace the extracted directory with the contents of the package.

    Extracts into a staging directory first so that an interruption never leaves
    a half-extracted result where the cache check would accept it.
    """
    staging_dir = extract_dir.with_name(extract_dir.name + _EXTRACT_STAGING_SUFFIX)
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True)
    try:
        with zipfile.ZipFile(zip_path) as package:
            package.extractall(staging_dir)
    except BaseException:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    staging_dir.rename(extract_dir)


def run_update_reference_fetch_genome(
    virus_family: str,
    output_dir: Path,
    *,
    force: bool = False,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Download the genomes fetch-metadata selected and return the path to their FASTA."""
    resolved = resolve_tools(_FETCH_GENOME_REQUIRED_TOOLS, tool_config)
    validate_virus_family(virus_family)

    accessions_path = _reference_accessions_path(output_dir, virus_family)
    validate_input_file(accessions_path, "Reference accession list")
    requested_accessions = _read_accessions(accessions_path)
    logger.info(
        "Reference accessions: %d (%s)",
        len(requested_accessions),
        accessions_path,
    )

    workdir = step_workdir(
        output_dir,
        UPDATE_REFERENCE_COMMAND,
        UPDATE_REFERENCE_FETCH_GENOME_SUBCOMMAND,
        virus_family,
    )
    workdir.mkdir(parents=True, exist_ok=True)

    package_path = workdir / REFERENCE_GENOMES_ZIP_FILENAME
    extract_dir = workdir / REFERENCE_GENOMES_DIR_NAME
    genome_fasta = workdir / REFERENCE_GENOME_FASTA_PATH_RELATIVE_TO_WORKDIR

    marker = step_done_marker(
        workdir,
        UPDATE_REFERENCE_COMMAND,
        UPDATE_REFERENCE_FETCH_GENOME_SUBCOMMAND,
    )
    step_outputs = [genome_fasta]
    sha256_path = workdir / _ACCESSIONS_SHA256_FILENAME

    if not force:
        if _accessions_changed(accessions_path, sha256_path):
            logger.info(
                "Reference accessions changed since last download; re-downloading",
            )
            marker.unlink(missing_ok=True)
            sha256_path.unlink(missing_ok=True)
            package_path.unlink(missing_ok=True)
        elif step_is_done(marker, step_outputs):
            logger.info(
                "Skipping %s (cached): %s",
                UPDATE_REFERENCE_FETCH_GENOME_SUBCOMMAND,
                workdir,
            )
            if not sha256_path.exists():
                _record_accessions_sha256(accessions_path, sha256_path)
            return genome_fasta

    # Dropped before the rest of the work restarts: leaving it would let a run interrupted
    # here be read as finished.
    marker.unlink(missing_ok=True)

    if not force and is_valid_package(package_path):
        logger.info("Skipping genome download (cached): %s", package_path)
    else:
        download_validated_package(
            lambda partial: download_virus_genomes_by_accessions(
                accessions_path,
                partial,
                datasets_cmd=resolved["datasets"],
            ),
            package_path,
        )

    _extract_package(package_path, extract_dir)
    _warn_missing_accessions(requested_accessions, genome_fasta)

    mark_step_done(marker, step_outputs)
    _record_accessions_sha256(accessions_path, sha256_path)
    logger.info("Fetch genome complete: %s", genome_fasta)
    return genome_fasta
