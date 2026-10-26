# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Fetch NCBI virus genome metadata for a virus family (datasets -> dataformat -> CSV)."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from viral_phylogeny_builder._commands import (
    UPDATE_REFERENCE_COMMAND,
    UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
)
from viral_phylogeny_builder._pipeline import (
    mark_step_done,
    should_skip,
    step_done_marker,
    step_is_done,
    step_workdir,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools
from viral_phylogeny_builder.update_reference.ncbi_datasets import (
    VIRUS_GENOME_METADATA_TSV_HEADER,
    convert_package_to_tsv,
    download_validated_package,
    download_virus_genome_metadata,
    is_valid_package,
)
from viral_phylogeny_builder.update_reference.reference_selection import (
    REFERENCE_ACCESSIONS_TXT_FILENAME,
    REFERENCE_METADATA_CSV_FILENAME,
    select_reference_metadata,
)
from viral_phylogeny_builder.update_reference.taxdump import (
    TAXDUMP_DIR_NAME,
    download_taxdump,
    has_taxdump,
    validate_taxdump_dir,
)

logger = logging.getLogger(__name__)

_METADATA_PACKAGE_FILENAME = "metadata.zip"
_METADATA_TSV_FILENAME = "virus_genome_metadata.tsv"
FETCH_METADATA_CSV_FILENAME = "virus_genome_metadata.csv"

_FETCH_METADATA_REQUIRED_TOOLS = ("datasets", "dataformat", "taxonkit")

# Checked alongside the marker: these two are what the later subcommands read.
_FETCH_METADATA_STEP_OUTPUT_FILENAMES = (
    REFERENCE_METADATA_CSV_FILENAME,
    REFERENCE_ACCESSIONS_TXT_FILENAME,
)


def _fetch_metadata_step_outputs(workdir: Path) -> list[Path]:
    return [workdir / filename for filename in _FETCH_METADATA_STEP_OUTPUT_FILENAMES]


def _count_csv_data_rows(csv_path: Path) -> int:
    """Count the data rows of an existing metadata CSV, excluding the header."""
    with csv_path.open(encoding="utf-8", newline="") as source:
        return sum(1 for _ in csv.reader(source)) - 1


def _convert_tsv_to_csv(tsv_path: Path, csv_path: Path) -> int:
    """Convert the metadata TSV to CSV row by row and return the data row count."""
    partial_path = csv_path.with_name(csv_path.name + ".partial")
    data_row_count = 0
    try:
        with (
            tsv_path.open(encoding="utf-8", newline="") as source,
            partial_path.open("w", encoding="utf-8", newline="") as dest,
        ):
            # csv.writer defaults to the excel dialect, which terminates rows with CRLF
            writer = csv.writer(dest, lineterminator="\n")
            reader = csv.reader(source, delimiter="\t")

            header = next(reader, None)
            if header is None:
                raise ValueError(f"Metadata TSV is empty: {tsv_path}")
            if tuple(header) != VIRUS_GENOME_METADATA_TSV_HEADER:
                raise ValueError(
                    f"Unexpected metadata TSV header: {header} "
                    f"(expected {list(VIRUS_GENOME_METADATA_TSV_HEADER)})",
                )
            writer.writerow(header)

            for row in reader:
                writer.writerow(row)
                data_row_count += 1

        if data_row_count == 0:
            raise ValueError(f"Metadata TSV has no records: {tsv_path}")

        partial_path.replace(csv_path)
    except BaseException:
        partial_path.unlink(missing_ok=True)
        raise

    return data_row_count


def run_update_reference_fetch_metadata(
    virus_family: str,
    output_dir: Path,
    *,
    exclude_virus_taxid: int,
    taxdump_dir: Path | None = None,
    update_taxdump: bool = False,
    force: bool = False,
    tool_config: dict[str, str] | None = None,
) -> None:
    """Fetch NCBI virus genome metadata and select the representative sequences.

    taxdump_dir defaults to the managed copy under output_dir, which is downloaded when
    missing or when update_taxdump is set; a supplied directory is only checked for the
    files taxonkit reads.
    """
    resolved = resolve_tools(_FETCH_METADATA_REQUIRED_TOOLS, tool_config)
    validate_virus_family(virus_family)

    workdir = step_workdir(
        output_dir,
        UPDATE_REFERENCE_COMMAND,
        UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
        virus_family,
    )
    workdir.mkdir(parents=True, exist_ok=True)

    package_path = workdir / _METADATA_PACKAGE_FILENAME
    if not force and is_valid_package(package_path):
        logger.info("Skipping metadata download (cached): %s", package_path)
        package_downloaded = False
    else:
        download_validated_package(
            lambda partial: download_virus_genome_metadata(
                virus_family,
                partial,
                datasets_cmd=resolved["datasets"],
            ),
            package_path,
        )
        package_downloaded = True

    # A newly downloaded package leaves any TSV/CSV derived from the previous one stale
    rebuild_derived = package_downloaded

    # After the package check rather than at the entry, so that a package NCBI has replaced
    # still reaches the steps that derive from it.
    marker = step_done_marker(
        workdir,
        UPDATE_REFERENCE_COMMAND,
        UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
    )
    step_outputs = _fetch_metadata_step_outputs(workdir)

    if not force and not rebuild_derived and step_is_done(marker, step_outputs):
        logger.info(
            "Skipping %s (cached): %s",
            UPDATE_REFERENCE_FETCH_METADATA_SUBCOMMAND,
            workdir,
        )
        return

    # Dropped before the rest of the work restarts: leaving it would let a run interrupted
    # here be read as finished, with outputs from two different runs side by side.
    marker.unlink(missing_ok=True)

    tsv_path = workdir / _METADATA_TSV_FILENAME
    if not should_skip(tsv_path, force=rebuild_derived, step_name="dataformat tsv virus-genome"):
        convert_package_to_tsv(package_path, tsv_path, dataformat_cmd=resolved["dataformat"])

    csv_path = workdir / FETCH_METADATA_CSV_FILENAME
    if should_skip(csv_path, force=rebuild_derived, step_name="TSV to CSV conversion"):
        data_row_count = _count_csv_data_rows(csv_path)
    else:
        data_row_count = _convert_tsv_to_csv(tsv_path, csv_path)
    logger.info("Metadata records: %d", data_row_count)

    if taxdump_dir is None:
        managed_taxdump_dir = output_dir / TAXDUMP_DIR_NAME
        if update_taxdump or not has_taxdump(managed_taxdump_dir):
            taxdump_dir = download_taxdump(managed_taxdump_dir)
        else:
            logger.info("Taxdump update skipped (use --update-taxdump to update)")
            taxdump_dir = managed_taxdump_dir
    else:
        validate_taxdump_dir(taxdump_dir)

    select_reference_metadata(
        csv_path,
        workdir,
        exclude_virus_taxid=exclude_virus_taxid,
        taxdump_dir=taxdump_dir,
        taxonkit_cmd=resolved["taxonkit"],
    )

    mark_step_done(marker, step_outputs)
    logger.info("Fetch metadata complete: %s", workdir)
