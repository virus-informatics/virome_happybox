# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Select one representative genome per virus taxid from the fetched NCBI metadata."""

from __future__ import annotations

import logging
from pathlib import Path

from viral_phylogeny_builder._virus_genome_metadata_columns import (
    METADATA_ACCESSION_COLUMN,
)
from viral_phylogeny_builder.update_reference.host_classification import (
    classify_host_lineage,
    read_host_manual_key,
)
from viral_phylogeny_builder.update_reference.lineage_metadata import (
    read_filtered_metadata,
)
from viral_phylogeny_builder.update_reference.taxonkit_lineage import (
    TAXONKIT_REFORMATTED_LINEAGE_GENUS_INDEX,
    run_taxonkit_lineage,
)

logger = logging.getLogger(__name__)

# Holds the taxid lists and the raw taxonkit output, which are inputs to the lineage columns
# rather than outputs of the step.
_TAXONKIT_TMP_DIR_NAME = "tmp_taxonkit"

UNKNOWN_VIRUS_GENUS = "Unknown"

# The workdir already carries the virus family, so the names drop the prefix the archive
# script needed to tell its per-family outputs apart inside one directory.
LINEAGE_CSV_FILENAME = "lineage.csv"
COMPLETE_LENGTH_FILTERED_CSV_FILENAME = "complete_length_filtered.csv"
REFERENCE_METADATA_CSV_FILENAME = "reference_metadata.csv"
REFERENCE_ACCESSIONS_TXT_FILENAME = "reference_accessions.txt"


def virus_genus_from_reformatted_lineage(reformatted_lineage: str) -> str:
    """Read the genus off one reformatted lineage, or UNKNOWN_VIRUS_GENUS when it is empty."""
    reformatted_lineage = reformatted_lineage.strip()
    if not reformatted_lineage:
        return UNKNOWN_VIRUS_GENUS

    genus = reformatted_lineage.split(";")[TAXONKIT_REFORMATTED_LINEAGE_GENUS_INDEX].strip()
    if not genus:
        return UNKNOWN_VIRUS_GENUS

    return genus.split()[0]


def completeness_score(completeness_upper: str) -> int:
    """Rank one stripped, upper-cased Completeness value for the reference sort, lowest first.

    A COMPLETE genome wins, then any other value NCBI gave, then a row that carries none.
    """
    if completeness_upper == "COMPLETE":
        return 0
    if completeness_upper:
        return 1
    return 2


def select_reference_metadata(
    metadata_csv: Path,
    workdir: Path,
    *,
    exclude_virus_taxid: int,
    taxdump_dir: Path,
    taxonkit_cmd: str = "taxonkit",
) -> None:
    """Write the four selection outputs into workdir.

    lineage.csv, complete_length_filtered.csv, reference_metadata.csv and
    reference_accessions.txt.

    Caching is the caller's, keyed on the whole subcommand.
    """
    metadata = read_filtered_metadata(metadata_csv, exclude_virus_taxid=exclude_virus_taxid)
    tmp_dir = workdir / _TAXONKIT_TMP_DIR_NAME

    virus_lineage_df = run_taxonkit_lineage(
        metadata["virus_taxid"].tolist(),
        taxdump_dir,
        tmp_dir,
        "virus",
        taxonkit_cmd=taxonkit_cmd,
    ).add_prefix("virus_")
    host_lineage_df = run_taxonkit_lineage(
        metadata["host_taxid"].tolist(),
        taxdump_dir,
        tmp_dir,
        "host",
        taxonkit_cmd=taxonkit_cmd,
    ).add_prefix("host_")

    metadata = metadata.merge(virus_lineage_df, on="virus_taxid", how="left")
    metadata["virus_genus"] = metadata["virus_reformatted_lineage"].map(
        virus_genus_from_reformatted_lineage,
    )
    metadata = metadata.merge(host_lineage_df, on="host_taxid", how="left")

    # The host lineage is split into nodes before the manual key sees it,
    # so a row that joined no lineage has to reach it as a string.
    # TODO: VII-150 drop the taxon_name and rank columns once the remaining steps are ported.
    # Neither side of those two pairs is read, and a missing value reaches the CSV as the
    # same empty field an empty string does.
    for column in (
        "virus_lineage",
        "virus_taxon_name",  # NOTE: unused later
        "virus_rank",  # NOTE: unused later
        "virus_reformatted_lineage",
        "host_lineage",
        "host_taxon_name",  # NOTE: unused later
        "host_rank",  # NOTE: unused later
        "host_reformatted_lineage",
    ):
        metadata[column] = metadata[column].fillna("")

    host_rules = read_host_manual_key()
    classified_hosts = [
        classify_host_lineage(host_lineage, host_rules) for host_lineage in metadata["host_lineage"]
    ]

    metadata["manual_host_class"] = [
        classification.manual_host_class for classification in classified_hosts
    ]
    metadata["manual_host_taxonomy"] = [
        classification.manual_host_taxonomy for classification in classified_hosts
    ]
    metadata["manual_host"] = metadata["manual_host_class"]
    metadata["host_priority"] = [
        classification.host_priority for classification in classified_hosts
    ]
    metadata["completeness_priority"] = metadata["completeness_upper"].map(completeness_score)

    complete = metadata.loc[metadata["completeness_upper"] == "COMPLETE"].copy()
    if complete.empty:
        observed = sorted(metadata["completeness_upper"].unique())
        raise ValueError(
            f"No COMPLETE records were found. Observed Completeness values: {observed}",
        )

    complete["complete_genome_length_median"] = complete.groupby("virus_taxid")[
        "length_numeric"
    ].transform("median")
    complete["length_threshold_50pct_complete_median"] = (
        complete["complete_genome_length_median"] * 0.5
    )

    length_filtered = complete.loc[
        complete["length_numeric"] >= complete["length_threshold_50pct_complete_median"]
    ].copy()
    logger.info(
        "COMPLETE rows: %d, after the length filter: %d",
        len(complete),
        len(length_filtered),
    )

    length_filtered = length_filtered.sort_values(
        by=[
            "virus_taxid",
            "host_priority",
            "completeness_priority",
            "length_numeric",
            "accession_base",
        ],
        ascending=[True, True, True, False, True],
    )

    reference_metadata = length_filtered.groupby("virus_taxid", as_index=False).head(1).copy()
    logger.info("Reference rows: %d", len(reference_metadata))

    lineage_csv = workdir / LINEAGE_CSV_FILENAME
    complete_length_filtered_csv = workdir / COMPLETE_LENGTH_FILTERED_CSV_FILENAME
    reference_metadata_csv = workdir / REFERENCE_METADATA_CSV_FILENAME
    reference_accessions_txt = workdir / REFERENCE_ACCESSIONS_TXT_FILENAME

    metadata.to_csv(lineage_csv, index=False)
    length_filtered.to_csv(complete_length_filtered_csv, index=False)
    reference_metadata.to_csv(reference_metadata_csv, index=False)

    # The accession keeps its version here, unlike accession_base: this list is what the
    # genome download is given.
    reference_metadata[
        METADATA_ACCESSION_COLUMN
    ].str.strip().drop_duplicates().sort_values().to_csv(
        reference_accessions_txt,
        index=False,
        header=False,
    )

    logger.info(
        "COMPLETE median length summary:\n%s",
        complete["complete_genome_length_median"].describe().to_string(),
    )
    logger.info(
        "Filtered length summary:\n%s",
        length_filtered["length_numeric"].describe().to_string(),
    )

    logger.info("Lineage metadata CSV: %s", lineage_csv)
    logger.info("Length-filtered metadata CSV: %s", complete_length_filtered_csv)
    logger.info("Reference metadata CSV: %s", reference_metadata_csv)
    logger.info("Reference accessions: %s", reference_accessions_txt)
