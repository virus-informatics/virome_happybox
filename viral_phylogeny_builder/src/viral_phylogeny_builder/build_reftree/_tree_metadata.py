# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Prepare tree_metadata.csv that virome_tree_plot_common.R consumes."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from viral_phylogeny_builder._accession import strip_accession_version
from viral_phylogeny_builder._fasta import iter_fasta

logger = logging.getLogger(__name__)

TREE_METADATA_COLUMNS = ("Accession", "Genus", "manual_host", "Organism_Name", "GenBank_RefSeq")


def _extract_tip_labels(alignment_fasta: Path) -> list[str]:
    """Read FASTA headers and return unique tip labels (version-stripped, order-preserved)."""
    seen: set[str] = set()
    labels: list[str] = []
    for header, _seq in iter_fasta(alignment_fasta):
        record_id = header.split()[0]
        label = strip_accession_version(record_id)
        if label and label not in seen:
            labels.append(label)
            seen.add(label)
    return labels


def prepare_tree_metadata(
    alignment_fasta: Path,
    rdrp_metadata_csv: Path,
    output_csv: Path,
) -> None:
    """Build the tree_metadata.csv from alignment FASTA headers and RdRp metadata."""
    tip_labels = _extract_tip_labels(alignment_fasta)
    if not tip_labels:
        raise ValueError(f"No FASTA headers found in {alignment_fasta}")

    df_metadata = pd.read_csv(rdrp_metadata_csv, dtype=str).fillna("")
    df_metadata["accession_key"] = df_metadata["Accession"].map(strip_accession_version)

    acc_set = set(tip_labels)
    df_metadata = df_metadata[df_metadata["accession_key"].isin(acc_set)]

    # NCBI Datasets uses different column names from NCBI Virus
    df_metadata = df_metadata.rename(
        columns={
            "virus_genus": "Genus",
            "Virus Name": "Organism_Name",
            "Source database": "GenBank_RefSeq",
        },
    )

    df_out = df_metadata[list(TREE_METADATA_COLUMNS) + ["accession_key"]].drop_duplicates()
    df_out = df_out.drop_duplicates(subset=["accession_key"], keep="first")
    df_out = df_out.drop(columns=["accession_key"])

    missing = acc_set - set(df_out["Accession"].map(strip_accession_version))
    if missing:
        logger.warning(
            "%d tip labels not found in metadata (e.g. %s)",
            len(missing),
            ", ".join(sorted(missing)[:5]),
        )

    df_out.to_csv(output_csv, index=False)
    logger.info("Tree metadata: %d rows written to %s", len(df_out), output_csv)
