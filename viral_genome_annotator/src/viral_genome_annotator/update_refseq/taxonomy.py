# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Resolve NCBI taxids to Family/Genus/Species via taxonkit."""

from __future__ import annotations

import csv
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# taxonkit reformat -f format placeholders: {f}=Family, {g}=Genus, {s}=Species
# See: taxonkit reformat --help
_TAXONKIT_REFORMAT_TEMPLATE = "{f}\t{g}\t{s}"

# Output of `taxonkit lineage | taxonkit reformat -I 1 -f "{f}\t{g}\t{s}"`:
#   col0: taxid  col1: lineage  col2: family  col3: genus  col4: species
# Verified with:
#   echo "1980427" | taxonkit lineage | taxonkit reformat -I 1 -f "{f}\t{g}\t{s}"
#   => 1980427\tViruses;...;Emaravirus fici\tFimoviridae\tEmaravirus\tEmaravirus fici
_TAXONKIT_OUTPUT_COLUMNS = 5


@dataclass
class TaxonomyEntry:
    """Resolved taxonomy for a single taxid."""

    taxid: int
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None


def _run_taxonkit(
    taxids: list[int],
    taxonkit_data_dir: Path,
) -> str:
    taxid_input = "\n".join(str(t) for t in taxids)
    lineage = subprocess.run(
        ["taxonkit", "lineage", "--data-dir", str(taxonkit_data_dir)],
        input=taxid_input,
        capture_output=True,
        text=True,
        check=True,
    )
    reformat = subprocess.run(
        [
            "taxonkit",
            "reformat",
            "--data-dir",
            str(taxonkit_data_dir),
            "-I",
            "1",
            "-f",
            _TAXONKIT_REFORMAT_TEMPLATE,
        ],
        input=lineage.stdout,
        capture_output=True,
        text=True,
        check=True,
    )
    return reformat.stdout


def _parse_taxonkit_output(stdout: str) -> dict[int, TaxonomyEntry]:
    entries: dict[int, TaxonomyEntry] = {}
    for line in stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < _TAXONKIT_OUTPUT_COLUMNS:
            continue
        try:
            taxid = int(parts[0])
        except ValueError:
            logger.warning("taxonkit output has non-integer taxid: %r", parts[0])
            continue
        family = parts[2] or None
        genus = parts[3] or None
        species = parts[4] or None
        entries[taxid] = TaxonomyEntry(
            taxid=taxid,
            family=family,
            genus=genus,
            species=species,
        )
    return entries


def resolve_taxonomy(
    taxids: set[int],
    taxonkit_data_dir: Path,
) -> dict[int, TaxonomyEntry]:
    """Resolve a set of taxids to Family/Genus/Species via taxonkit.

    Returns a dict keyed by taxid.  Taxids that taxonkit cannot resolve
    are included with all-None fields.
    """
    if not taxids:
        return {}
    sorted_taxids = sorted(taxids)
    logger.info("Resolving %d unique taxids via taxonkit", len(sorted_taxids))
    stdout = _run_taxonkit(sorted_taxids, taxonkit_data_dir)
    entries = _parse_taxonkit_output(stdout)

    unresolved = taxids - entries.keys()
    for taxid in sorted(unresolved):
        logger.warning("taxonkit could not resolve taxid %d", taxid)
        entries[taxid] = TaxonomyEntry(taxid=taxid)

    return entries


def write_taxonomy_tsv(
    entries: dict[int, TaxonomyEntry],
    dest: Path,
) -> None:
    """Write taxonomy entries as a TSV file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    sorted_entries = sorted(entries.values(), key=lambda e: e.taxid)
    with dest.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "Family", "Genus", "Species"])
        for entry in sorted_entries:
            writer.writerow(
                [
                    entry.taxid,
                    entry.family or "",
                    entry.genus or "",
                    entry.species or "",
                ],
            )
    logger.info("Wrote %d entries to %s", len(sorted_entries), dest)
