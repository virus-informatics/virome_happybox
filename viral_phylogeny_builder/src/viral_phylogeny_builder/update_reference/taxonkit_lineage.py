# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Resolve taxids to lineages via taxonkit lineage + taxonkit reformat."""

from __future__ import annotations

import logging
import subprocess
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Status codes that mean the taxid has no lineage: -1 for one not found in the database, 0
# for a deleted one. Any other status code is a resolved taxid: the input taxid itself, or
# the new taxid when it was merged. See the Output section of `taxonkit lineage --help`.
_TAXONKIT_UNRESOLVED_STATUS_CODES = ("-1", "0")

# taxonkit reformat -f placeholders: {k}=superkingdom, {p}=phylum, {c}=class,
# {o}=order, {f}=family, {g}=genus, {s}=species.
# Same value as the taxonkit v0.20.0 default, pinned rather than left implicit because
# consumers address the reformatted lineage by position (the genus is the second-to-last
# field), so a change of the default would silently shift which rank they read.
_TAXONKIT_REFORMAT_RANK_FORMAT = "{k};{p};{c};{o};{f};{g};{s}"

# Derived from the format rather than written out, so that editing the format moves the genus
# position with it instead of leaving consumers reading the rank next to it.
TAXONKIT_REFORMATTED_LINEAGE_GENUS_INDEX = _TAXONKIT_REFORMAT_RANK_FORMAT.split(";").index("{g}")

# Columns of the two subcommands, in the order they document under Output: the input line,
# the status code, the lineage, the name, the rank, then the reformatted lineage that
# reformat appends. The rows are split positionally, so this order decides which field each
# column takes.
TAXONKIT_TAXID_COLUMN = "taxid"
TAXONKIT_STATUS_CODE_COLUMN = "status_code"
TAXONKIT_LINEAGE_COLUMN = "lineage"
TAXONKIT_TAXON_NAME_COLUMN = "taxon_name"
TAXONKIT_RANK_COLUMN = "rank"
TAXONKIT_REFORMATTED_LINEAGE_COLUMN = "reformatted_lineage"
TAXONKIT_OUTPUT_COLUMNS = (
    TAXONKIT_TAXID_COLUMN,
    TAXONKIT_STATUS_CODE_COLUMN,
    TAXONKIT_LINEAGE_COLUMN,
    TAXONKIT_TAXON_NAME_COLUMN,
    TAXONKIT_RANK_COLUMN,
    TAXONKIT_REFORMATTED_LINEAGE_COLUMN,
)


def _run_taxonkit(cmd: list[str], *, label: str, stdin_text: str | None = None) -> str:
    """Run one taxonkit subcommand and return its stdout."""
    result = subprocess.run(
        cmd,
        input=stdin_text,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        subcommand = cmd[1]
        raise RuntimeError(
            f"taxonkit {subcommand} failed for {label} (exit {result.returncode}): {result.stderr}",
        )
    return result.stdout


def _parse_taxonkit_output(output: str, *, prefix: str, lineage_file: Path) -> list[list[str]]:
    """Split the taxonkit output into rows, rejecting one that is not fully formed.

    A short row is not padded: an empty column would be indistinguishable from a lineage
    taxonkit returned as empty, which it does for every unresolved taxid.
    """
    rows = []
    for line_number, line in enumerate(output.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != len(TAXONKIT_OUTPUT_COLUMNS):
            raise RuntimeError(
                f"taxonkit output for {prefix} has {len(fields)} columns instead of"
                f" {len(TAXONKIT_OUTPUT_COLUMNS)} at line {line_number} of {lineage_file}",
            )
        rows.append(fields)
    return rows


def run_taxonkit_lineage(
    taxids: Iterable[str],
    taxdump_dir: Path,
    tmp_dir: Path,
    prefix: str,
    *,
    taxonkit_cmd: str = "taxonkit",
) -> pd.DataFrame:
    """Resolve taxids to lineages, one row per unique non-empty taxid.

    Raises RuntimeError when a taxonkit subcommand fails, when a row does not hold the
    documented number of columns, or when the row count does not match the taxid count.
    """
    tmp_dir.mkdir(parents=True, exist_ok=True)

    clean_taxids = sorted({stripped for taxid in taxids if (stripped := taxid.strip())})
    if not clean_taxids:
        return pd.DataFrame(columns=list(TAXONKIT_OUTPUT_COLUMNS), dtype=str)

    taxid_file = tmp_dir / f"{prefix}.taxids.txt"
    lineage_file = tmp_dir / f"{prefix}.lineage.tsv"

    taxid_file.write_text("\n".join(clean_taxids) + "\n")

    lineage_stdout = _run_taxonkit(
        [
            taxonkit_cmd,
            "lineage",
            "--show-status-code",
            "--show-name",
            "--show-rank",
            str(taxid_file),
            "--data-dir",
            str(taxdump_dir),
        ],
        label=prefix,
    )
    reformat_stdout = _run_taxonkit(
        [
            taxonkit_cmd,
            "reformat",
            "--taxid-field",
            "1",  # taxid filed is first index in taxonkit lineage
            "--format",
            _TAXONKIT_REFORMAT_RANK_FORMAT,  # Add default format in case of breaking change
            "--data-dir",
            str(taxdump_dir),
        ],
        label=prefix,
        stdin_text=lineage_stdout,
    )
    lineage_file.write_text(reformat_stdout)

    rows = _parse_taxonkit_output(reformat_stdout, prefix=prefix, lineage_file=lineage_file)
    # taxonkit echoes every input line back, so this is not expected to trigger; it catches a
    # taxid holding a newline, a blank input line, or output that stops short of the input.
    if len(rows) != len(clean_taxids):
        raise RuntimeError(
            f"taxonkit returned {len(rows)} rows for {len(clean_taxids)} taxids"
            f" of {prefix}, see {lineage_file}",
        )

    lineage_df = pd.DataFrame(rows, columns=list(TAXONKIT_OUTPUT_COLUMNS), dtype=str)

    # show unresolved counts
    unresolved = lineage_df.loc[
        lineage_df[TAXONKIT_STATUS_CODE_COLUMN].isin(_TAXONKIT_UNRESOLVED_STATUS_CODES)
    ]
    if not unresolved.empty:
        logger.warning(
            "%s: %d/%d taxids unresolved, see %s",
            prefix,
            len(unresolved),
            len(lineage_df),
            lineage_file,
        )

    return lineage_df
