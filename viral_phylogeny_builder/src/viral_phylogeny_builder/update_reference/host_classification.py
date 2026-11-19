# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Classify a host lineage with the bundled manual key and rank it for reference selection."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

UNKNOWN_HOST = "Unknown"

HOST_MANUAL_KEY_PATH = Path(__file__).parent / "host_manual_key.txt"

# Default separator taxonkit lineage puts between the nodes of one lineage.
_HOST_LINEAGE_NODE_SEPARATOR = ";"

_HOST_MANUAL_KEY_HOST_CLASS_COLUMN = "host_class"
_HOST_MANUAL_KEY_HOST_TAXONOMY_COLUMN = "host_taxonomy"
_HOST_MANUAL_KEY_INCLUDE_KEYWORD_COLUMN = "include_keyword"
_HOST_MANUAL_KEY_EXCLUDE_KEYWORD_COLUMN = "exclude_keyword"
_HOST_MANUAL_KEY_COLUMNS = (
    _HOST_MANUAL_KEY_HOST_CLASS_COLUMN,
    _HOST_MANUAL_KEY_HOST_TAXONOMY_COLUMN,
    _HOST_MANUAL_KEY_INCLUDE_KEYWORD_COLUMN,
    _HOST_MANUAL_KEY_EXCLUDE_KEYWORD_COLUMN,
)

# Would rather be built from host_manual_key.txt, but that file carries no priority column,
# and which host outranks which is a call for whoever supplies the key, not for this code.
# Host classes that rank a reference candidate, best first. Candidates whose host class is
# absent here fall through to _HOST_PRIORITY_BY_LINEAGE_KEYWORD.
_HOST_PRIORITY_BY_HOST_CLASS = {
    "Human": 0,
    "Mammalia": 1,
    "Aves": 2,
    "Arthropoda": 3,
}

# Keywords matched against the lineage nodes when the host class carries no priority of its
# own, first match wins. Acari and Flies reach a priority this way: both are arthropods, but
# the manual key gives them their own host class, absent from _HOST_PRIORITY_BY_HOST_CLASS.
_HOST_PRIORITY_BY_LINEAGE_KEYWORD_LOWER = (
    ("homo sapiens", 0),
    ("mammalia", 1),
    ("aves", 2),
    ("arthropoda", 3),
)


@dataclass(frozen=True, slots=True)
class HostManualKeyRule:
    """One row of host_manual_key.txt, with the keywords as the file spells them."""

    host_class: str
    host_taxonomy: str
    include_keyword: str
    exclude_keyword: str

    @property
    def include_keyword_lower(self) -> str:
        """The include keyword folded for the case-insensitive lineage match."""
        return self.include_keyword.lower()

    @property
    def exclude_keyword_lower(self) -> str:
        """The exclude keyword folded for the case-insensitive lineage match."""
        return self.exclude_keyword.lower()


@dataclass(frozen=True, slots=True)
class HostClassification:
    """Manual host class / taxonomy of one lineage and the priority derived from them."""

    manual_host_class: str
    manual_host_taxonomy: str
    host_priority: int


def read_host_manual_key(path: Path = HOST_MANUAL_KEY_PATH) -> tuple[HostManualKeyRule, ...]:
    """Read the manual host key TSV, keeping row order: the first matching rule wins.

    The path defaults to the copy bundled with the package; the parameter is for tests
    only, as the classification rules are not meant to be swapped per run.

    Raises ValueError when the header lacks a required column or no usable rule is left.
    """
    # keep_default_na leaves a keyword such as "NA" or "None" as text: the defaults would read
    # it as a missing value and the rule would be dropped as if it had no include keyword.
    manual_key = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    missing_columns = [
        column for column in _HOST_MANUAL_KEY_COLUMNS if column not in manual_key.columns
    ]
    if missing_columns:
        raise ValueError(
            f"Host manual key is missing columns: {missing_columns} (in {path})",
        )

    rules = []
    for _, row in manual_key.iterrows():
        include_keyword = row[_HOST_MANUAL_KEY_INCLUDE_KEYWORD_COLUMN].strip()
        # A rule without an include keyword would match every lineage
        if not include_keyword:
            continue
        rules.append(
            HostManualKeyRule(
                host_class=row[_HOST_MANUAL_KEY_HOST_CLASS_COLUMN].strip(),
                host_taxonomy=row[_HOST_MANUAL_KEY_HOST_TAXONOMY_COLUMN].strip(),
                include_keyword=include_keyword,
                exclude_keyword=row[_HOST_MANUAL_KEY_EXCLUDE_KEYWORD_COLUMN].strip(),
            ),
        )

    if not rules:
        raise ValueError(f"Host manual key holds no classification rule: {path}")

    return tuple(rules)


# Matched against the lineage as one string, a keyword lands inside a longer name:
# Homo inside Acanthomorphata, Aves inside the epithet of Ochlerotatus flavescens.
def _host_lineage_nodes_lower(host_lineage: str) -> frozenset[str]:
    return frozenset(
        stripped.lower()
        for node in host_lineage.split(_HOST_LINEAGE_NODE_SEPARATOR)
        if (stripped := node.strip())
    )


def _host_priority(host_class: str, host_lineage_nodes_lower: frozenset[str]) -> int:
    priority_of_host_class = _HOST_PRIORITY_BY_HOST_CLASS.get(host_class)
    if priority_of_host_class is not None:
        return priority_of_host_class

    for keyword_lower, priority_of_keyword in _HOST_PRIORITY_BY_LINEAGE_KEYWORD_LOWER:
        if keyword_lower in host_lineage_nodes_lower:
            return priority_of_keyword

    return len(_HOST_PRIORITY_BY_HOST_CLASS)


def classify_host_lineage(
    host_lineage: str,
    rules: tuple[HostManualKeyRule, ...],
) -> HostClassification:
    """Match one host lineage against the rules and score the host for reference selection.

    A lineage that matches no rule is classified as UNKNOWN_HOST; it still gets a priority,
    since the lineage keywords are matched independently of the manual host class.
    """
    host_lineage_nodes_lower = _host_lineage_nodes_lower(host_lineage)

    manual_host_class = UNKNOWN_HOST
    manual_host_taxonomy = UNKNOWN_HOST
    for rule in rules:
        if rule.include_keyword_lower not in host_lineage_nodes_lower:
            continue
        if rule.exclude_keyword_lower and rule.exclude_keyword_lower in host_lineage_nodes_lower:
            continue
        manual_host_class = rule.host_class
        manual_host_taxonomy = rule.host_taxonomy
        break

    return HostClassification(
        manual_host_class=manual_host_class,
        manual_host_taxonomy=manual_host_taxonomy,
        host_priority=_host_priority(manual_host_class, host_lineage_nodes_lower),
    )
