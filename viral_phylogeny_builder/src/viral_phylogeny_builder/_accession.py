# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Accession handling shared by the steps that join on the versionless accession."""

from __future__ import annotations


def strip_accession_version(accession: str) -> str:
    """Reduce an accession to its versionless base (`NC_004718.3` -> `NC_004718`)."""
    # Split on a space rather than on any whitespace: the latter returns no token at all for
    # an empty accession, which the row filter is meant to drop rather than fail on
    return accession.strip().split(" ")[0].split(".", 1)[0]
