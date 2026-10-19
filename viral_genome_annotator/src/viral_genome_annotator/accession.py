# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Notation of NCBI accessions shared by update_refseq and visualize."""

from __future__ import annotations


def strip_accession_version(accession: str) -> str:
    """GCF_000839185.1 -> GCF_000839185, NC_001678.1 -> NC_001678"""
    dot_idx = accession.rfind(".")
    return accession[:dot_idx] if dot_idx != -1 else accession
