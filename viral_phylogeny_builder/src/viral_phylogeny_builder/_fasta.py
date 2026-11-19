# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""FASTA reading and writing shared by the steps that subset sequences."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

# EMBOSS getorf and the archive pipeline both wrap at 80 columns.
_FASTA_LINE_WIDTH = 80


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs from a FASTA file without loading all into memory.

    *header* is the full text after ``>``, including the description.
    """
    header: str | None = None
    seq_parts: list[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if header is not None:
        yield header, "".join(seq_parts)


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Read a FASTA file into {record ID: (header, sequence)}."""
    return {header.split()[0]: (header, seq) for header, seq in iter_fasta(path)}


def write_fasta_record(handle: TextIO, header: str, seq: str) -> None:
    """Write a single FASTA record with 80-column wrapping."""
    handle.write(f">{header}\n")
    for i in range(0, len(seq), _FASTA_LINE_WIDTH):
        handle.write(seq[i : i + _FASTA_LINE_WIDTH] + "\n")


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    """Write (header, sequence) pairs as a wrapped FASTA file."""
    with path.open("w") as handle:
        for header, seq in records:
            write_fasta_record(handle, header, seq)
