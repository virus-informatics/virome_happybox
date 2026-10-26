# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""MMseqs2 wrappers for the RdRp homology search."""

from __future__ import annotations

import shutil
from pathlib import Path

from viral_phylogeny_builder._pipeline import run_subprocess

MMSEQS_REQUIRED_TOOLS = ("mmseqs",)

_MMSEQS_SEARCH_EVALUE = "1e-5"
_MMSEQS_SEARCH_MAX_SEQS = "10"

# convertalis writes no header line, so the readers of the hits TSV depend on this order.
MMSEQS_HITS_COLUMNS = (
    "query",
    "target",
    "pident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
    "tlen",
    "qlen",
    "qcov",
    "theader",
)

_MMSEQS_TARGET_DB_NAME = "mmseqs_rdrp_db"
_MMSEQS_QUERY_DB_SUFFIX = ".orfs.query_db"
_MMSEQS_RESULT_DB_SUFFIX = ".orfs.result_db"
_MMSEQS_SEARCH_TMP_DIR_SUFFIX = ".mmseqs_search_tmp"


def _run_mmseqs_createdb(input_fasta: Path, db_prefix: Path, *, mmseqs_cmd: str = "mmseqs") -> None:
    run_subprocess(
        [mmseqs_cmd, "createdb", str(input_fasta), str(db_prefix)],
        "mmseqs createdb",
    )


def _run_mmseqs_search(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    tmp_dir: Path,
    *,
    threads: int,
    mmseqs_cmd: str = "mmseqs",
) -> None:
    run_subprocess(
        [
            mmseqs_cmd,
            "search",
            str(query_db),
            str(target_db),
            str(result_db),
            str(tmp_dir),
            "-e",
            _MMSEQS_SEARCH_EVALUE,
            "--threads",
            str(threads),
            "--max-seqs",
            _MMSEQS_SEARCH_MAX_SEQS,
        ],
        "mmseqs search",
    )


def _run_mmseqs_convertalis(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    output_tsv: Path,
    *,
    mmseqs_cmd: str = "mmseqs",
) -> None:
    run_subprocess(
        [
            mmseqs_cmd,
            "convertalis",
            str(query_db),
            str(target_db),
            str(result_db),
            str(output_tsv),
            "--format-output",
            ",".join(MMSEQS_HITS_COLUMNS),
        ],
        "mmseqs convertalis",
    )


def run_mmseqs_search_to_tsv(
    query_fasta: Path,
    target_fasta: Path,
    output_tsv: Path,
    work_dir: Path,
    *,
    virus_family: str,
    threads: int,
    mmseqs_cmd: str = "mmseqs",
) -> None:
    """Search query_fasta against target_fasta and write the hits to output_tsv.

    work_dir holds the MMseqs2 databases. They are intermediates, not outputs to keep.
    virus_family only names them, following the archive pipeline.
    """
    # mmseqs exits non-zero when a database file is already in place, reporting
    # "result_db.dbtype exists already!". The directory is recreated rather than reused.
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    target_db = work_dir / _MMSEQS_TARGET_DB_NAME
    query_db = work_dir / f"{virus_family}{_MMSEQS_QUERY_DB_SUFFIX}"
    result_db = work_dir / f"{virus_family}{_MMSEQS_RESULT_DB_SUFFIX}"
    search_tmp = work_dir / f"{virus_family}{_MMSEQS_SEARCH_TMP_DIR_SUFFIX}"

    _run_mmseqs_createdb(target_fasta, target_db, mmseqs_cmd=mmseqs_cmd)
    # createindex is deliberately left out.
    # The MMseqs2 user guide recommends it only when the target database is reused for
    # several searches; each virus family searches its curated RdRp database once.
    # Measured on Bornaviridae: it cost 1.5 s and an 859 MB index file to save 1.0 s,
    # and returned identical hits.
    _run_mmseqs_createdb(query_fasta, query_db, mmseqs_cmd=mmseqs_cmd)
    _run_mmseqs_search(
        query_db,
        target_db,
        result_db,
        search_tmp,
        threads=threads,
        mmseqs_cmd=mmseqs_cmd,
    )
    _run_mmseqs_convertalis(query_db, target_db, result_db, output_tsv, mmseqs_cmd=mmseqs_cmd)
