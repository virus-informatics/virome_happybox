# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Extract RdRp sequences from the reference genomes (getorf -> MMseqs2 -> best hit)."""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from viral_phylogeny_builder._accession import strip_accession_version
from viral_phylogeny_builder._commands import (
    BUILD_REFTREE_COMMAND,
    BUILD_REFTREE_SEARCH_SUBCOMMAND,
)
from viral_phylogeny_builder._fasta import read_fasta, write_fasta
from viral_phylogeny_builder._mmseqs import (
    MMSEQS_HITS_COLUMNS,
    MMSEQS_REQUIRED_TOOLS,
    run_mmseqs_search_to_tsv,
)
from viral_phylogeny_builder._orf import ORF_REQUIRED_TOOLS, run_getorf
from viral_phylogeny_builder._pipeline import (
    mark_step_done,
    should_skip,
    staged_output,
    step_done_marker,
    step_is_done,
    step_workdir,
    validate_input_file,
    validate_step_output,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools
from viral_phylogeny_builder._virus_genome_metadata_columns import (
    METADATA_ACCESSION_COLUMN,
)

logger = logging.getLogger(__name__)

BUILD_REFTREE_SEARCH_ORFS_FILENAME_SUFFIX = ".orfs.faa"
BUILD_REFTREE_SEARCH_HITS_FILENAME_SUFFIX = ".rdrp.mmseqs.tsv"
BUILD_REFTREE_SEARCH_RDRP_FILENAME_SUFFIX = ".rdrp.faa"
BUILD_REFTREE_SEARCH_RDRP_METADATA_FILENAME_SUFFIX = ".rdrp_metadata.csv"

# Holds the MMseqs2 databases, which are intermediates of the search rather than outputs.
_MMSEQS_WORK_DIR_NAME = "tmp_mmseqs"

# The values the archive pipeline passed on the command line.
_MIN_QCOV = 0.1
_MIN_ALNLEN = 100

_HITS_NUMERIC_COLUMNS = ["pident", "alnlen", "evalue", "bits", "qcov", "qlen"]


def build_reftree_search_orfs_filename(virus_family: str) -> str:
    """Name of the ORF FASTA search writes for one virus family.

    The name repeats the virus family the workdir already carries so that the file stays
    identifiable once it is copied out of the output layout.
    """
    return f"{virus_family}{BUILD_REFTREE_SEARCH_ORFS_FILENAME_SUFFIX}"


def build_reftree_search_hits_filename(virus_family: str) -> str:
    """Name of the MMseqs2 hits TSV search writes for one virus family."""
    return f"{virus_family}{BUILD_REFTREE_SEARCH_HITS_FILENAME_SUFFIX}"


def build_reftree_search_rdrp_filename(virus_family: str) -> str:
    """Name of the RdRp FASTA search writes for one virus family."""
    return f"{virus_family}{BUILD_REFTREE_SEARCH_RDRP_FILENAME_SUFFIX}"


def build_reftree_search_rdrp_metadata_filename(virus_family: str) -> str:
    """Name of the RdRp metadata CSV search writes for one virus family."""
    return f"{virus_family}{BUILD_REFTREE_SEARCH_RDRP_METADATA_FILENAME_SUFFIX}"


def accession_from_getorf_query(query: str) -> str:
    """Recover the genome record ID that getorf built one ORF ID from."""
    q = str(query).strip()

    # EMBOSS getorf headers often preserve the original sequence ID at the start.
    # Examples:
    # NC_045512.2_1
    # NC_045512.2:1-1000
    # lcl|NC_045512.2_cds_...
    m = re.search(r"([A-Z]{1,4}_[0-9]+\.[0-9]+|[A-Z]{1,4}[0-9]+\.[0-9]+)", q)
    if m:
        return m.group(1)

    base = q.split()[0].split(":")[0].split("_")[0]
    return base


def extract_best_rdrp_hits(
    orf_fasta: Path,
    hits_tsv: Path,
    metadata_csv: Path,
    out_fasta: Path,
    out_metadata: Path,
) -> None:
    """Pick one best RdRp ORF per genome and write it with the joined reference metadata."""
    hits = pd.read_csv(hits_tsv, sep="\t", header=None, dtype=str)
    hits.columns = list(MMSEQS_HITS_COLUMNS)

    for col in _HITS_NUMERIC_COLUMNS:
        hits[col] = pd.to_numeric(hits[col], errors="coerce")

    hits = hits.loc[hits["alnlen"].fillna(0) >= _MIN_ALNLEN].copy()
    hits = hits.loc[hits["qcov"].fillna(0) >= _MIN_QCOV].copy()

    if hits.empty:
        raise ValueError("No RdRp hits remained after filtering.")

    hits["accession"] = hits["query"].map(accession_from_getorf_query)
    hits["accession_base"] = hits["accession"].map(strip_accession_version)

    hits = hits.sort_values(
        by=["accession_base", "evalue", "bits", "alnlen", "qcov"],
        ascending=[True, True, False, False, False],
    )

    best = hits.groupby("accession_base", as_index=False).head(1).copy()

    orf_records = read_fasta(orf_fasta)

    fasta_out: list[tuple[str, str]] = []
    for _, row in best.iterrows():
        query = str(row["query"])
        if query not in orf_records:
            continue

        _, seq = orf_records[query]
        node = strip_accession_version(row["accession"])
        fasta_out.append((node, seq))

    if not fasta_out:
        raise ValueError("No ORF sequences matched best-hit query IDs.")

    with staged_output(out_fasta) as staging_path:
        write_fasta(fasta_out, staging_path)
        validate_step_output(staging_path, "RdRp FASTA", expect_fasta=True)

    meta = pd.read_csv(metadata_csv, dtype=str).fillna("")
    meta["accession_base"] = meta[METADATA_ACCESSION_COLUMN].map(strip_accession_version)

    best_meta = best.merge(meta, on="accession_base", how="left", suffixes=("_hit", ""))

    with staged_output(out_metadata) as staging_path:
        best_meta.to_csv(staging_path, index=False)
        validate_step_output(staging_path, "RdRp metadata CSV")

    logger.info("Best RdRp hits: %d", len(best))
    logger.info("RdRp FASTA records: %d", len(fasta_out))


def run_build_reftree_search(
    virus_family: str,
    genome_fasta: Path,
    reference_metadata_csv: Path,
    rdrp_db: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    force: bool = False,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Run build-reftree search and return the path to the extracted RdRp FASTA."""
    resolved = resolve_tools(ORF_REQUIRED_TOOLS + MMSEQS_REQUIRED_TOOLS, tool_config)
    validate_virus_family(virus_family)
    validate_input_file(genome_fasta, "Reference genome FASTA")
    # The best-hit extraction is what reads this one.
    # Checked up front so that a missing file is reported before getorf runs on the genomes.
    validate_input_file(reference_metadata_csv, "Reference metadata CSV")
    validate_input_file(rdrp_db, "Curated RdRp FASTA")

    workdir = step_workdir(
        output_dir,
        BUILD_REFTREE_COMMAND,
        BUILD_REFTREE_SEARCH_SUBCOMMAND,
        virus_family,
    )
    workdir.mkdir(parents=True, exist_ok=True)

    orfs_output = workdir / build_reftree_search_orfs_filename(virus_family)
    hits_output = workdir / build_reftree_search_hits_filename(virus_family)
    rdrp_output = workdir / build_reftree_search_rdrp_filename(virus_family)
    rdrp_metadata_output = workdir / build_reftree_search_rdrp_metadata_filename(virus_family)
    step_outputs = [orfs_output, hits_output, rdrp_output, rdrp_metadata_output]

    marker = step_done_marker(workdir, BUILD_REFTREE_COMMAND, BUILD_REFTREE_SEARCH_SUBCOMMAND)
    if not force and step_is_done(marker, step_outputs):
        logger.info("Skipping %s (cached): %s", BUILD_REFTREE_SEARCH_SUBCOMMAND, workdir)
        return rdrp_output

    marker.unlink(missing_ok=True)
    if force:
        for output in step_outputs:
            output.unlink(missing_ok=True)

    if not should_skip(orfs_output, force=force, step_name="getorf"):
        with staged_output(orfs_output) as staging_path:
            run_getorf(genome_fasta, staging_path, getorf_cmd=resolved["getorf"])
            validate_step_output(staging_path, "getorf output", expect_fasta=True)

    if not should_skip(hits_output, force=force, step_name="mmseqs search"):
        with staged_output(hits_output) as staging_path:
            run_mmseqs_search_to_tsv(
                orfs_output,
                rdrp_db,
                staging_path,
                workdir / _MMSEQS_WORK_DIR_NAME,
                virus_family=virus_family,
                threads=threads,
                mmseqs_cmd=resolved["mmseqs"],
            )
            validate_step_output(staging_path, "mmseqs search output")

    logger.info("Extracting best RdRp ORFs")
    extract_best_rdrp_hits(
        orfs_output,
        hits_output,
        reference_metadata_csv,
        rdrp_output,
        rdrp_metadata_output,
    )

    mark_step_done(marker, step_outputs)
    logger.info("Build-reftree search complete: %s", rdrp_output)
    logger.info("RdRp metadata: %s", rdrp_metadata_output)
    return rdrp_output
