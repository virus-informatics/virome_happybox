# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Query sequence RdRp search (getorf + BLASTP)."""

from __future__ import annotations

import csv
import logging
import re
import tempfile
from pathlib import Path

from viral_phylogeny_builder._blastp import (
    BLASTP_OUTFMT_COLUMNS,
    BLASTP_REQUIRED_TOOLS,
    run_blastp,
    run_makeblastdb,
    sanitize_evalue_tag,
)
from viral_phylogeny_builder._commands import QUERY_COMMAND, QUERY_SEARCH_SUBCOMMAND
from viral_phylogeny_builder._fasta import iter_fasta, write_fasta_record
from viral_phylogeny_builder._orf import ORF_REQUIRED_TOOLS, run_getorf
from viral_phylogeny_builder._pipeline import (
    run_subprocess,
    should_skip,
    staged_output,
    step_workdir,
    validate_input_file,
    validate_step_output,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools

logger = logging.getLogger(__name__)

QUERY_SEARCH_INPUT_ORFS_FILENAME = "input.orf.faa"
_BLASTP_DB_DIR_NAME = "blastp_db"
_BLASTP_DB_PREFIX_NAME = "rdrp_db"

_EVALUE_COLUMN_INDEX = BLASTP_OUTFMT_COLUMNS.index("evalue")
_QSEQID_COLUMN_INDEX = BLASTP_OUTFMT_COLUMNS.index("qseqid")

_HIT_EXTRACTION_REQUIRED_TOOLS = ("seqkit",)

_NAME_MAP_COLUMNS = ("node", "contig", "orf_id")


def _infer_contig_id(orf_header: str) -> str:
    """Infer the source contig ID from a getorf ORF header.

    getorf names ORFs as ``{contig}_{N}`` where N is a 1-based ORF index.
    Coordinate suffixes like ``:start-end`` are also stripped.
    """
    name = orf_header.strip().split()[0]
    name = name.split(":")[0]
    name = re.sub(r"_[0-9]+$", "", name)
    name = re.sub(r"\.[0-9]+$", "", name)
    return name


def _renumber_by_contig(
    hit_fasta: Path,
    output_fasta: Path,
    output_name_map: Path,
) -> None:
    """Renumber hit ORFs per contig and write the renamed FASTA + name map."""
    contig_counts: dict[str, int] = {}

    with output_fasta.open("w") as fout, output_name_map.open("w", newline="") as fmap:
        writer = csv.writer(fmap, delimiter="\t")
        writer.writerow(_NAME_MAP_COLUMNS)

        for header, seq in iter_fasta(hit_fasta):
            orf_id = header.split()[0]
            contig = _infer_contig_id(orf_id)
            contig_counts[contig] = contig_counts.get(contig, 0) + 1

            if contig_counts[contig] == 1:
                node = contig
            else:
                node = f"{contig}_orf{contig_counts[contig]}"

            write_fasta_record(fout, node, seq)
            writer.writerow((node, contig, orf_id))


def _filter_hit_ids(blastp_tsv: Path, evalue_threshold: float) -> set[str]:
    """Extract unique qseqid values that pass the evalue threshold."""
    kept: set[str] = set()
    with blastp_tsv.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            evalue = float(row[_EVALUE_COLUMN_INDEX])
            if evalue <= evalue_threshold:
                kept.add(row[_QSEQID_COLUMN_INDEX])
    return kept


def _extract_hit_fasta(
    hit_ids: set[str],
    orf_fasta: Path,
    output_fasta: Path,
    *,
    seqkit_cmd: str = "seqkit",
) -> None:
    """Extract ORF sequences matching hit_ids via seqkit grep."""
    ids_file = Path(tempfile.mkstemp(suffix=".txt", dir=output_fasta.parent)[1])
    try:
        ids_file.write_text("\n".join(sorted(hit_ids)) + "\n")
        run_subprocess(
            [
                seqkit_cmd,
                "grep",
                "-f",
                str(ids_file),
                str(orf_fasta),
                "-o",
                str(output_fasta),
            ],
            "seqkit grep (hit extraction)",
        )
    finally:
        ids_file.unlink(missing_ok=True)


def run_query_search(
    virus_family: str,
    input_fasta: Path,
    rdrp_db: Path,
    output_dir: Path,
    *,
    blastp_evalue: float = 1e-5,
    threads: int = 4,
    force: bool = False,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Run query RdRp search and return the path to the hit ORF FASTA.

    Renaming will follow.
    """
    resolved = resolve_tools(
        ORF_REQUIRED_TOOLS + BLASTP_REQUIRED_TOOLS + _HIT_EXTRACTION_REQUIRED_TOOLS,
        tool_config,
    )
    validate_virus_family(virus_family)
    validate_input_file(input_fasta, "Query nucleotide FASTA")
    validate_input_file(rdrp_db, "RdRp protein DB FASTA")

    workdir = step_workdir(output_dir, QUERY_COMMAND, QUERY_SEARCH_SUBCOMMAND, virus_family)
    workdir.mkdir(parents=True, exist_ok=True)

    orfs_output = workdir / QUERY_SEARCH_INPUT_ORFS_FILENAME
    if not should_skip(orfs_output, force=force, step_name="getorf"):
        with staged_output(orfs_output) as staging_path:
            run_getorf(input_fasta, staging_path, getorf_cmd=resolved["getorf"])
            validate_step_output(staging_path, "getorf output", expect_fasta=True)

    # ── makeblastdb ──────────────────────────────────────────────
    db_dir = workdir / _BLASTP_DB_DIR_NAME
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / _BLASTP_DB_PREFIX_NAME

    if not should_skip(db_prefix.with_suffix(".pin"), force=force, step_name="makeblastdb"):
        run_makeblastdb(rdrp_db, db_prefix, makeblastdb_cmd=resolved["makeblastdb"])

    # ── blastp ────────────────────────────────────────────────────
    evalue_tag = sanitize_evalue_tag(blastp_evalue)
    blastp_tsv = workdir / f"{virus_family}.evalue_{evalue_tag}.blastp.tsv"

    if not should_skip(blastp_tsv, force=force, step_name="blastp"):
        with staged_output(blastp_tsv) as staging_path:
            run_blastp(
                orfs_output,
                db_prefix,
                staging_path,
                evalue=blastp_evalue,
                threads=threads,
                blastp_cmd=resolved["blastp"],
            )
            validate_step_output(staging_path, "blastp output")

    # ── hit extraction ───────────────────────────────────────────
    hit_fasta = workdir / f"{virus_family}.evalue_{evalue_tag}.fasta"

    if not should_skip(hit_fasta, force=force, step_name="hit extraction"):
        hit_ids = _filter_hit_ids(blastp_tsv, blastp_evalue)
        if not hit_ids:
            raise ValueError(
                f"No BLASTP hits passed evalue threshold {blastp_evalue}. "
                "Consider increasing --blastp-evalue.",
            )
        logger.info("BLASTP hit ORFs after evalue filter: %d", len(hit_ids))

        with staged_output(hit_fasta) as staging_path:
            _extract_hit_fasta(hit_ids, orfs_output, staging_path, seqkit_cmd=resolved["seqkit"])
            validate_step_output(staging_path, "hit FASTA", expect_fasta=True)

    # ── contig renumbering ──────────────────────────────────────
    renamed_fasta = workdir / f"{virus_family}.evalue_{evalue_tag}.renamed.fasta"
    name_map_tsv = workdir / f"{virus_family}.evalue_{evalue_tag}.name_map.tsv"

    if not should_skip(renamed_fasta, force=force, step_name="contig renumbering"):
        with (
            staged_output(renamed_fasta) as staging_fasta,
            staged_output(
                name_map_tsv,
            ) as staging_map,
        ):
            _renumber_by_contig(hit_fasta, staging_fasta, staging_map)
            validate_step_output(staging_fasta, "renamed FASTA", expect_fasta=True)
            validate_step_output(staging_map, "name map TSV")

    logger.info("Query search complete: %s", renamed_fasta)
    return renamed_fasta
