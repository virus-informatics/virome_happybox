# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Add query RdRp sequences to reference MSA (mafft --add -> trimal -> gap filter)."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

from viral_phylogeny_builder._alignment import (
    ALIGNMENT_REQUIRED_TOOLS,
    run_seqkit_rmdup,
    run_trimal,
)
from viral_phylogeny_builder._commands import QUERY_ALIGN_SUBCOMMAND, QUERY_COMMAND
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

_MAFFT_ADDED_FILENAME = "mafft_added.fasta"
_TRIMAL_FILENAME = "trimal.fasta"
_RMDUP_FILENAME = "rmdup.fasta"
QUERY_ALIGN_FINAL_FILENAME = "final.fasta"


def _run_mafft_add(
    query_fasta: Path,
    reference_alignment: Path,
    output_fasta: Path,
    *,
    threads: int,
    mafft_cmd: str = "mafft",
) -> None:
    stdout = run_subprocess(
        [
            mafft_cmd,
            "--add",
            str(query_fasta),
            "--keeplength",
            "--thread",
            str(threads),
            str(reference_alignment),
        ],
        "mafft --add",
    )
    output_fasta.write_text(stdout)


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header_line, sequence) tuples from a FASTA file."""
    sequences: list[tuple[str, str]] = []
    current_header: str | None = None
    current_seq: list[str] = []

    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, "".join(current_seq)))

    return sequences


def _get_alignment_length(fasta_path: Path, *, seqkit_cmd: str = "seqkit") -> tuple[int, int]:
    """Return (alignment_length, sequence_count) via seqkit fx2tab."""
    stdout = run_subprocess(
        [seqkit_cmd, "fx2tab", "-n", "-l", str(fasta_path)],
        "seqkit fx2tab",
    )
    lines = [line for line in stdout.strip().split("\n") if line]
    if not lines:
        raise ValueError(f"No sequences found in {fasta_path}")
    alignment_length = int(lines[0].split("\t")[1])
    return alignment_length, len(lines)


def _filter_gappy_sequences(
    input_fasta: Path,
    output_fasta: Path,
    *,
    gap_threshold_pct: int,
    seqkit_cmd: str = "seqkit",
) -> int:
    """Remove sequences exceeding the gap threshold. Returns retained count."""
    alignment_length, total_count = _get_alignment_length(input_fasta, seqkit_cmd=seqkit_cmd)
    min_nongap_length = int(alignment_length - (alignment_length * gap_threshold_pct / 100))

    logger.info(
        "Alignment length: %d, min non-gap residues: %d (gap threshold: %d%%)",
        alignment_length,
        min_nongap_length,
        gap_threshold_pct,
    )

    # seqkit seq -g -m: remove gaps then filter by min length
    filtered_stdout = run_subprocess(
        [seqkit_cmd, "seq", "-g", "-m", str(min_nongap_length), str(input_fasta)],
        "seqkit seq -g -m (gap filter)",
    )

    # seqkit seq -n: extract valid sequence names
    valid_ids_text = run_subprocess(
        [seqkit_cmd, "seq", "-n"],
        "seqkit seq -n (gap filter)",
        stdin=filtered_stdout,
    ).strip()
    if not valid_ids_text:
        raise ValueError(
            f"All sequences filtered at {gap_threshold_pct}% gap threshold. "
            "Consider increasing --gap-threshold-pct.",
        )

    retained = len(valid_ids_text.split("\n"))

    # seqkit grep -f: extract original sequences (with gaps) by valid IDs
    valid_ids_file = Path(tempfile.mkstemp(suffix=".txt", dir=output_fasta.parent)[1])
    try:
        valid_ids_file.write_text(valid_ids_text + "\n")
        run_subprocess(
            [
                seqkit_cmd,
                "grep",
                "-f",
                str(valid_ids_file),
                str(input_fasta),
                "-o",
                str(output_fasta),
            ],
            "seqkit grep (gap filter)",
        )
    finally:
        valid_ids_file.unlink(missing_ok=True)

    logger.info("Retained %d / %d sequences", retained, total_count)
    return retained


def run_query_align(
    virus_family: str,
    query_rdrp_fasta: Path,
    reference_alignment: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    force: bool = False,
    gap_threshold_pct: int = 70,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Run query alignment pipeline and return path to final alignment FASTA."""
    validate_virus_family(virus_family)
    resolved = resolve_tools(ALIGNMENT_REQUIRED_TOOLS, tool_config)
    validate_input_file(query_rdrp_fasta, "Query RdRp FASTA")
    validate_input_file(reference_alignment, "Reference alignment")

    workdir = step_workdir(output_dir, QUERY_COMMAND, QUERY_ALIGN_SUBCOMMAND, virus_family)
    workdir.mkdir(parents=True, exist_ok=True)

    mafft_output = workdir / _MAFFT_ADDED_FILENAME
    if not should_skip(mafft_output, force=force, step_name="mafft --add"):
        with staged_output(mafft_output) as staging_path:
            _run_mafft_add(
                query_rdrp_fasta,
                reference_alignment,
                staging_path,
                threads=threads,
                mafft_cmd=resolved["mafft"],
            )
            validate_step_output(staging_path, "mafft --add output", expect_fasta=True)

    trimal_output = workdir / _TRIMAL_FILENAME
    if not should_skip(trimal_output, force=force, step_name="trimal"):
        with staged_output(trimal_output) as staging_path:
            run_trimal(mafft_output, staging_path, trimal_cmd=resolved["trimal"])
            validate_step_output(staging_path, "trimal output", expect_fasta=True)

    rmdup_output = workdir / _RMDUP_FILENAME
    if not should_skip(rmdup_output, force=force, step_name="seqkit rmdup"):
        with staged_output(rmdup_output) as staging_path:
            run_seqkit_rmdup(trimal_output, staging_path, seqkit_cmd=resolved["seqkit"])
            validate_step_output(staging_path, "seqkit rmdup output", expect_fasta=True)

    final_output = workdir / QUERY_ALIGN_FINAL_FILENAME
    if not should_skip(final_output, force=force, step_name="gap filter"):
        with staged_output(final_output) as staging_path:
            retained = _filter_gappy_sequences(
                rmdup_output,
                staging_path,
                gap_threshold_pct=gap_threshold_pct,
                seqkit_cmd=resolved["seqkit"],
            )
            if retained == 0:
                raise ValueError(
                    f"All sequences filtered at {gap_threshold_pct}% gap threshold. "
                    "Consider increasing --gap-threshold-pct.",
                )
            validate_step_output(staging_path, "gap filter output", expect_fasta=True)

    logger.info("Query alignment complete: %s", final_output)
    return final_output
