# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Build reference MSA from RdRp sequences (mafft --auto -> trimal -> seqkit rmdup)."""

from __future__ import annotations

import logging
from pathlib import Path

from viral_phylogeny_builder._alignment import (
    ALIGNMENT_REQUIRED_TOOLS,
    run_seqkit_rmdup,
    run_trimal,
)
from viral_phylogeny_builder._commands import (
    BUILD_REFTREE_ALIGN_SUBCOMMAND,
    BUILD_REFTREE_COMMAND,
)
from viral_phylogeny_builder._pipeline import (
    mark_step_done,
    run_subprocess,
    should_skip,
    staged_output,
    step_done_marker,
    step_is_done,
    validate_input_file,
    validate_step_output,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools

logger = logging.getLogger(__name__)

BUILD_REFTREE_ALIGN_MAFFT_FILENAME_SUFFIX = ".mafft.faa"
BUILD_REFTREE_ALIGN_TRIMAL_FILENAME_SUFFIX = ".mafft.trimal.faa"
BUILD_REFTREE_ALIGN_FINAL_FILENAME_SUFFIX = ".final.faa"

_MAFFT_STDERR_SUFFIX = ".stderr"


def build_reftree_align_mafft_filename(virus_family: str) -> str:
    """Name of the MAFFT alignment align writes for one virus family.

    The name repeats the virus family the workdir already carries so that the file stays
    identifiable once it is copied out of the output layout.
    """
    return f"{virus_family}{BUILD_REFTREE_ALIGN_MAFFT_FILENAME_SUFFIX}"


def build_reftree_align_trimal_filename(virus_family: str) -> str:
    """Name of the trimAl alignment align writes for one virus family."""
    return f"{virus_family}{BUILD_REFTREE_ALIGN_TRIMAL_FILENAME_SUFFIX}"


def build_reftree_align_final_filename(virus_family: str) -> str:
    """Name of the final alignment align writes for one virus family."""
    return f"{virus_family}{BUILD_REFTREE_ALIGN_FINAL_FILENAME_SUFFIX}"


def _run_mafft_auto(
    input_fasta: Path,
    output_fasta: Path,
    *,
    threads: int,
    stderr_path: Path,
    mafft_cmd: str = "mafft",
) -> None:
    stdout = run_subprocess(
        [
            mafft_cmd,
            "--auto",
            "--thread",
            str(threads),
            str(input_fasta),
        ],
        "mafft --auto",
        stderr_path=stderr_path,
    )
    output_fasta.write_text(stdout)


def run_build_reftree_align(
    virus_family: str,
    rdrp_fasta: Path,
    rdrp_metadata_csv: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    force: bool = False,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Run build-reftree alignment pipeline and return path to final alignment FASTA."""
    resolved = resolve_tools(ALIGNMENT_REQUIRED_TOOLS, tool_config)
    validate_virus_family(virus_family)
    validate_input_file(rdrp_fasta, "RdRp FASTA")
    validate_input_file(rdrp_metadata_csv, "RdRp metadata CSV")

    workdir = output_dir / BUILD_REFTREE_COMMAND / BUILD_REFTREE_ALIGN_SUBCOMMAND / virus_family
    workdir.mkdir(parents=True, exist_ok=True)

    marker = step_done_marker(workdir, BUILD_REFTREE_COMMAND, BUILD_REFTREE_ALIGN_SUBCOMMAND)
    mafft_output = workdir / build_reftree_align_mafft_filename(virus_family)
    trimal_output = workdir / build_reftree_align_trimal_filename(virus_family)
    final_output = workdir / build_reftree_align_final_filename(virus_family)
    step_outputs = [mafft_output, final_output]
    all_cached_outputs = [mafft_output, trimal_output, final_output]

    if not force and step_is_done(marker, step_outputs):
        logger.info("Skipping %s (cached): %s", BUILD_REFTREE_ALIGN_SUBCOMMAND, workdir)
        return final_output

    marker.unlink(missing_ok=True)
    if force:
        for output in all_cached_outputs:
            output.unlink(missing_ok=True)

    if not should_skip(mafft_output, force=force, step_name="mafft --auto"):
        # kept next to the final output, not the staging file, so it survives a failed run
        mafft_stderr = mafft_output.with_name(mafft_output.name + _MAFFT_STDERR_SUFFIX)
        with staged_output(mafft_output) as staging_path:
            _run_mafft_auto(
                rdrp_fasta,
                staging_path,
                threads=threads,
                stderr_path=mafft_stderr,
                mafft_cmd=resolved["mafft"],
            )
            validate_step_output(staging_path, "mafft --auto output", expect_fasta=True)

    if not should_skip(trimal_output, force=force, step_name="trimal"):
        with staged_output(trimal_output) as staging_path:
            run_trimal(mafft_output, staging_path, trimal_cmd=resolved["trimal"])
            validate_step_output(staging_path, "trimal output", expect_fasta=True)

    if not should_skip(final_output, force=force, step_name="seqkit rmdup"):
        with staged_output(final_output) as staging_path:
            run_seqkit_rmdup(trimal_output, staging_path, seqkit_cmd=resolved["seqkit"])
            validate_step_output(staging_path, "seqkit rmdup output", expect_fasta=True)

    mark_step_done(marker, step_outputs)
    logger.info("Build-reftree alignment complete: %s", final_output)
    return final_output
