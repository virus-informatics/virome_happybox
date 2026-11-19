# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Shared pipeline-step helpers: tool checks, input validation, caching, subprocess."""

from __future__ import annotations

import logging
import subprocess
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)


_VIRUS_FAMILY_PATH_SEPARATORS = ("/", "\\")
_VIRUS_FAMILY_RESERVED_DIR_NAMES = (".", "..")


def validate_virus_family(virus_family: str) -> None:
    """Reject names that would not survive the output layout's {virus_family} directory.

    Rejecting rather than replacing keeps the directory name equal to what was requested:
    an absolute name would silently discard the parent directories, and "." would collapse
    the per-family directories onto each other.
    """
    if not virus_family or virus_family.strip() != virus_family:
        raise ValueError(f"Virus family name is empty or padded with whitespace: {virus_family!r}")
    if virus_family in _VIRUS_FAMILY_RESERVED_DIR_NAMES:
        raise ValueError(f"Virus family name must not be a relative path element: {virus_family!r}")
    for separator in _VIRUS_FAMILY_PATH_SEPARATORS:
        if separator in virus_family:
            raise ValueError(
                f"Virus family name must not contain {separator!r}: {virus_family!r}",
            )


def validate_input_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")


def _contains_fasta_header(path: Path) -> bool:
    with path.open() as handle:
        return any(line.startswith(">") for line in handle)


def validate_step_output(path: Path, label: str, *, expect_fasta: bool = False) -> None:
    """Fail at the step that produced the bad output, not at the next step that reads it."""
    if not path.exists():
        raise FileNotFoundError(f"{label} was not created: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")
    if expect_fasta and not _contains_fasta_header(path):
        raise ValueError(f"{label} has no FASTA header line: {path}")


_STEP_OUTPUT_STAGING_SUFFIX = ".tmp"


@contextmanager
def staged_output(final_output: Path) -> Iterator[Path]:
    """Yield a staging path that is renamed onto final_output only after the block succeeds.

    should_skip treats any non-empty file as a completed step, so a partial file written
    directly to the final path would be picked up as a cache hit on the next run.
    """
    staging_path = final_output.with_name(final_output.name + _STEP_OUTPUT_STAGING_SUFFIX)
    staging_path.unlink(missing_ok=True)
    try:
        yield staging_path
    except BaseException:
        staging_path.unlink(missing_ok=True)
        raise
    staging_path.replace(final_output)


def _is_complete_output(output: Path) -> bool:
    return output.exists() and output.stat().st_size > 0


def should_skip(output: Path, *, force: bool, step_name: str) -> bool:
    if _is_complete_output(output) and not force:
        logger.info("Skipping %s (cached): %s", step_name, output)
        return True
    return False


def step_workdir(output_dir: Path, command: str, subcommand: str, virus_family: str) -> Path:
    """Directory one subcommand keeps the outputs of one virus family in.

    Shared rather than rebuilt per module: a step that reads an earlier step's outputs has
    to resolve them through the same layout the earlier step writes them to.
    """
    return output_dir / command / subcommand / virus_family


_STEP_DONE_MARKER_SUFFIX = ".done"


def step_done_marker(workdir: Path, command: str, subcommand: str) -> Path:
    """Path of the marker recording that a subcommand ran to the end in this workdir.

    Named after the subcommand rather than left bare so that the file still says which step
    it belongs to once it is read outside its own directory.
    """
    return workdir / f"{command}_{subcommand}{_STEP_DONE_MARKER_SUFFIX}"


def step_is_done(marker: Path, outputs: Sequence[Path]) -> bool:
    """Report whether an earlier run finished the step and left its outputs behind.

    Outputs written directly cannot tell a finished run from one interrupted between them,
    which is what the marker is for. It is created empty, so its presence is the whole
    signal; testing it for a non-zero size would never see it. The outputs are checked
    alongside it so that losing one sends the step through again.
    """
    if not marker.exists():
        return False

    return all(_is_complete_output(output) for output in outputs)


def mark_step_done(marker: Path, outputs: Sequence[Path]) -> None:
    """Create the marker, refusing to when an output the later steps read is unusable."""
    for output in outputs:
        validate_step_output(output, output.name)

    marker.touch()


_SUBPROCESS_STDERR_LOG_TAIL_LINES = 80


def _tail_lines(text: str, max_lines: int) -> str:
    return "\n".join(text.rstrip().split("\n")[-max_lines:])


def _log_failed_subprocess(
    step_name: str,
    returncode: int,
    stderr: str,
    stderr_path: Path | None,
) -> None:
    """Report stderr on failure.

    CalledProcessError.__str__ omits the stderr attribute, so a bare traceback would
    lose the only explanation the tool produced.
    """
    logger.error("%s failed with exit status %d", step_name, returncode)
    if stderr:
        logger.error(
            "%s stderr (last %d lines):\n%s",
            step_name,
            _SUBPROCESS_STDERR_LOG_TAIL_LINES,
            _tail_lines(stderr, _SUBPROCESS_STDERR_LOG_TAIL_LINES),
        )
    if stderr_path is not None:
        logger.error("%s full stderr: %s", step_name, stderr_path)


def run_subprocess(
    cmd: list[str],
    step_name: str,
    *,
    stdin: str | None = None,
    stderr_path: Path | None = None,
) -> str:
    """Run a subprocess and return its stdout."""
    logger.info("Running %s", step_name)
    if stderr_path is not None:
        stderr_path.unlink(missing_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True, input=stdin)
    if stderr_path is not None:
        stderr_path.write_text(result.stderr)
    if result.returncode != 0:
        _log_failed_subprocess(step_name, result.returncode, result.stderr, stderr_path)
        raise subprocess.CalledProcessError(
            result.returncode,
            cmd,
            result.stdout,
            result.stderr,
        )
    if result.stderr:
        logger.debug("%s stderr:\n%s", step_name, result.stderr.rstrip())
    return result.stdout
