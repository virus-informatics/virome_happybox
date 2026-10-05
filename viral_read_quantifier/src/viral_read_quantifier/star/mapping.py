# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

_STAR_OUTFILTER_MULTIMAP_NMAX = 1
_STAR_OUTFILTER_MISMATCH_N_OVER_READLMAX = 0.05
_STAR_OUTFILTER_MATCH_NMIN_OVER_LREAD = 0.80
_STAR_OUTFILTER_SCORE_MIN_OVER_LREAD = 0.80
_STAR_ALIGN_INTRON_MAX = 1

_STAR_LOG_INPUT_READS_RE = re.compile(
    r"Number of input reads\s*\|\s*(\d+)",
)
_STAR_BAM_SORT_RAM_RE = re.compile(r"--limitBAMsortRAM\s+(\d+)")
# 1 TB — guard against parse errors; large enough for any realistic node
_LIMIT_BAM_SORT_RAM_MAX_BYTES = 1_000_000_000_000


class StarInputReadsZeroError(RuntimeError):
    """STAR completed successfully but Log.final.out reports 0 input reads."""


def _parse_suggested_bam_sort_ram(stderr: str) -> int | None:
    """Extract ``--limitBAMsortRAM`` value from STAR stderr.

    Returns ``None`` if the pattern is not found or the value is out of range.
    """
    m = _STAR_BAM_SORT_RAM_RE.search(stderr)
    if m is None:
        logger.warning("Could not find --limitBAMsortRAM suggestion in STAR stderr")
        return None
    try:
        value = int(m.group(1))
    except (ValueError, OverflowError):
        logger.warning("Failed to parse --limitBAMsortRAM value: %s", m.group(1))
        return None
    if value <= 0 or value > _LIMIT_BAM_SORT_RAM_MAX_BYTES:
        logger.warning("Suggested --limitBAMsortRAM value out of range: %d", value)
        return None
    return value


@dataclass(frozen=True)
class StarMapResult:
    bam: Path
    log_final: Path


def parse_star_log_input_reads(log_final: Path) -> int | None:
    """Parse ``Number of input reads`` from STAR ``Log.final.out``.

    Returns the count, or ``None`` when the file is missing / unreadable /
    does not contain the expected line.
    """
    try:
        text = log_final.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None
    m = _STAR_LOG_INPUT_READS_RE.search(text)
    if m is None:
        return None
    return int(m.group(1))


def _bam_output_exists(result: StarMapResult) -> bool:
    return result.bam.is_file() and result.bam.stat().st_size > 0 and result.log_final.is_file()


def _should_skip(result: StarMapResult) -> bool:
    if not _bam_output_exists(result):
        return False
    input_reads = parse_star_log_input_reads(result.log_final)
    if input_reads is None or input_reads == 0:
        logger.warning(
            "Existing STAR output has 0 or unparseable input reads; will re-run: %s",
            result.log_final,
        )
        return False
    return True


def run_star_mapping(
    *,
    run_id: str,
    fastp_r1: Path,
    fastp_r2: Path | None,
    index_dir: Path,
    out_dir: Path,
    star_bin: Path,
    threads: int,
    force: bool,
    limit_bam_sort_ram_bytes: int | None = None,
) -> StarMapResult:
    """Run STAR mapping. Returns paths to sorted BAM and Log.final.out."""
    sample_dir = out_dir / run_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    prefix = str(sample_dir / f"{run_id}.")
    bam = Path(f"{prefix}Aligned.sortedByCoord.out.bam")
    log_final = Path(f"{prefix}Log.final.out")

    result = StarMapResult(bam=bam, log_final=log_final)

    if not force and _should_skip(result):
        logger.info("STAR output already exists, skipping: %s", run_id)
        return result

    star_legacy_tmp = sample_dir / f"{run_id}._STARtmp"
    if star_legacy_tmp.is_dir():
        shutil.rmtree(star_legacy_tmp)

    star_tmp = sample_dir / f"{run_id}.STARtmp.{os.getpid()}"
    if star_tmp.is_dir():
        shutil.rmtree(star_tmp)

    read_files = [str(fastp_r1)]
    if fastp_r2 is not None:
        read_files.append(str(fastp_r2))

    cmd = [
        str(star_bin),
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(index_dir),
        "--readFilesIn",
        *read_files,
        "--readFilesCommand",
        "zcat",
        "--outFileNamePrefix",
        prefix,
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        "--outSAMattributes",
        "NH",
        "HI",
        "AS",
        "nM",
        "NM",
        "MD",
        "--outFilterMultimapNmax",
        str(_STAR_OUTFILTER_MULTIMAP_NMAX),
        "--outFilterMismatchNoverReadLmax",
        str(_STAR_OUTFILTER_MISMATCH_N_OVER_READLMAX),
        "--outFilterMatchNminOverLread",
        str(_STAR_OUTFILTER_MATCH_NMIN_OVER_LREAD),
        "--outFilterScoreMinOverLread",
        str(_STAR_OUTFILTER_SCORE_MIN_OVER_LREAD),
        "--alignIntronMax",
        str(_STAR_ALIGN_INTRON_MAX),
        "--outTmpDir",
        str(star_tmp),
    ]

    if limit_bam_sort_ram_bytes is not None:
        cmd.extend(["--limitBAMsortRAM", str(limit_bam_sort_ram_bytes)])

    logger.info("Running STAR mapping for %s", run_id)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    if proc.returncode == 102 and limit_bam_sort_ram_bytes is None:
        suggested = _parse_suggested_bam_sort_ram(proc.stderr)
        if suggested is not None:
            logger.warning(
                "STAR BAM sorting failed due to insufficient memory. "
                "Retrying with --limitBAMsortRAM %d (suggested by STAR)",
                suggested,
            )
            if star_tmp.is_dir():
                shutil.rmtree(star_tmp)
            cmd.extend(["--limitBAMsortRAM", str(suggested)])
            subprocess.run(cmd, check=True)
        else:
            if proc.stderr:
                logger.error("STAR stderr:\n%s", proc.stderr.rstrip())
            raise subprocess.CalledProcessError(
                proc.returncode,
                cmd,
                output=proc.stdout,
                stderr=proc.stderr,
            )
    elif proc.returncode != 0:
        if proc.stderr:
            logger.error("STAR stderr:\n%s", proc.stderr.rstrip())
        raise subprocess.CalledProcessError(
            proc.returncode,
            cmd,
            output=proc.stdout,
            stderr=proc.stderr,
        )

    input_reads = parse_star_log_input_reads(result.log_final)
    if input_reads is None or input_reads == 0:
        raise StarInputReadsZeroError(
            f"STAR completed but Log.final.out reports 0 (or unparseable) input reads: "
            f"{result.log_final}",
        )

    return result


def quickcheck_bam(bam: Path) -> None:
    """Run samtools quickcheck on a BAM file. Raises on failure."""
    result = subprocess.run(
        ["samtools", "quickcheck", str(bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"samtools quickcheck failed for {bam}: {result.stderr.strip()}")


def index_bam(bam: Path, threads: int) -> None:
    """Create BAM index (.bai) if it does not exist."""
    bai = Path(f"{bam}.bai")
    if bai.is_file():
        return
    subprocess.run(
        ["samtools", "index", "-@", str(threads), str(bam)],
        check=True,
    )
