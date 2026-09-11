# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import json
import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

_FASTP_QUALITY_PHRED = 20
_FASTP_MIN_READ_LENGTH = 30


@dataclass(frozen=True)
class FastpResult:
    trimmed_r1: Path
    trimmed_r2: Path | None
    json_report: Path


def _outputs_exist(result: FastpResult) -> bool:
    if not result.trimmed_r1.is_file() or result.trimmed_r1.stat().st_size == 0:
        return False
    if result.trimmed_r2 is not None:
        if not result.trimmed_r2.is_file() or result.trimmed_r2.stat().st_size == 0:
            return False
    return True


def run_fastp(
    *,
    run_id: str,
    raw_r1: Path,
    raw_r2: Path | None,
    out_dir: Path,
    threads: int,
    force: bool,
) -> FastpResult:
    """Run fastp QC. Returns paths to trimmed FASTQ and JSON report."""
    sample_dir = out_dir / run_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    json_report = sample_dir / f"{run_id}.fastp.json"
    html_report = sample_dir / f"{run_id}.fastp.html"

    if raw_r2 is not None:
        trimmed_r1 = sample_dir / f"{run_id}_R1.fastp.fastq.gz"
        trimmed_r2 = sample_dir / f"{run_id}_R2.fastp.fastq.gz"
    else:
        trimmed_r1 = sample_dir / f"{run_id}.fastp.fastq.gz"
        trimmed_r2 = None

    result = FastpResult(trimmed_r1=trimmed_r1, trimmed_r2=trimmed_r2, json_report=json_report)

    if not force and _outputs_exist(result):
        logger.info("fastp output already exists, skipping: %s", run_id)
        return result

    cmd = ["fastp", "-i", str(raw_r1)]
    if raw_r2 is not None:
        cmd.extend(["-I", str(raw_r2), "-o", str(trimmed_r1), "-O", str(trimmed_r2)])
        cmd.append("--detect_adapter_for_pe")
    else:
        cmd.extend(["-o", str(trimmed_r1)])

    cmd.extend(
        [
            "--qualified_quality_phred",
            str(_FASTP_QUALITY_PHRED),
            "--length_required",
            str(_FASTP_MIN_READ_LENGTH),
            "-y",
            "-3",
            "-x",
            "--thread",
            str(threads),
            "--html",
            str(html_report),
            "--json",
            str(json_report),
        ],
    )

    logger.info("Running fastp for %s", run_id)
    subprocess.run(cmd, check=True)

    return result


def read_fastp_after_filtering_total_reads(json_report: Path) -> int:
    """Read total_reads (individual read count) from fastp after_filtering summary."""
    with json_report.open(encoding="utf-8") as f:
        data = json.load(f)

    total = data["summary"]["after_filtering"]["total_reads"]
    if not isinstance(total, int) or total <= 0:
        raise ValueError(
            f"Invalid after_filtering.total_reads in fastp JSON: {total} ({json_report})",
        )
    return total
