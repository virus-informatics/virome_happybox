# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

# SN keys of `samtools stats`. Definitions are quoted from samtools-stats.1 (1.23.1).
_SN_KEY_RAW_TOTAL_SEQUENCES = "raw total sequences"
_SN_KEY_NON_PRIMARY_ALIGNMENTS = "non-primary alignments"
_SN_KEY_SUPPLEMENTARY_ALIGNMENTS = "supplementary alignments"
_SN_KEY_READS_UNMAPPED = "reads unmapped"


@dataclass(frozen=True)
class BamReadStats:
    """Read counts read from the SN section of `samtools stats`."""

    # "total number of reads in a file, excluding supplementary and secondary reads.
    #  Same number reported by samtools view -c -F 0x900."
    raw_total_sequences: int
    # "number of secondary reads (flag 0x100 (256) set)"
    non_primary_alignments: int
    # "number of supplementary reads (flag 0x800 (2048) set)"
    supplementary_alignments: int
    # "number of unmapped reads (flag 0x4 is set)"
    reads_unmapped: int


def _parse_sn_section(stats_output: str) -> dict[str, str]:
    """Parse the SN section of `samtools stats` output into {key: value}."""
    sn_values: dict[str, str] = {}
    for line in stats_output.splitlines():
        if not line.startswith("SN\t"):
            continue
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        sn_values[cols[1].rstrip(":")] = cols[2]
    return sn_values


def _require_sn_int(sn_values: dict[str, str], key: str, bam: Path) -> int:
    raw = sn_values.get(key)
    if raw is None:
        raise ValueError(f"samtools stats output has no SN '{key}' for {bam}")
    try:
        return int(raw)
    except ValueError:
        raise ValueError(f"Invalid SN '{key}' value '{raw}' in samtools stats for {bam}") from None


def read_bam_read_stats(bam: Path, save_to: Path) -> BamReadStats:
    """Run `samtools stats` on a BAM, save the raw output, and return its read counts."""
    # NOTE: `samtools stats` is used because its "raw total sequences" is documented as the
    # primary line count (same as `samtools view -c -F 0x900`) and the whole breakdown comes
    # from a single pass. `samtools idxstats` cannot exclude secondary alignments and needs an
    # indexed, coordinate-sorted BAM, which the input BAM is not guaranteed to be;
    # `samtools flagstat` has no documented equivalence to -F 0x900; `samtools view -c -F 0x900`
    # would need one extra pass per breakdown item.
    result = subprocess.run(
        ["samtools", "stats", str(bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"samtools stats failed for {bam}: {result.stderr.strip()}")
    save_to.write_text(result.stdout)

    sn_values = _parse_sn_section(result.stdout)
    stats = BamReadStats(
        raw_total_sequences=_require_sn_int(sn_values, _SN_KEY_RAW_TOTAL_SEQUENCES, bam),
        non_primary_alignments=_require_sn_int(sn_values, _SN_KEY_NON_PRIMARY_ALIGNMENTS, bam),
        supplementary_alignments=_require_sn_int(sn_values, _SN_KEY_SUPPLEMENTARY_ALIGNMENTS, bam),
        reads_unmapped=_require_sn_int(sn_values, _SN_KEY_READS_UNMAPPED, bam),
    )
    if stats.raw_total_sequences <= 0:
        raise ValueError(f"No primary reads in {bam} (raw total sequences = 0)")

    logger.info(
        "samtools stats %s: raw total sequences=%d, non-primary alignments=%d, "
        "supplementary alignments=%d, reads unmapped=%d",
        bam,
        stats.raw_total_sequences,
        stats.non_primary_alignments,
        stats.supplementary_alignments,
        stats.reads_unmapped,
    )
    return stats
