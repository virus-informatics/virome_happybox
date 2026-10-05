# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
import subprocess
from io import TextIOWrapper
from pathlib import Path

from viral_read_quantifier.models import LibraryLayout, Strandedness

logger = logging.getLogger(__name__)


def _rpm(mapped: int, total: int) -> float:
    return mapped / total * 1_000_000.0 if total > 0 else 0.0


# SAM FLAG bit fields (decode with: samtools flags <value>)
_FLAG_REVERSE = 0x10  # 16:  read mapped to reverse strand
_FLAG_READ1 = 0x40  # 64:  first in pair
_FLAG_READ2 = 0x80  # 128: second in pair
_FLAG_READ1_REVERSE = 0x50  # 80:  first in pair + reverse strand
_FLAG_READ2_REVERSE = 0x90  # 144: second in pair + reverse strand

_EXCLUDE_NON_PRIMARY = 0x904  # 2308: UNMAP | SECONDARY | SUPPLEMENTARY
_EXCLUDE_NON_PRIMARY_AND_REVERSE = 0x914  # 2324: above + REVERSE

_PE_FLAG_TABLE: dict[Strandedness, dict[str, list[tuple[int, int]]]] = {
    Strandedness.FR_FIRSTSTRAND: {
        "sense": [
            (_FLAG_READ1_REVERSE, _EXCLUDE_NON_PRIMARY),
            (_FLAG_READ2, _EXCLUDE_NON_PRIMARY_AND_REVERSE),
        ],
        "antisense": [
            (_FLAG_READ1, _EXCLUDE_NON_PRIMARY_AND_REVERSE),
            (_FLAG_READ2_REVERSE, _EXCLUDE_NON_PRIMARY),
        ],
    },
    Strandedness.FR_SECONDSTRAND: {
        "sense": [
            (_FLAG_READ1, _EXCLUDE_NON_PRIMARY_AND_REVERSE),
            (_FLAG_READ2_REVERSE, _EXCLUDE_NON_PRIMARY),
        ],
        "antisense": [
            (_FLAG_READ1_REVERSE, _EXCLUDE_NON_PRIMARY),
            (_FLAG_READ2, _EXCLUDE_NON_PRIMARY_AND_REVERSE),
        ],
    },
}

_SE_FLAG_TABLE: dict[Strandedness, dict[str, tuple[int | None, int]]] = {
    Strandedness.FR_FIRSTSTRAND: {
        "sense": (_FLAG_REVERSE, _EXCLUDE_NON_PRIMARY),
        "antisense": (None, _EXCLUDE_NON_PRIMARY_AND_REVERSE),
    },
    Strandedness.FR_SECONDSTRAND: {
        "sense": (None, _EXCLUDE_NON_PRIMARY_AND_REVERSE),
        "antisense": (_FLAG_REVERSE, _EXCLUDE_NON_PRIMARY),
    },
}


def _samtools_view(
    src: Path,
    dest: Path,
    *,
    include_flags: int | None,
    exclude_flags: int,
    threads: int,
) -> None:
    cmd = ["samtools", "view", "-@", str(threads), "-b", "-o", str(dest)]
    if include_flags is not None:
        cmd.extend(["-f", str(include_flags)])
    cmd.extend(["-F", str(exclude_flags), str(src)])
    subprocess.run(cmd, check=True)


def _samtools_merge(parts: list[Path], dest: Path, threads: int) -> None:
    subprocess.run(
        ["samtools", "merge", "-f", "-@", str(threads), str(dest)] + [str(p) for p in parts],
        check=True,
    )
    for p in parts:
        p.unlink(missing_ok=True)


def _samtools_index(bam: Path, threads: int) -> None:
    subprocess.run(["samtools", "index", "-@", str(threads), str(bam)], check=True)


def _samtools_idxstats(bam: Path, save_to: Path) -> list[tuple[str, int, int]]:
    result = subprocess.run(
        ["samtools", "idxstats", str(bam)],
        capture_output=True,
        text=True,
        check=True,
    )
    save_to.write_text(result.stdout)
    segments = []
    for line in result.stdout.strip().splitlines():
        cols = line.split("\t")
        if cols[0] != "*":
            segments.append((cols[0], int(cols[1]), int(cols[2])))
    return segments


def _stream_depth(
    run_id: str,
    bam: Path,
    strand_label: str,
    out: TextIOWrapper,
    save_to: Path,
) -> None:
    proc = subprocess.Popen(
        ["samtools", "depth", "-aa", "-d", "0", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    with save_to.open("w") as raw:
        for line in proc.stdout:
            raw.write(line)
            cols = line.rstrip("\n").split("\t")
            out.write(f"{run_id}\t{cols[0]}\t{cols[1]}\t{strand_label}\t{cols[2]}\n")
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools depth failed for {bam}")


def _extract_total_bam(
    bam: Path,
    out_dir: Path,
    run_id: str,
    threads: int,
) -> Path:
    # Not in _PE/_SE_FLAG_TABLE: unstranded uses the same filter for SE and PE.
    total_bam = out_dir / f"{run_id}.total.bam"
    _samtools_view(
        bam,
        total_bam,
        include_flags=None,
        exclude_flags=_EXCLUDE_NON_PRIMARY,
        threads=threads,
    )
    _samtools_index(total_bam, threads)
    return total_bam


def _extract_stranded_bams(
    bam: Path,
    out_dir: Path,
    run_id: str,
    strandedness: Strandedness,
    library_layout: LibraryLayout,
    threads: int,
) -> tuple[Path, Path]:
    sense_bam = out_dir / f"{run_id}.sense.bam"
    antisense_bam = out_dir / f"{run_id}.antisense.bam"

    if library_layout == LibraryLayout.PAIRED:
        flags = _PE_FLAG_TABLE[strandedness]
        for strand_label, dest in [("sense", sense_bam), ("antisense", antisense_bam)]:
            parts = []
            for i, (inc, exc) in enumerate(flags[strand_label]):
                part = out_dir / f"{run_id}.{strand_label}.part{i}.bam"
                _samtools_view(bam, part, include_flags=inc, exclude_flags=exc, threads=threads)
                parts.append(part)
            _samtools_merge(parts, dest, threads)
    else:
        flags = _SE_FLAG_TABLE[strandedness]
        for strand_label, dest in [("sense", sense_bam), ("antisense", antisense_bam)]:
            inc, exc = flags[strand_label]
            _samtools_view(bam, dest, include_flags=inc, exclude_flags=exc, threads=threads)

    _samtools_index(sense_bam, threads)
    _samtools_index(antisense_bam, threads)
    return sense_bam, antisense_bam


def _write_total_counts(
    run_id: str,
    segments: list[tuple[str, int, int]],
    total_input_reads: int,
    strandedness: Strandedness,
    out: TextIOWrapper,
) -> None:
    confidence = 1 if strandedness != Strandedness.UNSTRANDED else 0
    action = "sense_antisense" if strandedness != Strandedness.UNSTRANDED else "total_only"
    for seg_name, seg_len, mapped in segments:
        rpm = _rpm(mapped, total_input_reads)
        out.write(
            f"{run_id}\t{seg_name}\t{seg_len}\ttotal\t{mapped}\t"
            f"{total_input_reads}\t{rpm}\t{strandedness.value}\t{confidence}\t{action}\n",
        )


def _write_stranded_counts(
    run_id: str,
    sense_segments: list[tuple[str, int, int]],
    antisense_segments: list[tuple[str, int, int]],
    total_input_reads: int,
    strandedness: Strandedness,
    out: TextIOWrapper,
) -> None:
    antisense_map = {name: mapped for name, _, mapped in antisense_segments}
    for seg_name, seg_len, sense_mapped in sense_segments:
        antisense_mapped = antisense_map.get(seg_name, 0)
        sense_rpm = _rpm(sense_mapped, total_input_reads)
        anti_rpm = _rpm(antisense_mapped, total_input_reads)
        out.write(
            f"{run_id}\t{seg_name}\t{seg_len}\tsense\t{sense_mapped}\t"
            f"{total_input_reads}\t{sense_rpm}\t{strandedness.value}\t1\tsense_antisense\n",
        )
        out.write(
            f"{run_id}\t{seg_name}\t{seg_len}\tantisense\t{antisense_mapped}\t"
            f"{total_input_reads}\t{anti_rpm}\t{strandedness.value}\t1\tsense_antisense\n",
        )


_COUNTS_HEADER = (
    "sample\tsegment\tlength\tstrand\tmapped_reads\t"
    "total_input_reads\tRPM\tstrandedness\tstrandedness_confidence\taction\n"
)
_DEPTH_HEADER = "sample\tsegment\tposition\tstrand\tcoverage\n"


def split_and_quantify(
    *,
    run_id: str,
    bam: Path,
    strandedness: Strandedness,
    library_layout: LibraryLayout,
    total_input_reads: int,
    out_dir: Path,
    threads: int,
) -> tuple[Path, Path]:
    """Split BAM by strand and quantify. Returns (counts_tsv, depth_tsv)."""
    counts_tsv = out_dir / f"{run_id}.strand_counts.tsv"
    depth_tsv = out_dir / f"{run_id}.strand_depth.tsv"

    if strandedness == Strandedness.UNSTRANDED:
        total_bam = _extract_total_bam(bam, out_dir, run_id, threads)
        segments = _samtools_idxstats(total_bam, out_dir / f"{run_id}.total.idxstats.tsv")

        with counts_tsv.open("w") as f:
            f.write(_COUNTS_HEADER)
            _write_total_counts(run_id, segments, total_input_reads, strandedness, f)

        with depth_tsv.open("w") as f:
            f.write(_DEPTH_HEADER)
            _stream_depth(run_id, total_bam, "total", f, out_dir / f"{run_id}.total.depth.tsv")
    else:
        sense_bam, antisense_bam = _extract_stranded_bams(
            bam,
            out_dir,
            run_id,
            strandedness,
            library_layout,
            threads,
        )
        sense_stats = _samtools_idxstats(sense_bam, out_dir / f"{run_id}.sense.idxstats.tsv")
        antisense_stats = _samtools_idxstats(
            antisense_bam,
            out_dir / f"{run_id}.antisense.idxstats.tsv",
        )

        with counts_tsv.open("w") as f:
            f.write(_COUNTS_HEADER)
            _write_stranded_counts(
                run_id,
                sense_stats,
                antisense_stats,
                total_input_reads,
                strandedness,
                f,
            )

        with depth_tsv.open("w") as f:
            f.write(_DEPTH_HEADER)
            _stream_depth(run_id, sense_bam, "sense", f, out_dir / f"{run_id}.sense.depth.tsv")
            _stream_depth(
                run_id,
                antisense_bam,
                "antisense",
                f,
                out_dir / f"{run_id}.antisense.depth.tsv",
            )

    return counts_tsv, depth_tsv
