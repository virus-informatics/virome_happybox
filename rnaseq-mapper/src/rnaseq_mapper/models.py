# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class Strandedness(Enum):
    UNSTRANDED = "fr-unstranded"
    FR_FIRSTSTRAND = "fr-firststrand"
    FR_SECONDSTRAND = "fr-secondstrand"


class LibraryLayout(Enum):
    SINGLE = "SINGLE"
    PAIRED = "PAIRED"


@dataclass(frozen=True)
class SampleEntry:
    run_id: str
    library_layout: LibraryLayout
    strandedness: Strandedness

    def fastq_paths(self, fastq_dir: Path) -> tuple[Path, Path | None]:
        if self.library_layout == LibraryLayout.PAIRED:
            return (
                fastq_dir / f"{self.run_id}_1.fastq.gz",
                fastq_dir / f"{self.run_id}_2.fastq.gz",
            )
        return (fastq_dir / f"{self.run_id}.fastq.gz", None)

    def alignment_path(self, alignment_dir: Path) -> Path:
        bam = alignment_dir / f"{self.run_id}.bam"
        if bam.is_file():
            return bam
        sam = alignment_dir / f"{self.run_id}.sam"
        if sam.is_file():
            return sam
        raise FileNotFoundError(f"No .bam or .sam found for {self.run_id} in {alignment_dir}")
