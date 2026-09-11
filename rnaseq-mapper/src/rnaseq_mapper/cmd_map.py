# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import csv
import logging
import subprocess
from pathlib import Path

from rnaseq_mapper.external import require_command
from rnaseq_mapper.fastp import FastpResult, read_fastp_after_filtering_total_reads, run_fastp
from rnaseq_mapper.models import LibraryLayout, SampleEntry, Strandedness
from rnaseq_mapper.samtools_stats import read_bam_read_stats
from rnaseq_mapper.star.binary import check_star_version, resolve_star_binary
from rnaseq_mapper.star.mapping import StarMapResult, index_bam, quickcheck_bam, run_star_mapping
from rnaseq_mapper.strand_split import split_and_quantify

logger = logging.getLogger(__name__)

_MERGED_OUTPUT_PREFIX = "all_samples"
MERGED_STRAND_COUNTS_BASENAME = f"{_MERGED_OUTPUT_PREFIX}.strand_counts.tsv"
MERGED_STRAND_DEPTH_BASENAME = f"{_MERGED_OUTPUT_PREFIX}.strand_depth.tsv"


def _parse_sra_info(sra_info_csv: Path) -> dict[str, LibraryLayout]:
    """Parse SraRunInfo.csv. Returns {Run: LibraryLayout}."""
    result: dict[str, LibraryLayout] = {}
    with sra_info_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            run_id = row["Run"].strip()
            layout_raw = row["LibraryLayout"].strip().upper()
            try:
                layout = LibraryLayout(layout_raw)
            except ValueError:
                raise ValueError(
                    f"Unsupported LibraryLayout '{layout_raw}' for {run_id} in {sra_info_csv}",
                ) from None

            result[run_id] = layout
    if not result:
        raise ValueError(f"No samples found in {sra_info_csv}")
    return result


def _parse_strandedness_info(strandedness_csv: Path) -> dict[str, Strandedness]:
    """Parse strandedness CSV. Returns {Run: Strandedness}."""
    result: dict[str, Strandedness] = {}
    with strandedness_csv.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            run_id = row["Run"].strip()
            strand_raw = row["strandedness"].strip()
            try:
                strand = Strandedness(strand_raw)
            except ValueError:
                raise ValueError(
                    f"Unsupported strandedness '{strand_raw}' for {run_id} in {strandedness_csv}. "
                    f"Allowed: {', '.join(s.value for s in Strandedness)}",
                ) from None
            result[run_id] = strand
    return result


def _load_samples(
    sra_info_csv: Path,
    strandedness_csv: Path,
) -> list[SampleEntry]:
    """Join sra_info and strandedness CSVs on Run. Returns list of SampleEntry."""
    layouts = _parse_sra_info(sra_info_csv)
    strandedness_map = _parse_strandedness_info(strandedness_csv)

    missing = set(layouts.keys()) - set(strandedness_map.keys())
    if missing:
        raise ValueError(f"Runs missing from strandedness CSV: {', '.join(sorted(missing))}")

    samples = []
    for run_id, layout in layouts.items():
        samples.append(
            SampleEntry(
                run_id=run_id,
                library_layout=layout,
                strandedness=strandedness_map[run_id],
            ),
        )
    return samples


def _validate_fastq(sample: SampleEntry, fastq_dir: Path) -> tuple[Path, Path | None]:
    """Resolve and validate FASTQ paths for a sample."""
    r1, r2 = sample.fastq_paths(fastq_dir)
    if not r1.is_file():
        raise FileNotFoundError(f"FASTQ R1 not found: {r1}")
    if r2 is not None and not r2.is_file():
        raise FileNotFoundError(f"FASTQ R2 not found: {r2}")
    return r1, r2


def _merge_tsvs(per_sample_files: list[Path], dest: Path) -> None:
    """Concatenate per-sample TSVs into a single file, keeping one header."""
    with dest.open("w") as out:
        for i, src in enumerate(per_sample_files):
            with src.open() as f:
                for line_no, line in enumerate(f):
                    if i > 0 and line_no == 0:
                        continue
                    out.write(line)


def _process_sample(
    sample: SampleEntry,
    fastq_dir: Path,
    fastp_dir: Path,
    star_dir: Path,
    index_dir: Path,
    star_bin: Path,
    threads: int,
    force: bool,
    limit_bam_sort_ram_bytes: int | None = None,
) -> tuple[Path, Path]:
    """Process a single sample through the full pipeline. Returns (counts_tsv, depth_tsv)."""
    r1, r2 = _validate_fastq(sample, fastq_dir)

    fastp_result: FastpResult = run_fastp(
        run_id=sample.run_id,
        raw_r1=r1,
        raw_r2=r2,
        out_dir=fastp_dir,
        threads=threads,
        force=force,
    )

    star_result: StarMapResult = run_star_mapping(
        run_id=sample.run_id,
        fastp_r1=fastp_result.trimmed_r1,
        fastp_r2=fastp_result.trimmed_r2,
        index_dir=index_dir,
        out_dir=star_dir,
        star_bin=star_bin,
        threads=threads,
        force=force,
        limit_bam_sort_ram_bytes=limit_bam_sort_ram_bytes,
    )

    quickcheck_bam(star_result.bam)
    index_bam(star_result.bam, threads)

    total_reads = read_fastp_after_filtering_total_reads(fastp_result.json_report)

    sample_star_dir = star_dir / sample.run_id
    counts_tsv, depth_tsv = split_and_quantify(
        run_id=sample.run_id,
        bam=star_result.bam,
        strandedness=sample.strandedness,
        library_layout=sample.library_layout,
        total_input_reads=total_reads,
        out_dir=sample_star_dir,
        threads=threads,
    )

    return counts_tsv, depth_tsv


def run_map(
    *,
    sra_info_csv: Path,
    strandedness_info_csv: Path,
    index_dir: Path,
    fastq_dir: Path,
    output_dir: Path,
    star_bin: Path | None,
    threads: int,
    force: bool,
    limit_bam_sort_ram_bytes: int | None = None,
) -> None:
    """Run fastp QC, STAR mapping, strand separation, and quantification."""
    if not index_dir.is_dir():
        raise FileNotFoundError(f"STAR index directory not found: {index_dir}")
    if not fastq_dir.is_dir():
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")

    resolved_star_bin = resolve_star_binary(star_bin)
    check_star_version(resolved_star_bin)
    require_command("samtools")
    require_command("fastp")

    samples = _load_samples(sra_info_csv, strandedness_info_csv)
    logger.info("Loaded %d sample(s) from CSV", len(samples))

    genome_name = index_dir.name
    fastp_dir = output_dir / "fastp"
    star_dir = output_dir / "star" / genome_name
    summary_dir = output_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    counts_files: list[Path] = []
    depth_files: list[Path] = []
    failed: list[tuple[str, str]] = []

    for sample in samples:
        logger.info(
            "Processing sample: %s (layout=%s, strandedness=%s)",
            sample.run_id,
            sample.library_layout.value,
            sample.strandedness.value,
        )
        try:
            counts_tsv, depth_tsv = _process_sample(
                sample,
                fastq_dir,
                fastp_dir,
                star_dir,
                index_dir,
                resolved_star_bin,
                threads,
                force,
                limit_bam_sort_ram_bytes=limit_bam_sort_ram_bytes,
            )
            counts_files.append(counts_tsv)
            depth_files.append(depth_tsv)
        except Exception:
            logger.warning("Sample %s failed, skipping", sample.run_id, exc_info=True)
            failed.append((sample.run_id, "see log for details"))

    if counts_files:
        all_counts = summary_dir / MERGED_STRAND_COUNTS_BASENAME
        all_depth = summary_dir / MERGED_STRAND_DEPTH_BASENAME
        _merge_tsvs(counts_files, all_counts)
        _merge_tsvs(depth_files, all_depth)
        logger.info("Merged counts: %s", all_counts)
        logger.info("Merged depth:  %s", all_depth)

    if failed:
        logger.warning(
            "Failed samples (%d/%d): %s",
            len(failed),
            len(samples),
            ", ".join(run_id for run_id, _ in failed),
        )

    if not counts_files and failed:
        raise RuntimeError(f"All {len(failed)} sample(s) failed. No output produced.")


def _convert_sam_to_sorted_bam(sam: Path, out_dir: Path, threads: int) -> Path:
    """Convert SAM to coordinate-sorted BAM via samtools sort."""
    bam = out_dir / f"{sam.stem}.bam"
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)],
        check=True,
    )
    return bam


def _process_alignment_sample(
    sample: SampleEntry,
    alignment_dir: Path,
    bam_dir: Path,
    threads: int,
) -> tuple[Path, Path]:
    """Process a single pre-aligned sample. Returns (counts_tsv, depth_tsv)."""
    alignment_path = sample.alignment_path(alignment_dir)
    sample_dir = bam_dir / sample.run_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    if alignment_path.suffix.lower() == ".sam":
        bam = _convert_sam_to_sorted_bam(alignment_path, sample_dir, threads)
    else:
        bam = alignment_path

    quickcheck_bam(bam)

    bam_stats = read_bam_read_stats(bam, sample_dir / f"{sample.run_id}.stats.txt")
    # Warn unconditionally: whether unmapped reads are in the BAM depends on how it was
    # produced, and that cannot be told apart from a BAM that simply has none.
    logger.warning(
        "%s: RPM denominator=%d = reads kept in %s (unmapped=%d), not the pre-mapping total. "
        "For a comparable RPM, regenerate the BAM with unmapped reads retained "
        "(ex. STAR --outSAMunmapped Within).",
        sample.run_id,
        bam_stats.raw_total_sequences,
        bam,
        bam_stats.reads_unmapped,
    )

    counts_tsv, depth_tsv = split_and_quantify(
        run_id=sample.run_id,
        bam=bam,
        strandedness=sample.strandedness,
        library_layout=sample.library_layout,
        total_input_reads=bam_stats.raw_total_sequences,
        out_dir=sample_dir,
        threads=threads,
    )
    return counts_tsv, depth_tsv


def run_map_from_alignments(
    *,
    sra_info_csv: Path,
    strandedness_info_csv: Path,
    alignment_dir: Path,
    output_dir: Path,
    threads: int,
    force: bool,
) -> None:
    """Run strand separation and quantification from pre-aligned BAM/SAM files."""
    if not alignment_dir.is_dir():
        raise FileNotFoundError(f"Alignment directory not found: {alignment_dir}")

    require_command("samtools")

    samples = _load_samples(sra_info_csv, strandedness_info_csv)
    logger.info("Loaded %d sample(s) from CSV", len(samples))

    for sample in samples:
        sample.alignment_path(alignment_dir)

    bam_dir = output_dir / "bam"
    summary_dir = output_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    counts_files: list[Path] = []
    depth_files: list[Path] = []
    failed: list[tuple[str, str]] = []

    for sample in samples:
        logger.info(
            "Processing alignment sample: %s (layout=%s, strandedness=%s)",
            sample.run_id,
            sample.library_layout.value,
            sample.strandedness.value,
        )
        try:
            counts_tsv, depth_tsv = _process_alignment_sample(
                sample,
                alignment_dir,
                bam_dir,
                threads,
            )
            counts_files.append(counts_tsv)
            depth_files.append(depth_tsv)
        except Exception:
            logger.warning("Sample %s failed, skipping", sample.run_id, exc_info=True)
            failed.append((sample.run_id, "see log for details"))

    if counts_files:
        all_counts = summary_dir / MERGED_STRAND_COUNTS_BASENAME
        all_depth = summary_dir / MERGED_STRAND_DEPTH_BASENAME
        _merge_tsvs(counts_files, all_counts)
        _merge_tsvs(depth_files, all_depth)
        logger.info("Merged counts: %s", all_counts)
        logger.info("Merged depth:  %s", all_depth)

    if failed:
        logger.warning(
            "Failed samples (%d/%d): %s",
            len(failed),
            len(samples),
            ", ".join(run_id for run_id, _ in failed),
        )

    if not counts_files and failed:
        raise RuntimeError(f"All {len(failed)} sample(s) failed. No output produced.")
