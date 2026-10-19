# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Orchestrate blast-result parsing, coverage, annotation, and Plotly visualization."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from viral_genome_annotator.accession import strip_accession_version
from viral_genome_annotator.visualize.annotation import extract_cds_annotations
from viral_genome_annotator.visualize.blast_result import (
    BlastHit,
    extract_tsv_stem,
    parse_blast_tsv,
)
from viral_genome_annotator.visualize.coverage import (
    AssemblyCoverage,
    aggregate_assembly_coverage,
    compute_accession_coverage,
    filter_by_covered_bases,
    write_coverage_csv,
)
from viral_genome_annotator.visualize.plot import write_html_figures
from viral_genome_annotator.visualize.plot_groups import build_plot_groups, group_metas_by_assembly
from viral_genome_annotator.visualize.region_metadata import (
    GroupTaxonomy,
    ResolvedRegionMeta,
    _build_version_stripped_index,
    build_assembly_to_nc_map,
    build_assembly_to_taxid_map,
    build_nc_to_assembly_map,
    extract_region_metadata,
    load_taxonomy_lookup,
    resolve_assembly_for_accession,
    resolve_group_taxonomy,
    resolve_region_metadata,
)

logger = logging.getLogger(__name__)

_BLAST_TSV_GLOB = "*.tblastx.tsv"
_ANNOTATIONS_GFF3_NAME = "annotations.gff3"
_TAXONOMY_TSV_NAME = "taxonomy.tsv"
_ASSEMBLY_SUBDIR = "assembly"
_ASSEMBLY_MAP_JSON_FILENAME = "assembly_map.json"
_VISUALIZE_SUBDIR = "visualize"
_FIGURES_SUBDIR = "figures"


def _write_qseqids_tsv(
    all_hits: list[BlastHit],
    sra_id: str,
    output_path: Path,
) -> None:
    """Write unique qseqid list with sra_id TSV."""
    seen: set[str] = set()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["qseqid", "sra"])
        for h in all_hits:
            if h.qseqid not in seen:
                seen.add(h.qseqid)
                writer.writerow([h.qseqid, sra_id])

    logger.info("Wrote qseqids TSV (%d rows): %s", len(seen), output_path)


def _write_sseqids_tsv(hit_sseqids: set[str], output_path: Path) -> None:
    """Write unique sseqid list TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for sseqid in sorted(hit_sseqids):
            writer.writerow([sseqid])

    logger.info("Wrote sseqids TSV (%d rows): %s", len(hit_sseqids), output_path)


def _write_missing_gff(missing: set[str], output_path: Path) -> None:
    """Write accessions not found in GFF3."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        for sseqid in sorted(missing):
            f.write(f"{sseqid}\n")

    logger.info("Wrote missing GFF accessions (%d): %s", len(missing), output_path)


_BEST_HIT_TSV_COLUMNS = [
    "sra",
    "qseqid",
    "sseqid",
    "sseqid_base",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "start",
    "end",
    "ref_Assembly",
    "accession_taxid",
    "assembly_taxid",
    "ref_Family",
    "ref_Genus",
    "ref_Host",
    "ref_Segment",
    "track_label",
    "covered_bases_total",
    "genome_coverage_pct",
    "best_pident",
]


def _make_track_label(sseqid: str, segment: str | None) -> str:
    if segment:
        return f"{sseqid} | segment {segment}"
    return sseqid


def _write_best_hit_tsv(
    all_hits: list[BlastHit],
    filtered_coverages: list[AssemblyCoverage],
    resolved_meta: list[ResolvedRegionMeta],
    group_taxonomy: dict[str, GroupTaxonomy],
    sra_id: str,
    output_path: Path,
) -> None:
    """Write blast hits for assemblies that passed min_covered_bases filter."""
    meta_map: dict[str, ResolvedRegionMeta] = {m.sseqid: m for m in resolved_meta}
    filtered_group_keys: set[str] = {ac.group_key for ac in filtered_coverages}
    coverage_map: dict[str, AssemblyCoverage] = {ac.group_key: ac for ac in filtered_coverages}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(_BEST_HIT_TSV_COLUMNS)

        for hit in all_hits:
            meta = meta_map.get(hit.sseqid)
            if meta is None:
                continue
            if meta.group_key not in filtered_group_keys:
                continue

            cov = coverage_map.get(meta.group_key)
            if cov is None:
                continue

            tax = group_taxonomy.get(meta.group_key)

            writer.writerow(
                [
                    sra_id,
                    hit.qseqid,
                    hit.sseqid,
                    strip_accession_version(hit.sseqid),
                    hit.pident,
                    hit.length,
                    hit.mismatch,
                    hit.gapopen,
                    hit.qstart,
                    hit.qend,
                    hit.sstart,
                    hit.send,
                    hit.evalue,
                    hit.bitscore,
                    hit.qlen,
                    hit.slen,
                    min(hit.sstart, hit.send),
                    max(hit.sstart, hit.send),
                    meta.group_key,
                    meta.accession_taxid,
                    tax.assembly_taxid if tax and tax.assembly_taxid is not None else "",
                    tax.family if tax else "",
                    tax.genus if tax else "",
                    meta.host or "",
                    meta.segment or "",
                    _make_track_label(hit.sseqid, meta.segment),
                    cov.covered_bases_total,
                    cov.genome_coverage_pct,
                    cov.best_pident,
                ],
            )
            row_count += 1

    logger.info("Wrote best-hit TSV (%d rows): %s", row_count, output_path)


def discover_blast_tsv_files(workdir: Path) -> list[Path]:
    """Glob ``{workdir}/tblastx/*.tblastx.tsv`` and return sorted paths."""
    tblastx_dir = workdir / "tblastx"
    if not tblastx_dir.is_dir():
        raise FileNotFoundError(f"tblastx directory not found: {tblastx_dir}")
    found = sorted(tblastx_dir.glob(_BLAST_TSV_GLOB))
    if not found:
        raise FileNotFoundError(f"No {_BLAST_TSV_GLOB} files in {tblastx_dir}")
    logger.info("Discovered %d blast TSV file(s) in %s", len(found), tblastx_dir)
    return found


def _resolve_metadata(
    hit_sseqids: set[str],
    gff3_path: Path,
    taxonomy_path: Path,
    assembly_map_path: Path,
) -> tuple[list[ResolvedRegionMeta], dict[str, GroupTaxonomy]]:
    """Load reference files and resolve metadata + taxonomy for hit accessions."""
    accession_to_assembly = build_nc_to_assembly_map(assembly_map_path)
    assembly_to_accessions = build_assembly_to_nc_map(assembly_map_path)
    stripped_index = _build_version_stripped_index(accession_to_assembly)

    hit_assemblies: set[str] = set()
    for sseqid in hit_sseqids:
        assembly = resolve_assembly_for_accession(
            sseqid,
            accession_to_assembly,
            stripped_index,
        )
        if assembly is not None:
            hit_assemblies.add(assembly)

    target_sseqids = set(hit_sseqids)
    for assembly in hit_assemblies:
        sibling_accessions = assembly_to_accessions.get(assembly, [])
        target_sseqids.update(sibling_accessions)

    regions = extract_region_metadata(gff3_path, target_sseqids)
    taxonomy = load_taxonomy_lookup(taxonomy_path)
    resolved_meta = resolve_region_metadata(regions, taxonomy, accession_to_assembly)

    assembly_to_taxid = build_assembly_to_taxid_map(assembly_map_path)
    group_taxonomy = resolve_group_taxonomy(resolved_meta, assembly_to_taxid, taxonomy)

    return resolved_meta, group_taxonomy


def run_visualize(
    blast_tsv_path: Path,
    refseq_dir: Path,
    output_dir: Path,
    *,
    min_covered_bases: int = 500,
    max_track_count: int = 15,
    title: str | None = None,
) -> list[Path]:
    """Run the full visualize pipeline for a single BLAST TSV and return written HTML paths."""
    gff3_path = refseq_dir / _ANNOTATIONS_GFF3_NAME
    taxonomy_path = refseq_dir / _TAXONOMY_TSV_NAME
    assembly_map_path = refseq_dir / _ASSEMBLY_SUBDIR / _ASSEMBLY_MAP_JSON_FILENAME

    if not gff3_path.exists():
        raise FileNotFoundError(f"GFF3 not found: {gff3_path}")
    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Taxonomy TSV not found: {taxonomy_path}")
    if not assembly_map_path.exists():
        raise FileNotFoundError(f"Assembly map not found: {assembly_map_path}")

    # --- 1. Parse BLAST TSV ---
    # blast.py:infer_sra_id_from_fasta derived from FASTA header, used as TSV filename
    sra_id = extract_tsv_stem(blast_tsv_path)
    all_hits = parse_blast_tsv(blast_tsv_path)

    if not all_hits:
        logger.warning("No BLAST hits found in %s", blast_tsv_path)
        return []

    logger.info("Loaded %d BLAST hits from %s", len(all_hits), blast_tsv_path)

    visualize_dir = output_dir / _VISUALIZE_SUBDIR
    hit_sseqids = {h.sseqid for h in all_hits}

    metadata_dir = visualize_dir / "metadata"
    _write_qseqids_tsv(all_hits, sra_id, metadata_dir / f"{sra_id}.tblastx.qseqids.tsv")
    _write_sseqids_tsv(hit_sseqids, metadata_dir / f"{sra_id}.tblastx.sseqids.tsv")

    # --- 2. Resolve metadata ---
    resolved_meta, group_taxonomy = _resolve_metadata(
        hit_sseqids,
        gff3_path,
        taxonomy_path,
        assembly_map_path,
    )

    missing_gff = hit_sseqids - {m.sseqid for m in resolved_meta}
    if missing_gff:
        logger.warning("%d sseqids not found in GFF3", len(missing_gff))
        _write_missing_gff(
            missing_gff,
            metadata_dir / f"{sra_id}.missing_gff_accessions.txt",
        )

    # --- 3. Per-accession coverage ---
    accession_coverages = compute_accession_coverage(all_hits)

    # --- 4. Assembly-level coverage + CSV output ---
    assembly_coverages = aggregate_assembly_coverage(accession_coverages, resolved_meta)
    coverage_csv_path = visualize_dir / f"{sra_id}.coverage.csv"
    write_coverage_csv(assembly_coverages, group_taxonomy, sra_id, coverage_csv_path)

    # --- 5. min_covered_bases filter ---
    filtered_coverages = filter_by_covered_bases(assembly_coverages, min_covered_bases)

    if not filtered_coverages:
        logger.warning(
            "No assemblies passed min_covered_bases=%d filter",
            min_covered_bases,
        )
        return []

    # --- 6. blast_best_hit.tsv ---
    best_hit_tsv_path = visualize_dir / f"{sra_id}.blast_best_hit.tsv"
    _write_best_hit_tsv(
        all_hits,
        filtered_coverages,
        resolved_meta,
        group_taxonomy,
        sra_id,
        best_hit_tsv_path,
    )

    # --- 7. PlotGroup assembly ---
    metas_by_assembly = group_metas_by_assembly(resolved_meta)
    cds_annotations = extract_cds_annotations(gff3_path, {m.sseqid for m in resolved_meta})
    plot_groups = build_plot_groups(
        all_hits,
        filtered_coverages,
        metas_by_assembly,
        group_taxonomy,
        cds_annotations,
        sra_id,
        max_track_count=max_track_count,
    )

    # --- 8. Render HTML figures ---
    figures_dir = visualize_dir / _FIGURES_SUBDIR
    html_paths = write_html_figures(plot_groups, figures_dir, title=title)

    logger.info(
        "Visualize pipeline: %d figure(s), coverage CSV at %s, best-hit TSV at %s",
        len(html_paths),
        coverage_csv_path,
        best_hit_tsv_path,
    )
    return html_paths
