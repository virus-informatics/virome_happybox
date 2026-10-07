# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Select the genus-best Assembly per figure and assemble PlotGroups."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

from viral_genome_annotator.visualize.plot import PlotGroup, TrackData

if TYPE_CHECKING:
    from viral_genome_annotator.visualize.annotation import CdsAnnotation
    from viral_genome_annotator.visualize.blast_result import BlastHit
    from viral_genome_annotator.visualize.coverage import AssemblyCoverage
    from viral_genome_annotator.visualize.region_metadata import GroupTaxonomy, ResolvedRegionMeta

logger = logging.getLogger(__name__)


def select_genus_best(
    filtered_coverages: list[AssemblyCoverage],
    group_taxonomy: dict[str, GroupTaxonomy],
    sra_id: str,
) -> pd.DataFrame:
    """Keep the Assembly with the highest covered_bases_total per genus."""
    rows = []
    for ac in filtered_coverages:
        tax = group_taxonomy[ac.group_key]
        rows.append(
            {
                "sra": sra_id,
                "ref_Assembly": ac.group_key,
                "ref_Family": tax.family,
                "ref_Genus": tax.genus,
                "assembly_taxid": tax.assembly_taxid,
                "covered_bases_total": ac.covered_bases_total,
                "genome_coverage_pct": ac.genome_coverage_pct,
                "best_pident": ac.best_pident,
            },
        )
    if not rows:
        return pd.DataFrame()

    # Without dtype=object, assembly_taxid mixes int with None and pandas turns it into NaN
    assemblies = pd.DataFrame(rows, dtype=object)

    # Stable sort keeps the input order among ties, matching slice_max(with_ties = FALSE)
    ranked = assemblies.sort_values("covered_bases_total", ascending=False, kind="stable")
    best = ranked.groupby("ref_Genus").head(1)

    # NOTE: genus="" means taxonomy.tsv lacked the taxid; distinct viruses may be collapsed
    empty_genus_rows = ranked[ranked["ref_Genus"] == ""]
    if len(empty_genus_rows) > 1:
        kept = empty_genus_rows.iloc[0]["ref_Assembly"]
        for _, row in empty_genus_rows.iloc[1:].iterrows():
            logger.warning(
                "Assembly %s (covered_bases=%s) dropped in favour of %s: both have unknown genus",
                row["ref_Assembly"],
                row["covered_bases_total"],
                kept,
            )

    return best.sort_values(["ref_Family", "ref_Genus"], kind="stable").reset_index(drop=True)


def group_metas_by_assembly(
    resolved_meta: list[ResolvedRegionMeta],
) -> dict[str, list[ResolvedRegionMeta]]:
    """Index accession metadata by Assembly.

    Accessions with no hit belong here too, so that one hit anywhere on the Assembly
    puts every one of its segments on a track.
    """
    metas_by_assembly: dict[str, list[ResolvedRegionMeta]] = {}
    for meta in resolved_meta:
        metas_by_assembly.setdefault(meta.group_key, []).append(meta)
    return metas_by_assembly


def build_plot_groups(
    all_hits: list[BlastHit],
    filtered_coverages: list[AssemblyCoverage],
    metas_by_assembly: dict[str, list[ResolvedRegionMeta]],
    group_taxonomy: dict[str, GroupTaxonomy],
    cds_annotations: list[CdsAnnotation],
    sra_id: str,
    *,
    max_track_count: int,
) -> list[PlotGroup]:
    """Assemble one PlotGroup per genus-best Assembly."""
    combos = select_genus_best(filtered_coverages, group_taxonomy, sra_id)

    hits_by_sseqid: dict[str, list[BlastHit]] = {}
    for hit in all_hits:
        hits_by_sseqid.setdefault(hit.sseqid, []).append(hit)

    cds_by_sseqid: dict[str, list[CdsAnnotation]] = {}
    for cds in cds_annotations:
        cds_by_sseqid.setdefault(cds.sseqid, []).append(cds)

    groups: list[PlotGroup] = []
    for _, combo in combos.iterrows():
        assembly_metas = metas_by_assembly.get(combo["ref_Assembly"], [])
        if not assembly_metas:
            continue

        # Without dtype=object, ref_Segment and ref_Host mix str with None and pandas
        # turns the None into NaN, which is truthy
        track_tbl = pd.DataFrame(
            [
                {
                    "sseqid": m.sseqid,
                    "ref_Segment": m.segment,
                    "slen": m.slen,
                    "accession_taxid": m.accession_taxid,
                    "ref_Host": m.host,
                }
                for m in assembly_metas
            ],
            dtype=object,
        )
        track_tbl["segment_order"] = (
            track_tbl["ref_Segment"].map({"L": 1, "M": 2, "S": 3}).fillna(100)
        )
        track_tbl = track_tbl.sort_values(["segment_order", "sseqid"])

        if len(track_tbl) > max_track_count:
            logger.warning(
                "[%s] %s: showing %d of %d reference sequences."
                " Pass --max-track-count %d to show all",
                combo["sra"],
                combo["ref_Assembly"],
                max_track_count,
                len(track_tbl),
                len(track_tbl),
            )

        track_tbl = track_tbl.head(max_track_count)

        if track_tbl.empty:
            continue

        tracks: list[TrackData] = []
        for _, track in track_tbl.iterrows():
            sseqid = track["sseqid"]
            hits = hits_by_sseqid.get(sseqid, [])
            # Plotly legend follows add_trace order; sort by genomic coordinate
            accession_cds_annotations = sorted(
                cds_by_sseqid.get(sseqid, []),
                key=lambda c: c.min_start,
            )
            tracks.append(
                TrackData(
                    sseqid=sseqid,
                    slen=track["slen"],
                    segment=track["ref_Segment"],
                    accession_taxid=track["accession_taxid"],
                    host=track["ref_Host"],
                    hits=hits,
                    cds_annotations=accession_cds_annotations,
                ),
            )

        groups.append(
            PlotGroup(
                tsv_stem=combo["sra"],
                assembly=combo["ref_Assembly"],
                assembly_taxid=combo["assembly_taxid"],
                family=combo["ref_Family"],
                genus=combo["ref_Genus"],
                tracks=tracks,
                coverage_pct=combo["genome_coverage_pct"],
                covered_bases=combo["covered_bases_total"],
                best_pident=combo["best_pident"],
            ),
        )

    logger.info("Built %d plot groups", len(groups))
    return groups
