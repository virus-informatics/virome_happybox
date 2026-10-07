# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Compute coverage statistics from BLAST hits: per-accession and per-Assembly."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import TYPE_CHECKING

from viral_genome_annotator.visualize.blast_result import BlastHit

if TYPE_CHECKING:
    from viral_genome_annotator.visualize.region_metadata import GroupTaxonomy, ResolvedRegionMeta

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Interval utilities
# ---------------------------------------------------------------------------


def merge_intervals_length(starts: list[int], ends: list[int]) -> int:
    """Merge overlapping intervals and return total covered length."""
    if not starts:
        return 0
    s = [min(a, b) for a, b in zip(starts, ends)]
    e = [max(a, b) for a, b in zip(starts, ends)]
    s, e = zip(*sorted(zip(s, e))) if s else ([], [])
    cur_s = s[0]
    cur_e = e[0]
    total = 0
    for i in range(1, len(s)):
        if s[i] <= cur_e + 1:
            cur_e = max(cur_e, e[i])
        else:
            total += cur_e - cur_s + 1
            cur_s = s[i]
            cur_e = e[i]
    return total + cur_e - cur_s + 1


# ---------------------------------------------------------------------------
# Per-accession (sseqid) coverage
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AccessionCoverage:
    """Coverage of a single accession within one TSV."""

    sseqid: str
    covered_bases: int
    n_hits: int
    best_pident: float
    mean_pident: float
    min_evalue: float
    max_bitscore: float


def compute_accession_coverage(hits: list[BlastHit]) -> list[AccessionCoverage]:
    """Group hits by sseqid and compute coverage statistics."""
    groups: dict[str, list[BlastHit]] = {}
    for hit in hits:
        groups.setdefault(hit.sseqid, []).append(hit)

    results: list[AccessionCoverage] = []
    for sseqid, group_hits in sorted(groups.items()):
        sstarts: list[int] = []
        sends: list[int] = []
        pidents: list[float] = []
        evalues: list[float] = []
        bitscores: list[float] = []
        for h in group_hits:
            sstarts.append(h.sstart)
            sends.append(h.send)
            pidents.append(h.pident)
            evalues.append(h.evalue)
            bitscores.append(h.bitscore)

        results.append(
            AccessionCoverage(
                sseqid=sseqid,
                covered_bases=merge_intervals_length(sstarts, sends),
                n_hits=len(group_hits),
                best_pident=max(pidents),
                mean_pident=round(mean(pidents), 2),
                min_evalue=min(evalues),
                max_bitscore=max(bitscores),
            ),
        )

    logger.info("Computed coverage for %d sseqid groups", len(results))
    return results


# ---------------------------------------------------------------------------
# Per-Assembly coverage
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AssemblyCoverage:
    """Coverage aggregated to the Assembly level within one TSV."""

    group_key: str
    genome_len: int
    covered_bases_total: int
    n_hits_total: int
    best_pident: float
    mean_pident: float
    min_evalue: float
    max_bitscore: float

    @property
    def genome_coverage_pct(self) -> float:
        if self.genome_len <= 0:
            return 0.0
        return round(self.covered_bases_total / self.genome_len * 100.0, 2)


def aggregate_assembly_coverage(
    accession_coverages: list[AccessionCoverage],
    region_meta: list[ResolvedRegionMeta],
) -> list[AssemblyCoverage]:
    """Aggregate accession-level coverage to Assembly (group_key) level."""
    meta_map: dict[str, ResolvedRegionMeta] = {m.sseqid: m for m in region_meta}

    acc_by_group: dict[str, list[AccessionCoverage]] = {}
    for ac in accession_coverages:
        meta = meta_map.get(ac.sseqid)
        if meta is None:
            continue
        acc_by_group.setdefault(meta.group_key, []).append(ac)

    results: list[AssemblyCoverage] = []
    for group_key, accs in sorted(acc_by_group.items()):
        covered_bases: list[int] = []
        n_hits: list[int] = []
        best_pidents: list[float] = []
        mean_pidents: list[float] = []
        min_evalues: list[float] = []
        max_bitscores: list[float] = []
        for ac in accs:
            covered_bases.append(ac.covered_bases)
            n_hits.append(ac.n_hits)
            best_pidents.append(ac.best_pident)
            mean_pidents.append(ac.mean_pident)
            min_evalues.append(ac.min_evalue)
            max_bitscores.append(ac.max_bitscore)

        genome_len = meta_map[accs[0].sseqid].genome_len

        results.append(
            AssemblyCoverage(
                group_key=group_key,
                genome_len=genome_len,
                covered_bases_total=sum(covered_bases),
                n_hits_total=sum(n_hits),
                best_pident=max(best_pidents),
                mean_pident=round(mean(mean_pidents), 2),
                min_evalue=min(min_evalues),
                max_bitscore=max(max_bitscores),
            ),
        )

    logger.info("Aggregated coverage for %d group_key groups", len(results))
    return results


# ---------------------------------------------------------------------------
# min_covered_bases filter
# ---------------------------------------------------------------------------


def filter_by_covered_bases(
    assembly_coverages: list[AssemblyCoverage],
    min_covered_bases: int,
) -> list[AssemblyCoverage]:
    """Return Assemblies whose covered_bases_total >= min_covered_bases."""
    passed = [ac for ac in assembly_coverages if ac.covered_bases_total >= min_covered_bases]
    logger.info(
        "min_covered_bases=%d: %d / %d assemblies passed",
        min_covered_bases,
        len(passed),
        len(assembly_coverages),
    )
    return passed


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

_COVERAGE_CSV_COLUMNS = [
    "sra_id",
    "Assembly",
    "assembly_taxid",
    "Genus",
    "Family",
    "genome_len",
    "covered_bases_total",
    "genome_coverage_pct",
    "n_hits_total",
    "best_pident",
    "mean_pident",
    "min_evalue",
    "max_bitscore",
]


def write_coverage_csv(
    assembly_coverages: list[AssemblyCoverage],
    group_taxonomy: dict[str, GroupTaxonomy],
    sra_id: str,
    output_path: Path,
) -> None:
    """Write Assembly-level coverage to CSV (all Assemblies, no min_covered_bases filter)."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(_COVERAGE_CSV_COLUMNS)

        for ac in assembly_coverages:
            taxonomy = group_taxonomy.get(ac.group_key)
            writer.writerow(
                [
                    sra_id,
                    ac.group_key,
                    "" if taxonomy is None else (taxonomy.assembly_taxid or ""),
                    "" if taxonomy is None else taxonomy.genus,
                    "" if taxonomy is None else taxonomy.family,
                    ac.genome_len,
                    ac.covered_bases_total,
                    ac.genome_coverage_pct,
                    ac.n_hits_total,
                    ac.best_pident,
                    ac.mean_pident,
                    ac.min_evalue,
                    ac.max_bitscore,
                ],
            )

    logger.info("Wrote coverage CSV (%d rows): %s", len(assembly_coverages), output_path)
