# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from viral_read_quantifier.visualize.coverage import plot_strand_coverage
from viral_read_quantifier.visualize.rpm import plot_strand_rpm

logger = logging.getLogger(__name__)

_OUTPUT_FORMATS = ("svg", "pdf")


def _save_figure(fig: plt.Figure, img_dir: Path, stem: str) -> None:
    try:
        for fmt in _OUTPUT_FORMATS:
            path = img_dir / f"{stem}.{fmt}"
            fig.savefig(str(path), format=fmt)
            logger.info("Wrote %s: %s", fmt.upper(), path)
    finally:
        plt.close(fig)


def run_visualize(
    *,
    depth_tsv: Path,
    counts_tsv: Path,
    output_dir: Path,
) -> None:
    """Read merged TSVs and generate per-sample images (SVG + PDF)."""
    if not depth_tsv.is_file():
        raise FileNotFoundError(f"Depth TSV not found: {depth_tsv}")
    if not counts_tsv.is_file():
        raise FileNotFoundError(f"Counts TSV not found: {counts_tsv}")

    depth = pd.read_csv(depth_tsv, sep="\t")
    counts = pd.read_csv(counts_tsv, sep="\t")

    if depth.empty and counts.empty:
        raise ValueError("Both depth and counts TSVs are empty — nothing to plot")

    segment_order = list(dict.fromkeys(depth["segment"].tolist()))
    sample_order = list(dict.fromkeys(depth["sample"].tolist()))

    if not sample_order and not counts.empty:
        sample_order = list(dict.fromkeys(counts["sample"].tolist()))
    if not segment_order and not counts.empty:
        segment_order = list(dict.fromkeys(counts["segment"].tolist()))

    img_dir = output_dir / "img"
    img_dir.mkdir(parents=True, exist_ok=True)

    for sample in sample_order:
        sample_depth = depth.loc[depth["sample"] == sample]
        if not sample_depth.empty:
            fig = plot_strand_coverage(sample, sample_depth, segment_order)
            _save_figure(fig, img_dir, f"{sample}.strand_coverage")

        sample_counts = counts.loc[counts["sample"] == sample]
        if not sample_counts.empty:
            fig = plot_strand_rpm(sample, sample_counts, segment_order)
            _save_figure(fig, img_dir, f"{sample}.strand_rpm")

    logger.info(
        "Generated images for %d sample(s) (SVG + PDF) in %s",
        len(sample_order),
        img_dir,
    )
