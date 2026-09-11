# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

from rnaseq_mapper.visualize._style import (
    STRAND_COLORS,
    STRAND_FALLBACK_COLOR,
    resolve_strand_order,
)


def plot_strand_coverage(
    sample_name: str,
    sample_depth: pd.DataFrame,
    segment_order: list[str],
) -> plt.Figure:
    """Generate a strand coverage figure for one sample.

    Args:
        sample_name: sample identifier (used in suptitle).
        sample_depth: DataFrame filtered to this sample,
            columns: segment, position, strand, coverage.
        segment_order: ordered list of segment names.

    Returns:
        The matplotlib Figure. Caller is responsible for saving and closing.
    """
    strand_order = resolve_strand_order(list(dict.fromkeys(sample_depth["strand"].tolist())))
    n_segments = len(segment_order)
    fig_height = max(3.5, 2.2 * n_segments)
    fig, axes = plt.subplots(
        n_segments,
        1,
        figsize=(11, fig_height),
        squeeze=False,
        constrained_layout=True,
    )
    axes = axes[:, 0]

    for ax, segment in zip(axes, segment_order):
        seg_df = sample_depth.loc[sample_depth["segment"] == segment]
        for strand in strand_order:
            strand_df = seg_df.loc[seg_df["strand"] == strand]
            if strand_df.empty:
                continue
            x = strand_df["position"].to_numpy()
            y = strand_df["coverage"].to_numpy().astype(np.float64)
            if strand == "antisense":
                y = -y
            ax.fill_between(
                x,
                0,
                y,
                color=STRAND_COLORS.get(strand, STRAND_FALLBACK_COLOR),
                alpha=0.85,
                linewidth=0,
            )

        ax.axhline(0, color="#333333", linewidth=0.6)
        ax.set_title(segment, loc="left", fontsize=10)
        ax.set_xlabel("Genome position")
        ax.set_ylabel("Coverage")
        # Antisense coverage is plotted downward but values are inherently positive.
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):g}"))
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(sample_name, fontsize=13)
    present_strands = set(sample_depth["strand"].tolist())
    fig.legend(
        handles=[
            Patch(facecolor=STRAND_COLORS.get(s, STRAND_FALLBACK_COLOR), label=s)
            for s in strand_order
            if s in present_strands
        ],
        loc="upper right",
        frameon=False,
    )

    return fig
