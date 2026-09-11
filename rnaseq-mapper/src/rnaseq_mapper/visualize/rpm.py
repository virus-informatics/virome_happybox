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


def plot_strand_rpm(
    sample_name: str,
    sample_counts: pd.DataFrame,
    segment_order: list[str],
) -> plt.Figure:
    """Generate an RPM bar chart figure for one sample.

    Args:
        sample_name: sample identifier (used in title).
        sample_counts: DataFrame filtered to this sample,
            columns: segment, strand, RPM.
        segment_order: ordered list of segment names.

    Returns:
        The matplotlib Figure. Caller is responsible for saving and closing.
    """
    strand_order = resolve_strand_order(list(dict.fromkeys(sample_counts["strand"].tolist())))
    n_segments = len(segment_order)
    x = np.arange(n_segments)

    fig, ax = plt.subplots(figsize=(max(3.6, n_segments * 1.2), 4.5), constrained_layout=True)

    for strand in strand_order:
        y_values = []
        for segment in segment_order:
            row = sample_counts.loc[
                (sample_counts["segment"] == segment) & (sample_counts["strand"] == strand)
            ]
            rpm = 0.0 if row.empty else float(row["RPM"].iloc[0])
            if strand == "antisense":
                rpm = -rpm
            y_values.append(rpm)
        ax.bar(
            x,
            y_values,
            color=STRAND_COLORS.get(strand, STRAND_FALLBACK_COLOR),
            alpha=0.88,
            width=0.72,
        )

    ax.axhline(0, color="#333333", linewidth=0.6)
    ax.set_title(sample_name, fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(segment_order, rotation=35, ha="right")
    ax.set_ylabel("RPM")
    # Antisense RPM is plotted downward but values are inherently positive.
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):g}"))
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    present_strands = set(sample_counts["strand"].tolist())
    ax.legend(
        handles=[
            Patch(facecolor=STRAND_COLORS.get(s, STRAND_FALLBACK_COLOR), label=s)
            for s in strand_order
            if s in present_strands
        ],
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
    )

    return fig
