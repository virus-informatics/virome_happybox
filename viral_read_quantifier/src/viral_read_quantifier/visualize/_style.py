# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

STRAND_COLORS: dict[str, str] = {
    "sense": "#2B6CB0",
    "antisense": "#C2410C",
    "total": "#4B5563",
}
STRAND_FALLBACK_COLOR = "#6B7280"

_PREFERRED_STRAND_ORDER = ["sense", "antisense", "total"]


def resolve_strand_order(strands: list[str]) -> list[str]:
    ordered = [s for s in _PREFERRED_STRAND_ORDER if s in strands]
    ordered.extend(s for s in strands if s not in ordered)
    return ordered
