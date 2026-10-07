# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

"""Generate interactive Plotly HTML figures for genome annotation visualization."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import plotly.graph_objects as go
from plotly.colors import sequential, unlabel_rgb
from plotly.subplots import make_subplots

from viral_genome_annotator.visualize.annotation import CdsAnnotation
from viral_genome_annotator.visualize.blast_result import BlastHit

logger = logging.getLogger(__name__)

_PIDENT_COLOR_MIN = 50.0
_PIDENT_COLOR_MAX = 100.0
_PIDENT_PALETTE: list[str] = [*sequential.PuRd[2:9], "rgb(0,0,0)"]
_PIDENT_COLORSCALE: list[tuple[float, str]] = [
    (index / (len(_PIDENT_PALETTE) - 1), color) for index, color in enumerate(_PIDENT_PALETTE)
]

_CDS_BOX_Y_BOTTOM = 0.60
_GENOME_LINE_Y = (_CDS_BOX_Y_BOTTOM + 1.0) / 2
_PIDENT_RECT_HEIGHT = 0.20
_PIDENT_RECT_Y_BOTTOM = (_CDS_BOX_Y_BOTTOM - _PIDENT_RECT_HEIGHT) / 2
_PIDENT_RECT_Y_TOP = _PIDENT_RECT_Y_BOTTOM + _PIDENT_RECT_HEIGHT
_CDS_LABEL_FONT_SIZE_PX = 20
_CDS_LABEL_CHAR_WIDTH_FACTOR = 0.6
_CDS_LABEL_PLOT_WIDTH_PX = 700
_CDS_LEGEND_ENTRY_HEIGHT_PX = 8
_CDS_LEGEND_OFFSET_PX = 110
_MARGIN_TOP_PX = 100
_MARGIN_BOTTOM_BASE_PX = 120
_MARGIN_LEFT_PX = 60
_MARGIN_RIGHT_PX = 40
_TRACK_CONTENT_HEIGHT_PX = 200
_TRACK_SPACING_PX = 60
_CDS_LABEL_TOGGLE_JS_PATH = Path(__file__).parent / "cds_label_toggle.js"
_TEMPLATE_VARS_START_MARKER = "/* __TEMPLATE_VARS_START__ */"
_TEMPLATE_VARS_END_MARKER = "/* __TEMPLATE_VARS_END__ */"


@dataclass(frozen=True, slots=True)
class TrackData:
    sseqid: str
    slen: int
    segment: str | None
    accession_taxid: int
    host: str | None
    hits: list[BlastHit]
    cds_annotations: list[CdsAnnotation]


@dataclass(frozen=True, slots=True)
class PlotGroup:
    """One figure = one input TSV × one group.

    A group is an Assembly, or a single accession that NCBI Assembly does not cover.
    """

    tsv_stem: str
    assembly: str
    assembly_taxid: int | None
    family: str | None
    genus: str | None
    tracks: list[TrackData]
    coverage_pct: float
    covered_bases: int
    best_pident: float


def _track_label(sseqid: str, segment: str | None) -> str:
    if segment:
        return f"{sseqid} | segment {segment}"
    return sseqid


def _cds_label_fits(product: str, cds_width_bp: int, slen: int) -> bool:
    """Estimate whether a product label fits inside a CDS box at initial zoom."""
    bp_per_px = slen / _CDS_LABEL_PLOT_WIDTH_PX
    char_width_bp = bp_per_px * _CDS_LABEL_FONT_SIZE_PX * _CDS_LABEL_CHAR_WIDTH_FACTOR
    label_width_bp = len(product) * char_width_bp
    return label_width_bp <= cds_width_bp


def _norm_pident(pident: float) -> float:
    clamped = max(_PIDENT_COLOR_MIN, min(_PIDENT_COLOR_MAX, pident))
    return (clamped - _PIDENT_COLOR_MIN) / (_PIDENT_COLOR_MAX - _PIDENT_COLOR_MIN)


def _pident_to_color(pident: float) -> str:
    norm = _norm_pident(pident)
    scale = _PIDENT_COLORSCALE
    for i in range(len(scale) - 1):
        lo_pos, lo_color = scale[i]
        hi_pos, hi_color = scale[i + 1]
        if lo_pos <= norm <= hi_pos:
            t = (norm - lo_pos) / (hi_pos - lo_pos) if hi_pos != lo_pos else 0.0
            lo_rgb = unlabel_rgb(lo_color)
            hi_rgb = unlabel_rgb(hi_color)
            r = round(lo_rgb[0] + t * (hi_rgb[0] - lo_rgb[0]))
            g = round(lo_rgb[1] + t * (hi_rgb[1] - lo_rgb[1]))
            b = round(lo_rgb[2] + t * (hi_rgb[2] - lo_rgb[2]))
            return f"rgb({r},{g},{b})"
    return "rgb(0,0,0)"


def _add_track(
    fig: go.Figure,
    track: TrackData,
    row: int,
    *,
    coverage_label: str | None = None,
) -> None:
    x_ref = "x" if row == 1 else f"x{row}"

    fig.add_trace(
        go.Scatter(
            x=[0, track.slen],
            y=[_GENOME_LINE_Y, _GENOME_LINE_Y],
            mode="lines",
            line=dict(color="black", width=2),
            showlegend=False,
            hoverinfo="skip",
        ),
        row=row,
        col=1,
    )

    for hit in sorted(track.hits, key=lambda h: h.pident):
        sstart = min(hit.sstart, hit.send)
        send = max(hit.sstart, hit.send)
        color = _pident_to_color(hit.pident)

        fig.add_trace(
            go.Scatter(
                x=[sstart, send, send, sstart, sstart],
                y=[
                    _PIDENT_RECT_Y_BOTTOM,
                    _PIDENT_RECT_Y_BOTTOM,
                    _PIDENT_RECT_Y_TOP,
                    _PIDENT_RECT_Y_TOP,
                    _PIDENT_RECT_Y_BOTTOM,
                ],
                fill="toself",
                fillcolor=color,
                line=dict(width=0),
                mode="lines",
                showlegend=False,
                hoveron="fills",
                hovertemplate=(
                    f"sstart: {sstart}<br>send: {send}<br>"
                    f"pident: {hit.pident:.1f}%<br>"
                    f"evalue: {hit.evalue:.2e}<br>"
                    f"qseqid: {hit.qseqid}"
                    "<extra></extra>"
                ),
            ),
            row=row,
            col=1,
        )

    for cds in track.cds_annotations:
        xs: list[float | None] = []
        ys: list[float | None] = []
        for start, end in cds.parts:
            if xs:
                xs.append(None)
                ys.append(None)
            xs.extend([start, end, end, start, start])
            ys.extend([_CDS_BOX_Y_BOTTOM, _CDS_BOX_Y_BOTTOM, 1.0, 1.0, _CDS_BOX_Y_BOTTOM])

        ranges_str = ", ".join(f"{s:,}–{e:,}" for s, e in cds.parts)
        hover_text = f"{cds.product}<br>{ranges_str} bp" + (
            f" | gene: {cds.gene}" if cds.gene else ""
        )

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                fill="toself",
                fillcolor="rgb(216,216,216)",
                line=dict(color="rgb(102,102,102)", width=1),
                mode="lines",
                name=cds.product,
                showlegend=True,
                legend="legend2",
                legendgroup=f"cds_{track.sseqid}_{cds.cds_id}",
                # hovertemplate is ignored on fill="toself" traces; use hoverinfo+text instead
                hoverinfo="text",
                text=hover_text,
            ),
            row=row,
            col=1,
        )

        largest_start, largest_end = max(cds.parts, key=lambda p: p[1] - p[0])
        label_visible = _cds_label_fits(cds.product, largest_end - largest_start, track.slen)
        fig.add_trace(
            go.Scatter(
                x=[(largest_start + largest_end) / 2],
                y=[(_CDS_BOX_Y_BOTTOM + 1.0) / 2],
                mode="text",
                text=[cds.product],
                textfont=dict(size=_CDS_LABEL_FONT_SIZE_PX),
                visible=True if label_visible else "legendonly",
                showlegend=False,
                hoverinfo="skip",
                meta=["cds_label", cds.product, cds.cds_id],
            ),
            row=row,
            col=1,
        )

    if coverage_label:
        y_domain_ref = "y domain" if row == 1 else f"y{row} domain"
        fig.add_annotation(
            text=coverage_label,
            x=track.slen * 0.99,
            y=1.0,
            xref=x_ref,
            yref=y_domain_ref,
            xanchor="right",
            yanchor="bottom",
            yshift=2,
            showarrow=False,
            font=dict(size=20),
        )


def render_plot_group(group: PlotGroup, *, title: str | None = None) -> go.Figure:
    """Render a single PlotGroup to a Plotly Figure."""
    n_tracks = len(group.tracks)

    if n_tracks == 0:
        fig = go.Figure()
        fig.update_layout(title="No tracks")
        return fig

    cds_legend_extra = (
        sum(len(t.cds_annotations) for t in group.tracks) * _CDS_LEGEND_ENTRY_HEIGHT_PX
    )
    plot_area_height = (
        n_tracks * _TRACK_CONTENT_HEIGHT_PX + max(0, n_tracks - 1) * _TRACK_SPACING_PX
    )
    figure_height = _MARGIN_TOP_PX + plot_area_height + _MARGIN_BOTTOM_BASE_PX + cds_legend_extra
    if n_tracks > 1:
        vertical_spacing = _TRACK_SPACING_PX / plot_area_height
    else:
        vertical_spacing = 0.0

    fig = make_subplots(
        rows=n_tracks,
        cols=1,
        shared_xaxes=False,
        vertical_spacing=vertical_spacing,
    )

    coverage_label = (
        f"coverage: {group.coverage_pct:.2f}% ({group.covered_bases}bp)"
        f" | best pident: {group.best_pident:.2f}"
    )

    for i, track in enumerate(group.tracks):
        row = i + 1
        label = _track_label(track.sseqid, track.segment)

        _add_track(
            fig,
            track,
            row,
            coverage_label=coverage_label if i == 0 else None,
        )

        x_range = [0, track.slen]
        is_bottom = row == n_tracks
        fig.update_xaxes(
            range=x_range,
            row=row,
            col=1,
            title_text="position on reference (sstart-send)" if is_bottom else None,
            title_font=dict(size=26) if is_bottom else None,
            tickfont=dict(size=24),
            linecolor="rgb(102,102,102)",
            linewidth=3,
            mirror=False,
            showgrid=False,
            showline=True,
        )
        fig.update_yaxes(
            range=[0, 1],
            row=row,
            col=1,
            tickvals=[(_PIDENT_RECT_Y_BOTTOM + _PIDENT_RECT_Y_TOP) / 2],
            ticktext=[label],
            tickfont=dict(size=26),
            linecolor="rgb(102,102,102)",
            linewidth=3,
            mirror=False,
            showgrid=False,
            showline=True,
        )

    if title is None:
        family = group.family or "unknown_family"
        genus = group.genus or "unknown_genus"
        title = f"{group.tsv_stem} | {family}, {genus} | {group.assembly}"

    margin_bottom = _MARGIN_BOTTOM_BASE_PX + cds_legend_extra

    fig.update_layout(
        title=dict(
            text=title,
            font=dict(size=28),
            x=0,
            xref="container",
            xanchor="left",
        ),
        height=figure_height,
        hovermode="closest",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(
            l=_MARGIN_LEFT_PX,
            r=_MARGIN_RIGHT_PX,
            t=_MARGIN_TOP_PX,
            b=margin_bottom,
        ),
        legend2=dict(
            orientation="h",
            y=-_CDS_LEGEND_OFFSET_PX / plot_area_height,
            yanchor="top",
            x=0.5,
            xanchor="center",
            font=dict(size=16),
            title=dict(text="CDS annotations", font=dict(size=18)),
        ),
    )

    _add_colorbar(fig, n_tracks=n_tracks)

    return fig


def _build_label_toggle_js(fig: go.Figure, div_id: str) -> str:
    """Build HTML buttons + per-CDS checkboxes to toggle CDS label traces."""
    import json
    from html import escape as html_escape

    cds_indices: dict[str, list[int]] = {}
    cds_display_names: dict[str, str] = {}
    auto_visible: dict[int, bool] = {}
    for i, t in enumerate(fig.data):
        meta = getattr(t, "meta", None)
        if isinstance(meta, (list, tuple)) and len(meta) >= 3 and meta[0] == "cds_label":
            product = str(meta[1])
            cds_id = str(meta[2])
            cds_indices.setdefault(cds_id, []).append(i)
            cds_display_names.setdefault(cds_id, product)
            auto_visible[i] = t.visible is True or t.visible is None

    if not cds_indices:
        return ""

    all_label_indices = sorted(auto_visible.keys())
    cds_ids = list(cds_indices.keys())

    cds_map_json = json.dumps(cds_indices)
    auto_vis_json = json.dumps({str(k): v for k, v in auto_visible.items()})
    all_indices_json = json.dumps(all_label_indices)
    cds_ids_json = json.dumps(cds_ids)

    checkboxes = []
    for idx, cds_id in enumerate(cds_ids):
        indices = cds_indices[cds_id]
        initially_checked = any(auto_visible.get(i, False) for i in indices)
        checked_attr = "checked" if initially_checked else ""
        label_html = html_escape(cds_display_names[cds_id])
        checkboxes.append(
            f'<label style="margin-right:8px;font-size:13px;cursor:pointer;">'
            f'<input type="checkbox" class="cds-cb-{div_id}" '
            f'data-idx="{idx}" {checked_attr} '
            f'onchange="toggleCds_{div_id}(this)"> '
            f"{label_html}</label>",
        )
    checkboxes_str = "\n    ".join(checkboxes)
    n_cds = len(cds_ids)

    js_template = _CDS_LABEL_TOGGLE_JS_PATH.read_text()
    const_block = (
        f'  const divId = "{div_id}";\n'
        f"  const cdsMap = {cds_map_json};\n"
        f"  const autoVis = {auto_vis_json};\n"
        f"  const allIndices = {all_indices_json};\n"
        f"  const cdsIds = {cds_ids_json};\n"
        f'  const cbClass = "cds-cb-{div_id}";'
    )
    start = js_template.index(_TEMPLATE_VARS_START_MARKER)
    end = js_template.index(_TEMPLATE_VARS_END_MARKER) + len(_TEMPLATE_VARS_END_MARKER)
    js_code = js_template[:start] + const_block + js_template[end:]

    return f"""
<div style="margin:8px 0 4px 60px;">
  <span style="font-size:13px;font-weight:bold;margin-right:8px;">CDS Labels:</span>
  <button onclick="restyleLabels_{div_id}('auto')" style="padding:4px 12px;font-size:13px;cursor:pointer;">Auto</button>
  <button onclick="restyleLabels_{div_id}('all')" style="padding:4px 12px;font-size:13px;cursor:pointer;">All</button>
  <button onclick="restyleLabels_{div_id}('off')" style="padding:4px 12px;font-size:13px;cursor:pointer;">Off</button>
</div>
<details open style="margin:0 0 8px 60px;">
  <summary style="font-size:13px;cursor:pointer;">CDS ({n_cds})</summary>
  <div style="margin-top:4px;">
    {checkboxes_str}
  </div>
</details>
<script>
{js_code}
</script>
"""


def _add_colorbar(fig: go.Figure, *, n_tracks: int) -> None:
    colorbar_len = 1.2 / n_tracks
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                colorscale=_PIDENT_COLORSCALE,
                cmin=_PIDENT_COLOR_MIN,
                cmax=_PIDENT_COLOR_MAX,
                colorbar=dict(
                    title=dict(text="pident", font=dict(size=28)),
                    tickfont=dict(size=24),
                    thickness=40,
                    len=colorbar_len,
                    x=1.02,
                    yanchor="top",
                    y=1.0,
                ),
                showscale=True,
            ),
            showlegend=False,
            hoverinfo="skip",
        ),
    )


def _sanitize_filename(text: str) -> str:
    return "".join(c if c.isalnum() or c in ("_", "-", ".") else "_" for c in text)


def write_html_figures(
    groups: list[PlotGroup],
    output_dir: Path,
    *,
    title: str | None = None,
) -> list[Path]:
    """Write one HTML file per PlotGroup. Returns list of written paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    seen_filenames: dict[str, str] = {}

    for group in groups:
        if not group.tracks:
            continue

        fig = render_plot_group(group, title=title)
        coverage_tag = f"{group.coverage_pct:.2f}"
        stem_safe = _sanitize_filename(group.tsv_stem)
        family_safe = _sanitize_filename(group.family or "unknown_family")
        assembly_safe = _sanitize_filename(group.assembly)
        filename = f"coverage_{coverage_tag}_{stem_safe}_{family_safe}_{assembly_safe}.html"

        old_key = f"coverage_{coverage_tag}_{family_safe}_{assembly_safe}.html"
        if old_key in seen_filenames:
            logger.warning(
                "Multiple samples hit the same family/assembly: "
                "tsv_stem=%s and %s both produce %s — "
                "resolved by including tsv_stem in filename",
                group.tsv_stem,
                seen_filenames[old_key],
                old_key,
            )
        seen_filenames.setdefault(old_key, group.tsv_stem)

        out_path = output_dir / filename

        div_id = filename.replace(".", "_").replace("-", "_")
        html_str = fig.to_html(
            include_plotlyjs=True,
            full_html=True,
            config={"toImageButtonOptions": {"format": "svg"}},
            div_id=div_id,
        )
        toggle_html = _build_label_toggle_js(fig, div_id)
        if toggle_html:
            div_tag = f'<div id="{div_id}"'
            html_str = html_str.replace(div_tag, toggle_html + div_tag)

        out_path.write_text(html_str)
        written.append(out_path)
        logger.info("Wrote %s", out_path)

    logger.info("Wrote %d HTML figures to %s", len(written), output_dir)
    return written
