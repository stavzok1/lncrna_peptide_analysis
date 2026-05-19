"""
Assemble main-text **combined** figures (Fig 1–6) at publication resolution.

Reads ``figures/combined_figure_manifest.yaml``, loads per-panel assets from
``figures/publication/`` (preferred) or ``figures/``, and writes **PDF + 300 dpi TIFF**
under ``figures/paper_tifs/`` (configurable).

External panels (Fig **2A** PowerPoint scheme, Fig **4B** lncBook) go under
``figures/external/`` as PNG or TIFF (300 dpi+). Marked ``optional: true`` in the manifest;
missing optional panels render a labeled placeholder so layout can be checked before
artwork is ready.

Prerequisites::

    python export_publication_figures.py
    # place exports:
    #   figures/external/fig2a_scheme.png
    #   figures/external/fig4b_lncbook.png

Usage::

    python manuscript/assemble_combined_manuscript_figures.py
    python manuscript/assemble_combined_manuscript_figures.py --only figure1 figure3
    python manuscript/assemble_combined_manuscript_figures.py --dpi 300 --width-in 11
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from repo_paths import FIGURES, FIGURES_PUBLICATION
from figure_export import COLOR_HALFTONE_TIFF_DPI, save_figure_bundle

import argparse

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import yaml
from matplotlib import image as mpimg
from PIL import Image

MANIFEST = FIGURES / "combined_figure_manifest.yaml"
EXTERNAL = FIGURES / "external"
PAPER_OUT = FIGURES / "paper_tifs"


def _load_rgba(path: Path) -> np.ndarray:
    im = Image.open(path)
    if im.mode != "RGBA":
        im = im.convert("RGBA")
    return np.asarray(im) / 255.0


def _trim_rgba_whitespace(
    arr: np.ndarray,
    *,
    rgb_thresh: float = 0.985,
    pad: int = 4,
) -> np.ndarray:
    """Crop near-white borders (common on exported lncBook / slide PNGs)."""
    if arr.ndim != 3 or arr.shape[2] < 3:
        return arr
    rgb = arr[:, :, :3]
    alpha = arr[:, :, 3] if arr.shape[2] >= 4 else np.ones((arr.shape[0], arr.shape[1]), dtype=arr.dtype)
    white = (
        (rgb[:, :, 0] > rgb_thresh)
        & (rgb[:, :, 1] > rgb_thresh)
        & (rgb[:, :, 2] > rgb_thresh)
        & (alpha > 0.5)
    )
    content = ~white
    if not np.any(content):
        return arr
    ys, xs = np.where(content)
    y0, y1 = int(ys.min()), int(ys.max()) + 1
    x0, x1 = int(xs.min()), int(xs.max()) + 1
    h, w = arr.shape[:2]
    y0 = max(0, y0 - pad)
    y1 = min(h, y1 + pad)
    x0 = max(0, x0 - pad)
    x1 = min(w, x1 + pad)
    if y1 - y0 < 16 or x1 - x0 < 16:
        return arr
    return arr[y0:y1, x0:x1].copy()


def _normalize_fontweight(weight: str | int) -> str | int:
    """Use numeric weights so 'normal' is not rendered as synthetic bold."""
    if isinstance(weight, int):
        return weight
    w = str(weight).strip().lower()
    if w in ("normal", "regular", "book", "", "none"):
        return 400
    if w in ("bold", "semibold", "demibold"):
        return 700
    if w in ("light", "thin"):
        return 300
    return weight


def resolve_panel_path(
    source: str,
    *,
    prefer: str,
    figures_root: Path,
    publication_root: Path,
) -> Path | None:
    """Resolve stem under publication/, figures/, or figures/external/."""
    rel = Path(source)
    if rel.is_absolute() and rel.exists():
        return rel
    stem = rel.name
    parent = rel.parent
    candidates: list[Path] = []
    if str(parent).startswith("external"):
        for ext in (".png", ".tif", ".tiff", ".pdf"):
            candidates.append(figures_root / parent / f"{stem}{ext}")
    if prefer == "publication":
        for ext in (".tif", ".tiff", ".png"):
            candidates.append(publication_root / parent / f"{stem}{ext}")
            if parent == Path("."):
                candidates.append(publication_root / f"{stem}{ext}")
    for ext in (".png", ".tif", ".tiff"):
        candidates.append(figures_root / parent / f"{stem}{ext}")
        if parent == Path("."):
            candidates.append(figures_root / f"{stem}{ext}")
    for p in candidates:
        if p.exists():
            return p
    return None


def _placeholder_rgba(label: str, note: str) -> np.ndarray:
    h, w = 480, 960
    arr = np.ones((h, w, 4), dtype=np.float64)
    arr[:, :, :3] *= 0.94
    arr[:, :, 3] = 1.0
    return arr


def _panel_load(
    panel: dict,
    *,
    prefer: str,
    figures_root: Path,
    publication_root: Path,
) -> tuple[np.ndarray, float]:
    """Return RGBA array and width/height aspect ratio."""
    source = str(panel["source"])
    optional = bool(panel.get("optional", False))
    path = resolve_panel_path(
        source, prefer=prefer, figures_root=figures_root, publication_root=publication_root
    )
    if path is None:
        if optional:
            arr = _placeholder_rgba("", source)
            return arr, arr.shape[1] / arr.shape[0]
        raise FileNotFoundError(
            f"Missing panel source '{source}'. Run export_publication_figures.py or check manifest."
        )
    arr = _load_rgba(path)
    if bool(panel.get("trim_whitespace")):
        arr_t = _trim_rgba_whitespace(arr)
        if arr_t.shape[0] >= 16 and arr_t.shape[1] >= 16:
            arr = arr_t
    h, w = arr.shape[0], arr.shape[1]
    return arr, w / h


def _panel_aspect(
    panel: dict,
    *,
    prefer: str,
    figures_root: Path,
    publication_root: Path,
) -> float:
    _, aspect = _panel_load(
        panel, prefer=prefer, figures_root=figures_root, publication_root=publication_root
    )
    return aspect


LABEL_PAD_X = 0.012
LABEL_PAD_Y = 0.010
# Panel letters above axes top edge (transAxes); va='bottom' keeps ink outside the raster.
OUTSIDE_LABEL_TX = -0.015
OUTSIDE_LABEL_TY = 1.0


def _layout_float(layout: dict, defaults: dict, key: str, fallback: float) -> float:
    if key in layout and layout[key] is not None:
        return float(layout[key])
    if key in defaults and defaults[key] is not None:
        return float(defaults[key])
    return fallback


@dataclass
class PanelLabelSpec:
    ax: plt.Axes
    label: str
    col: int | None = None
    row: int | None = None
    position: str = "outside_top"


def _panel_image_anchor(panel: dict, *, fill_cells: bool) -> str:
    """Matplotlib axes anchor for imshow (N = top-align panels in a shared grid row)."""
    raw = str(panel.get("anchor", "")).strip().lower()
    if raw in ("n", "north", "top"):
        return "N"
    if raw in ("c", "center", "centre", ""):
        return "N" if fill_cells else "C"
    return raw.upper()[:1] if raw else ("N" if fill_cells else "C")


def _place_all_panel_labels(
    fig: plt.Figure,
    specs: list[PanelLabelSpec],
    *,
    fontsize: float,
    weight: str,
    pad_x: float = LABEL_PAD_X,
    pad_y: float = LABEL_PAD_Y,
    unify_outside_x: bool = False,
) -> None:
    """Place panel letters outside each axes (transAxes), not over the raster."""
    if not specs:
        return
    fig.canvas.draw()

    fw = _normalize_fontweight(weight)
    x_fig_unified: float | None = None
    if unify_outside_x:
        outside = [sp for sp in specs if sp.label and sp.position != "inside_top"]
        if outside:
            x0_m = min(sp.ax.get_position().x0 for sp in outside)
            x_fig_unified = max(0.0, x0_m - pad_x)

    for sp in specs:
        if not sp.label:
            continue
        if sp.position == "inside_top":
            sp.ax.text(
                0.02,
                0.98,
                sp.label,
                transform=sp.ax.transAxes,
                ha="left",
                va="top",
                fontsize=fontsize,
                fontweight=fw,
                color="black",
            )
            continue
        if x_fig_unified is not None:
            pos = sp.ax.get_position()
            y_fig = min(1.0, pos.y1 + pad_y)
            fig.text(
                x_fig_unified,
                y_fig,
                sp.label,
                transform=fig.transFigure,
                ha="left",
                va="bottom",
                fontsize=fontsize,
                fontweight=fw,
                color="black",
                clip_on=False,
            )
            continue
        sp.ax.text(
            OUTSIDE_LABEL_TX,
            OUTSIDE_LABEL_TY,
            sp.label,
            transform=sp.ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=fontsize,
            fontweight=fw,
            color="black",
            clip_on=False,
        )


def _draw_panel(
    ax,
    panel: dict,
    *,
    prefer: str,
    figures_root: Path,
    publication_root: Path,
    label_fontsize: float,
    label_weight: str,
    label_specs: list[PanelLabelSpec],
    fill_cells: bool = False,
) -> None:
    ax.set_axis_off()
    source = str(panel["source"])
    label = str(panel.get("label", ""))
    optional = bool(panel.get("optional", False))
    path = resolve_panel_path(
        source, prefer=prefer, figures_root=figures_root, publication_root=publication_root
    )
    if path is None and optional:
        arr = _placeholder_rgba(label, source)
        anchor = _panel_image_anchor(panel, fill_cells=fill_cells)
        imshow_kw = dict(interpolation="nearest")
        aspect = "auto" if fill_cells else "equal"
        ax.imshow(arr, aspect=aspect, **imshow_kw)
        ax.set_anchor(anchor)
        ax.text(
            0.5,
            0.5,
            f"External panel\n{source}\n(optional — add under figures/external/)",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
            color="#333333",
        )
    else:
        arr, _ = _panel_load(
            panel, prefer=prefer, figures_root=figures_root, publication_root=publication_root
        )
        anchor = _panel_image_anchor(panel, fill_cells=fill_cells)
        imshow_kw = dict(interpolation="nearest")
        aspect = "auto" if fill_cells else "equal"
        ax.imshow(arr, aspect=aspect, **imshow_kw)
        ax.set_anchor(anchor)
    col_raw = panel.get("col")
    col = int(col_raw) if col_raw is not None else None
    row_raw = panel.get("row")
    row = int(row_raw) if row_raw is not None else None
    label_specs.append(
        PanelLabelSpec(
            ax=ax,
            label=label,
            col=col,
            row=row,
            position=str(panel.get("label_position", "outside_top")),
        )
    )


def _horizontal_figsize(
    panels: list[dict],
    width_in: float,
    *,
    prefer: str,
    width_ratios: list[float] | None,
    gutter: float = 0.035,
) -> tuple[float, float, list[float]]:
    aspects = [
        _panel_aspect(p, prefer=prefer, figures_root=FIGURES, publication_root=FIGURES_PUBLICATION)
        for p in panels
    ]
    wr = width_ratios if width_ratios is not None else aspects
    total_wr = sum(wr)
    row_h = max(width_in * w / (total_wr * a) for w, a in zip(wr, aspects))
    height_in = row_h * (1.0 + gutter)
    return width_in, height_in, wr


def _grid_figsize(
    panels: list[dict],
    nrows: int,
    ncols: int,
    width_in: float,
    *,
    prefer: str,
    width_ratios: list[float] | None,
    height_ratios: list[float] | None,
    gutter: float = 0.045,
    apply_gutter: bool = True,
) -> tuple[float, float, list[float], list[float]]:
    wr = list(width_ratios) if width_ratios is not None else [1.0] * ncols
    hr = list(height_ratios) if height_ratios is not None else [1.0] * nrows
    wr_sum = sum(wr)
    hr_sum = sum(hr)
    col_widths = [width_in * w / wr_sum for w in wr]
    row_heights = [hr[i] / hr_sum * width_in * 0.9 / nrows for i in range(nrows)]

    for p in panels:
        r = int(p["row"])
        c = int(p["col"])
        rs = int(p.get("rowspan", 1))
        cs = int(p.get("colspan", 1))
        cw = sum(col_widths[c : c + cs])
        a = _panel_aspect(
            p, prefer=prefer, figures_root=FIGURES, publication_root=FIGURES_PUBLICATION
        )
        need_h = cw / a
        have_h = sum(row_heights[r : r + rs])
        if have_h < need_h and have_h > 0:
            scale = need_h / have_h
            for i in range(r, r + rs):
                row_heights[i] *= scale

    height_in = sum(row_heights) * ((1.0 + gutter) if apply_gutter else 1.0)
    return width_in, height_in, wr, row_heights


def _panel_col_width_in(
    panel: dict,
    col_widths_in: list[float],
) -> float:
    col = int(panel["col"])
    colspan = int(panel.get("colspan", 1))
    return sum(col_widths_in[col : col + colspan])


def _panel_height_at_width(
    panel: dict,
    width_in: float,
    *,
    prefer: str,
    figures_root: Path,
    publication_root: Path,
) -> float:
    aspect = _panel_aspect(
        panel, prefer=prefer, figures_root=figures_root, publication_root=publication_root
    )
    return width_in / aspect if aspect > 0 else width_in


def _grid_figsize_split(
    panels: list[dict],
    nrows: int,
    ncols: int,
    width_in: float,
    *,
    prefer: str,
    split_row: int,
    top_width_ratios: list[float],
    bottom_width_ratios: list[float],
    height_ratios: list[float] | None,
    gutter: float = 0.04,
) -> tuple[float, float, list[float], list[float], list[float], list[float], list[float]]:
    """Top/bottom sections use different column widths; heights fit image content (no letterbox gap)."""
    top_wr = list(top_width_ratios)
    bot_wr = list(bottom_width_ratios)
    top_cw = [width_in * w / sum(top_wr) for w in top_wr]
    bot_cw = [width_in * w / sum(bot_wr) for w in bot_wr]

    top_panels = [p for p in panels if int(p["row"]) < split_row]
    bot_panels = [p for p in panels if int(p["row"]) >= split_row]

    def _heights(ps: list[dict], col_widths: list[float]) -> list[float]:
        return [
            _panel_height_at_width(
                p,
                _panel_col_width_in(p, col_widths),
                prefer=prefer,
                figures_root=FIGURES,
                publication_root=FIGURES_PUBLICATION,
            )
            for p in ps
        ]

    row_heights_top: list[float] = []
    for r in range(split_row):
        row_ps = [p for p in top_panels if int(p["row"]) == r and int(p.get("rowspan", 1)) == 1]
        if row_ps:
            row_heights_top.append(max(_heights(row_ps, top_cw)))
        else:
            row_heights_top.append(0.0)
    if not row_heights_top or max(row_heights_top) <= 0:
        row_heights_top = [1.0] * split_row

    span_ps = [p for p in top_panels if int(p.get("rowspan", 1)) >= split_row]
    span_h = 0.0
    if span_ps:
        span_h = max(
            _panel_height_at_width(
                p,
                _panel_col_width_in(p, top_cw),
                prefer=prefer,
                figures_root=FIGURES,
                publication_root=FIGURES_PUBLICATION,
            )
            for p in span_ps
        )
    top_section_h = max(span_h, sum(row_heights_top))

    bot_row_ps = [p for p in bot_panels if int(p.get("rowspan", 1)) == 1]
    bot_hs = _heights(bot_row_ps, bot_cw) if bot_row_ps else [1.0]
    bot_section_h = max(bot_hs) if bot_hs else 1.0

    # Inner top row weights: B and C stack
    if split_row == 2 and row_heights_top[0] > 0 and row_heights_top[1] > 0:
        top_rh = [row_heights_top[0], row_heights_top[1]]
    else:
        top_rh = [top_section_h / split_row] * split_row
    rh_sum = sum(top_rh)
    if rh_sum > 0 and top_section_h > rh_sum:
        scale = top_section_h / rh_sum
        top_rh = [r * scale for r in top_rh]

    bot_rh = [bot_section_h]
    outer_top = top_section_h
    outer_bot = bot_section_h
    height_in = (outer_top + outer_bot) * (1.0 + gutter)
    return width_in, height_in, top_wr, top_rh, bot_wr, bot_rh, [outer_top, outer_bot]


def _vertical_figsize(
    panels: list[dict],
    width_in: float,
    *,
    prefer: str,
    height_ratios: list[float] | None,
    gutter: float = 0.035,
) -> tuple[float, float]:
    aspects = [
        _panel_aspect(p, prefer=prefer, figures_root=FIGURES, publication_root=FIGURES_PUBLICATION)
        for p in panels
    ]
    if height_ratios is not None:
        heights = [width_in / a * h for a, h in zip(aspects, height_ratios)]
    else:
        heights = [width_in / a for a in aspects]
    height_in = sum(heights) * (1.0 + gutter)
    return width_in, height_in


def _build_figure(spec: dict, defaults: dict) -> plt.Figure:
    layout = spec["layout"]
    ltype = layout["type"]
    width_in = float(spec.get("width_in", defaults.get("width_in", 11.0)))
    height_in = spec.get("height_in")
    dpi = int(spec.get("dpi", defaults.get("dpi", COLOR_HALFTONE_TIFF_DPI)))
    prefer = str(spec.get("prefer", defaults.get("prefer", "publication")))
    label_fontsize = float(spec.get("label_fontsize", defaults.get("label_fontsize", 16)))
    label_weight = str(spec.get("label_weight", defaults.get("label_weight", "normal")))

    label_specs: list[PanelLabelSpec] = []

    def _finish(
        fig: plt.Figure,
        *,
        grid_margins: bool = False,
        multipanel: bool = False,
        unify_outside_label_x: bool = False,
    ) -> plt.Figure:
        has_label = any(bool(sp.label) for sp in label_specs)
        if grid_margins:
            fig.subplots_adjust(left=0.072, right=0.994, top=0.982, bottom=0.016)
        elif multipanel and has_label:
            fig.subplots_adjust(left=0.056, right=0.994, top=0.983, bottom=0.010)
        elif multipanel:
            fig.subplots_adjust(left=0.056, right=0.994, top=0.988, bottom=0.010)
        _place_all_panel_labels(
            fig,
            label_specs,
            fontsize=label_fontsize,
            weight=label_weight,
            unify_outside_x=unify_outside_label_x,
        )
        return fig

    if ltype == "single":
        panels = layout["panels"]
        aspect = _panel_aspect(
            panels[0],
            prefer=prefer,
            figures_root=FIGURES,
            publication_root=FIGURES_PUBLICATION,
        )
        fig_h = float(height_in) if height_in is not None else width_in / aspect
        fig = plt.figure(figsize=(width_in, fig_h), dpi=dpi)
        ax = fig.add_subplot(111)
        _draw_panel(
            ax,
            panels[0],
            prefer=prefer,
            figures_root=FIGURES,
            publication_root=FIGURES_PUBLICATION,
            label_fontsize=label_fontsize,
            label_weight=label_weight,
            label_specs=label_specs,
        )
        return _finish(fig, multipanel=False)

    if ltype == "vertical":
        panels = layout["panels"]
        fgutter = _layout_float(layout, defaults, "figure_gutter", 0.035)
        v_hspace = _layout_float(layout, defaults, "vertical_hspace", 0.008)
        w_in, h_in = _vertical_figsize(
            panels,
            width_in,
            prefer=prefer,
            height_ratios=layout.get("height_ratios"),
            gutter=fgutter,
        )
        if height_in is not None:
            h_in = float(height_in)
        fig = plt.figure(figsize=(w_in, h_in), dpi=dpi)
        gs = gridspec.GridSpec(len(panels), 1, hspace=v_hspace)
        for i, p in enumerate(panels):
            ax = fig.add_subplot(gs[i, 0])
            p = dict(p)
            p.setdefault("col", 0)
            _draw_panel(
                ax,
                p,
                prefer=prefer,
                figures_root=FIGURES,
                publication_root=FIGURES_PUBLICATION,
                label_fontsize=label_fontsize,
                label_weight=label_weight,
                label_specs=label_specs,
            )
        return _finish(fig, multipanel=True, unify_outside_label_x=True)

    if ltype == "horizontal":
        panels = layout["panels"]
        wr_layout = layout.get("width_ratios")
        fgutter = _layout_float(layout, defaults, "figure_gutter", 0.035)
        h_wspace = _layout_float(layout, defaults, "horizontal_wspace", 0.025)
        w_in, h_in, wr = _horizontal_figsize(
            panels,
            width_in,
            prefer=prefer,
            width_ratios=wr_layout,
            gutter=fgutter,
        )
        if height_in is not None:
            h_in = float(height_in)
        fig = plt.figure(figsize=(w_in, h_in), dpi=dpi)
        gs = gridspec.GridSpec(1, len(panels), width_ratios=wr, wspace=h_wspace)
        for i, p in enumerate(panels):
            ax = fig.add_subplot(gs[0, i])
            p = dict(p)
            p.setdefault("col", i)
            _draw_panel(
                ax,
                p,
                prefer=prefer,
                figures_root=FIGURES,
                publication_root=FIGURES_PUBLICATION,
                label_fontsize=label_fontsize,
                label_weight=label_weight,
                label_specs=label_specs,
            )
        return _finish(fig, multipanel=True)

    if ltype == "grid":
        panels = layout["panels"]
        nrows = int(layout["nrows"])
        ncols = int(layout["ncols"])
        split_row = layout.get("split_row")
        top_wr_cfg = layout.get("top_width_ratios")
        bot_wr_cfg = layout.get("bottom_width_ratios")
        use_split = (
            split_row is not None
            and top_wr_cfg is not None
            and bot_wr_cfg is not None
        )

        if use_split:
            split_row = int(split_row)
            fgutter = _layout_float(layout, defaults, "figure_gutter", 0.04)
            split_outer = _layout_float(layout, defaults, "split_outer_hspace", 0.012)
            split_w = _layout_float(layout, defaults, "split_inner_wspace", 0.035)
            split_h = _layout_float(layout, defaults, "split_inner_hspace", 0.032)
            w_in, h_in, top_wr, top_rh, bot_wr, bot_rh, outer_h = _grid_figsize_split(
                panels,
                nrows,
                ncols,
                width_in,
                prefer=prefer,
                split_row=split_row,
                top_width_ratios=list(top_wr_cfg),
                bottom_width_ratios=list(bot_wr_cfg),
                height_ratios=layout.get("height_ratios"),
                gutter=fgutter,
            )
            fill_cells = bool(layout.get("fill_cells", True))
        else:
            fgutter = _layout_float(layout, defaults, "figure_gutter", 0.045)
            w_in, h_in, col_wr, row_hr = _grid_figsize(
                panels,
                nrows,
                ncols,
                width_in,
                prefer=prefer,
                width_ratios=layout.get("width_ratios"),
                height_ratios=layout.get("height_ratios"),
                gutter=fgutter,
            )
            gw = _layout_float(layout, defaults, "grid_wspace", 0.038)
            gh = _layout_float(layout, defaults, "grid_hspace", 0.032)
        if height_in is not None:
            h_in = float(height_in)
        fig = plt.figure(figsize=(w_in, h_in), dpi=dpi)

        if use_split:
            gs_outer = gridspec.GridSpec(
                2,
                1,
                figure=fig,
                height_ratios=outer_h,
                hspace=split_outer,
            )
            gs_top = gs_outer[0].subgridspec(
                split_row,
                ncols,
                width_ratios=top_wr,
                height_ratios=top_rh,
                wspace=split_w,
                hspace=split_h,
            )
            gs_bot = gs_outer[1].subgridspec(
                nrows - split_row,
                ncols,
                width_ratios=bot_wr,
                height_ratios=bot_rh,
                wspace=split_w,
                hspace=split_h,
            )
            for p in panels:
                row, col = int(p["row"]), int(p["col"])
                colspan = int(p.get("colspan", 1))
                rowspan = int(p.get("rowspan", 1))
                if row < split_row:
                    ax = fig.add_subplot(gs_top[row : row + rowspan, col : col + colspan])
                else:
                    local_row = row - split_row
                    ax = fig.add_subplot(
                        gs_bot[local_row : local_row + rowspan, col : col + colspan]
                    )
                _draw_panel(
                    ax,
                    p,
                    prefer=prefer,
                    figures_root=FIGURES,
                    publication_root=FIGURES_PUBLICATION,
                    label_fontsize=label_fontsize,
                    label_weight=label_weight,
                    label_specs=label_specs,
                    fill_cells=fill_cells,
                )
        else:
            gs = gridspec.GridSpec(
                nrows,
                ncols,
                figure=fig,
                width_ratios=col_wr,
                height_ratios=row_hr,
                wspace=gw,
                hspace=gh,
            )
            for p in panels:
                row, col = int(p["row"]), int(p["col"])
                colspan = int(p.get("colspan", 1))
                rowspan = int(p.get("rowspan", 1))
                ax = fig.add_subplot(gs[row : row + rowspan, col : col + colspan])
                _draw_panel(
                    ax,
                    p,
                    prefer=prefer,
                    figures_root=FIGURES,
                    publication_root=FIGURES_PUBLICATION,
                    label_fontsize=label_fontsize,
                    label_weight=label_weight,
                    label_specs=label_specs,
                )
        return _finish(fig, grid_margins=True)

    raise ValueError(f"Unknown layout type: {ltype}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--manifest",
        type=Path,
        default=MANIFEST,
        help=f"YAML layout file (default: {MANIFEST.relative_to(_REPO)}).",
    )
    ap.add_argument(
        "--only",
        nargs="*",
        metavar="FIG",
        help="Subset keys from manifest (e.g. figure1 figure4). Default: all.",
    )
    ap.add_argument("--output-dir", type=Path, default=PAPER_OUT)
    ap.add_argument("--width-in", type=float, default=None, help="Override default combined width (inches).")
    ap.add_argument("--dpi", type=int, default=None, help="Raster DPI for TIFF (default 300).")
    args = ap.parse_args()

    if not args.manifest.exists():
        raise SystemExit(f"Missing manifest: {args.manifest}")

    cfg = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    defaults = dict(cfg.get("defaults", {}))
    if args.width_in is not None:
        defaults["width_in"] = args.width_in
    if args.dpi is not None:
        defaults["dpi"] = args.dpi

    out_dir = args.output_dir
    if not out_dir.is_absolute():
        out_dir = FIGURES / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    EXTERNAL.mkdir(parents=True, exist_ok=True)

    figures: dict = cfg.get("figures", {})
    keys = sorted(figures.keys())
    if args.only:
        keys = [k for k in keys if k in set(args.only)]
        if not keys:
            raise SystemExit(f"No manifest keys matched --only {args.only}")

    for key in keys:
        spec = figures[key]
        stem = spec.get("output_stem", key)
        print(f"Assembling {key} -> {stem} ...")
        fig = _build_figure(spec, defaults)
        out_png = out_dir / f"{stem}.png"
        save_figure_bundle(
            fig,
            out_png,
            png_dpi=defaults.get("dpi", COLOR_HALFTONE_TIFF_DPI),
            publication_dir=None,
            bbox_inches="tight",
            pad_inches=0.028,
            facecolor="white",
        )
        stem_path = out_png.with_suffix("")
        pdf_kw = dict(bbox_inches="tight", pad_inches=0.028, facecolor="white")
        fig.savefig(stem_path.with_suffix(".pdf"), format="pdf", **pdf_kw)
        tif_dpi = int(defaults.get("dpi", COLOR_HALFTONE_TIFF_DPI))
        fig.savefig(
            stem_path.with_suffix(".tif"),
            format="tiff",
            dpi=tif_dpi,
            pil_kwargs={"compression": "tiff_lzw"},
            **pdf_kw,
        )
        plt.close(fig)
        print(f"  Wrote {stem_path}.{{png,pdf,tif}} @ {tif_dpi} dpi")


if __name__ == "__main__":
    main()
