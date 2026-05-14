"""
Publication exports: vector **PDF** plus **TIFF** raster at journal-style DPI.

- **PDF:** vector (matplotlib does not use ``dpi`` for vector output).
- **TIFF:** LZW compression. ``--publication-tiff-kind line`` → **1200 dpi** (line art).
  ``color`` → **300 dpi** (colour / halftone-style figures). Default is ``color``, which
  matches common journal guidance for coloured scatter / heatmap panels.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

LINE_ART_TIFF_DPI = 1200
COLOR_HALFTONE_TIFF_DPI = 300


def add_publication_args(ap: Any) -> None:
    ap.add_argument(
        "--publication-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Also write PDF + TIFF next to each PNG path, under DIR mirroring the path "
        "relative to figures/ when the PNG lives under figures/; otherwise flat under DIR.",
    )
    ap.add_argument(
        "--publication-tiff-kind",
        choices=("line", "color"),
        default="color",
        help="TIFF rasterization: line=1200 dpi; color=300 dpi (default: color).",
    )


def publication_output_dir(
    png_path: Path,
    publication_root: Path,
    *,
    figures_root: Path,
) -> Path:
    """Directory under ``publication_root`` that mirrors ``png_path``'s parent relative to ``figures_root``."""
    publication_root = Path(publication_root)
    png_path = Path(png_path)
    figures_root = Path(figures_root)
    try:
        rel_parent = png_path.resolve().relative_to(figures_root.resolve()).parent
    except ValueError:
        rel_parent = Path()
    out_dir = publication_root / rel_parent
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def save_figure_bundle(
    fig: Any,
    png_path: Path,
    *,
    png_dpi: float | None,
    publication_dir: Path | None = None,
    publication_tiff_kind: str = "color",
    figures_root: Path | None = None,
    **kwargs: Any,
) -> None:
    """
    Write primary PNG, then optional PDF + TIFF under ``publication_dir`` (when set).

    ``publication_dir`` is the **root** for mirrored layout; pass ``figures_root`` so paths
    under ``figures/`` preserve subfolders (e.g. ``tcga_matrix/peptide_fraction/``).
    """
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    png_kw = dict(kwargs)
    if png_dpi is not None:
        png_kw["dpi"] = png_dpi
    fig.savefig(png_path, **png_kw)
    if publication_dir is None:
        return
    if figures_root is None:
        raise ValueError("save_figure_bundle requires figures_root when publication_dir is set")
    pub_leaf = publication_output_dir(png_path, publication_dir, figures_root=figures_root)
    stem = pub_leaf / png_path.stem
    pdf_kw = {k: v for k, v in png_kw.items() if k != "dpi"}
    fig.savefig(stem.with_suffix(".pdf"), format="pdf", **pdf_kw)
    tif_dpi = LINE_ART_TIFF_DPI if publication_tiff_kind == "line" else COLOR_HALFTONE_TIFF_DPI
    tif_kw = {k: v for k, v in pdf_kw.items() if k not in ("format", "pil_kwargs")}
    fig.savefig(
        stem.with_suffix(".tif"),
        format="tiff",
        dpi=tif_dpi,
        pil_kwargs={"compression": "tiff_lzw"},
        **tif_kw,
    )
