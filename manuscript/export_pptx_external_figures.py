"""
Export external manuscript panels from a PowerPoint file at publication resolution.

Primary art source (pick one):

- **Combined PDF** (default if present): ``figures/ECCB-Figures OPEN to STAV (2a and 4b).pdf``
  — **page 1** → Fig 4B ``figures/external/fig4b_lncbook.*``
  — **page 2** → full Fig 2 ``figures/external/fig2_full.*`` (no crop)

- **PowerPoint** (fallback): ``figures/ECCB-Figures OPEN to STAV (2a and 4b).pptx``
  — same mapping as slides 1–2.

Optional: ``--fig2-pdf path`` = a **single-page** PDF for Fig 2 only (uses PPTX for 4B).

Usage::

    python manuscript/export_pptx_external_figures.py
    python manuscript/export_pptx_external_figures.py --eccb-pdf figures/your.pdf
    python manuscript/export_pptx_external_figures.py --fig2-pdf figures/external/fig2.pdf
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
import zipfile
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from repo_paths import FIGURES
from figure_export import COLOR_HALFTONE_TIFF_DPI

DEFAULT_PPTX = FIGURES / "ECCB-Figures OPEN to STAV (2a and 4b).pptx"
DEFAULT_ECCB_PDF = FIGURES / "ECCB-Figures OPEN to STAV (2a and 4b).pdf"
EXTERNAL = FIGURES / "external"


def _slide_size_inches(pptx: Path) -> tuple[float, float]:
    with zipfile.ZipFile(pptx) as z:
        xml = z.read("ppt/presentation.xml").decode("utf-8")
    m = re.search(r'sldSz[^>]+cx="(\d+)"[^>]+cy="(\d+)"', xml)
    if not m:
        return 13.333, 7.5
    cx, cy = int(m.group(1)), int(m.group(2))
    return cx / 914400, cy / 914400


def _export_slides_com(pptx: Path, out_dir: Path, dpi: int) -> list[Path]:
    import win32com.client  # type: ignore

    pptx = pptx.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    w_in, _h_in = _slide_size_inches(pptx)
    width_px = int(round(w_in * dpi))
    height_px = int(round(_h_in * dpi))

    app = win32com.client.Dispatch("PowerPoint.Application")
    try:
        app.Visible = 1
    except Exception:
        pass
    paths: list[Path] = []
    try:
        pres = app.Presentations.Open(str(pptx), WithWindow=False, ReadOnly=True)
        try:
            for i in range(1, pres.Slides.Count + 1):
                tmp = out_dir / f"_slide{i}_export.png"
                pres.Slides(i).Export(str(tmp), "PNG", width_px, height_px)
                paths.append(tmp)
        finally:
            pres.Close()
    finally:
        app.Quit()
    return paths


def _extract_embedded(pptx: Path, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(pptx) as z:
        media = sorted(
            n for n in z.namelist() if n.startswith("ppt/media/") and n.lower().endswith(".png")
        )
        paths: list[Path] = []
        for i, name in enumerate(media, start=1):
            dst = out_dir / f"_slide{i}_export.png"
            dst.write_bytes(z.read(name))
            paths.append(dst)
        return paths


def _require_fitz():
    try:
        import fitz  # PyMuPDF  # noqa: F401
    except ImportError as exc:
        raise SystemExit("Install PyMuPDF for PDF input: pip install pymupdf") from exc
    import fitz

    return fitz


def _rasterize_pdf_page(pdf_path: Path, page_zero_based: int, out_png: Path, dpi: int) -> None:
    fitz = _require_fitz()
    doc = fitz.open(pdf_path)
    try:
        if page_zero_based < 0 or page_zero_based >= doc.page_count:
            raise SystemExit(
                f"{pdf_path}: page {page_zero_based + 1} out of range (doc has {doc.page_count} page(s))"
            )
        page = doc[page_zero_based]
        zoom = dpi / 72.0
        mat = fitz.Matrix(zoom, zoom)
        pix = page.get_pixmap(matrix=mat, alpha=False)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        pix.save(str(out_png))
    finally:
        doc.close()


def _pdf_page_count(pdf_path: Path) -> int:
    fitz = _require_fitz()
    doc = fitz.open(pdf_path)
    try:
        return doc.page_count
    finally:
        doc.close()


def _rasterize_pdf(pdf_path: Path, out_png: Path, dpi: int) -> None:
    _rasterize_pdf_page(pdf_path, 0, out_png, dpi)


def _write_bundle(stem: Path, src_png: Path, *, publication_dir: Path | None, dpi: int) -> None:
    from PIL import Image

    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_png, stem.with_suffix(".png"))

    im = Image.open(stem.with_suffix(".png"))
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        im = bg
    elif im.mode != "RGB":
        im = im.convert("RGB")

    im.save(stem.with_suffix(".pdf"), "PDF", resolution=float(dpi))
    im.save(stem.with_suffix(".tif"), format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))

    if publication_dir is not None:
        publication_dir = Path(publication_dir)
        publication_dir.mkdir(parents=True, exist_ok=True)
        for ext in (".png", ".pdf", ".tif"):
            shutil.copy2(stem.with_suffix(ext), publication_dir / stem.with_suffix(ext).name)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pptx", type=Path, default=DEFAULT_PPTX)
    ap.add_argument("--dpi", type=int, default=COLOR_HALFTONE_TIFF_DPI)
    try:
        rel = DEFAULT_ECCB_PDF.relative_to(_REPO)
    except ValueError:
        rel = DEFAULT_ECCB_PDF
    ap.add_argument(
        "--eccb-pdf",
        type=Path,
        default=None,
        help="Two-page deck: page 1 → fig4b_lncbook, page 2 → fig2_full. Default: "
        f"{rel} if that file exists.",
    )
    ap.add_argument(
        "--fig2-pdf",
        type=Path,
        default=None,
        help="Single-page PDF for fig2_full only (Fig 4B still from PPTX slide 1).",
    )
    ap.add_argument("--output-dir", type=Path, default=EXTERNAL)
    ap.add_argument("--publication-dir", type=Path, default=FIGURES / "publication" / "external")
    args = ap.parse_args()

    out = args.output_dir
    pub = args.publication_dir
    tmp = out / "_ppt_export_tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)

    fig2_png = tmp / "fig2_full_src.png"
    pptx = args.pptx if args.pptx.is_absolute() else _REPO / args.pptx
    slide_pngs: list[Path] = []

    eccb = args.eccb_pdf
    if eccb is not None:
        eccb = eccb if eccb.is_absolute() else _REPO / eccb
    elif DEFAULT_ECCB_PDF.is_file():
        eccb = DEFAULT_ECCB_PDF

    if eccb is not None and eccb.is_file():
        _require_fitz()
        n = _pdf_page_count(eccb)
        if n < 2:
            raise SystemExit(f"{eccb} must have at least 2 pages (got {n}).")
        fig4_png = tmp / "eccb_page1.png"
        _rasterize_pdf_page(eccb, 0, fig4_png, args.dpi)
        fig4b = out / "fig4b_lncbook"
        _write_bundle(fig4b, fig4_png, publication_dir=pub, dpi=args.dpi)
        print(f"Fig 4B from PDF page 1: {eccb} @ {args.dpi} dpi")
        _rasterize_pdf_page(eccb, 1, fig2_png, args.dpi)
        print(f"Fig 2 full from PDF page 2: {eccb} @ {args.dpi} dpi")
    else:
        if pptx.is_file():
            try:
                slide_pngs = _export_slides_com(pptx, tmp, args.dpi)
                print(f"Exported {len(slide_pngs)} slide(s) via PowerPoint COM @ {args.dpi} dpi")
            except Exception as exc:
                print(f"PowerPoint COM failed ({exc!r}); using embedded images (lower resolution).")
                slide_pngs = _extract_embedded(pptx, tmp)
            if slide_pngs:
                fig4b = out / "fig4b_lncbook"
                _write_bundle(fig4b, slide_pngs[0], publication_dir=pub, dpi=args.dpi)
                print(f"Wrote {fig4b}.png (+ pdf/tif)")

        fig2_pdf = args.fig2_pdf
        if fig2_pdf is not None:
            fig2_pdf = fig2_pdf if fig2_pdf.is_absolute() else _REPO / fig2_pdf
        else:
            fig2_pdf = None
            for candidate in (out / "fig2.pdf", out / "fig2_full.pdf"):
                if candidate.is_file():
                    fig2_pdf = candidate
                    break

        if fig2_pdf is not None and fig2_pdf.is_file():
            _rasterize_pdf(fig2_pdf, fig2_png, args.dpi)
            print(f"Rasterized Fig 2 from PDF: {fig2_pdf} @ {args.dpi} dpi")
        elif len(slide_pngs) >= 2:
            shutil.copy2(slide_pngs[1], fig2_png)
            print("Fig 2 from PPTX slide 2 (full slide, no crop)")
        else:
            raise SystemExit(
                "Fig 2 / Fig 4B source missing: add ECCB PDF, or PPTX with 2 slides, or --fig2-pdf."
            )

    fig2 = out / "fig2_full"
    _write_bundle(fig2, fig2_png, publication_dir=pub, dpi=args.dpi)
    print(f"Wrote {fig2}.png (+ pdf/tif) — full Fig 2, no crop")

    shutil.rmtree(tmp, ignore_errors=True)
    print("\nNext: python manuscript/assemble_combined_manuscript_figures.py --only figure2 figure4")


if __name__ == "__main__":
    main()
