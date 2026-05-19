"""
Figure 4A: TIS vs Ribo-seq p-values for analyzed Tr-lncRNA micropeptides (MPs).

**Default (no flags):** every row in ``data/significant_lnc_peptides.tsv`` (~501 exportable
MPs; NetMHC / epitope cohort) with **both** ``TISPvalue`` and ``RiboPvalue`` ≤ 0.05.

Optional ``--cohort tr_lncrna`` re-filters via **GeneSymbol or GeneID** (Ensembl,
versionless; ``smprot_gene_match.py``). ``--cohort tcga_matrix`` uses
``data/smprot_filtered_tcga_expr_genes.tsv``.

Axes are **p-values** on log scales from **10⁻¹²** to **10⁻¹** on both dimensions (nothing
below 10⁻¹² is shown on the axes; ultra–small p-values still participate in ranking).
Shaded
regions: **bluish-green** (TIS ≤ 1e-4 and Ribo ≥ 1e-4), **sky blue** (Ribo ≤ 1e-4 and TIS ≥ 1e-4),
and a **reddish-purple** patch where **both** are ≤ 1e-4 (so LINC00326-type MPs are inside shaded
area). Outer extent and reference lines use **10⁻¹** on both axes.

Framed labels for **PTPRG-AS1**, **LINC00326**, **LINC00958**, and every gene with
≥``--frame-min-mps`` MPs (default 17) at each MP in a shaded region (extreme LINC00958:
highlight dot + arrow). Additional
**plain-text** labels mark the most extreme MPs by combined significance (see
``--top-extreme-labels``). **Violet-corner** lavender labels are **gene symbols only**
(sparse: top ``--top-corner-labels`` in the purple patch, then up to
``--extra-corner-gene-labels`` more among still-unlabeled corner MPs). MPs with **either**
raw p **&lt;** ``P_LO`` (10⁻¹²) are **skipped** for those automatic labels (they sit on
the axis clip) **except** the special **LINC00958** anchor, which may be below the
floor on TIS.

By default the figure is drawn on a **white** background with dark axis labels; use
``--no-light`` for the legacy black theme used in early drafts.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_REPO / "pipeline"), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES
from figure_export import add_publication_args, save_figure_bundle
from smprot_gene_match import peptide_rows_match_significant_lnc

ROOT = REPO_ROOT


import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

import figure_palettes as pal

DATA = ROOT / "data"
FIGURES = ROOT / "figures"
DEFAULT_TSV = DATA / "significant_lnc_peptides.tsv"
TCGA_MATRIX_TSV = DATA / "smprot_filtered_tcga_expr_genes.tsv"
# Analyzed / NetMHC cohort (~501); guard against accidental TCGA-matrix (~2606) defaults.
MANUSCRIPT_COHORT_MAX_MPS = 550

P_BOTH = 0.05
STRONG = 1e-4
P_LO = 1e-12  # lower axis limit (both axes); do not display p < this.
P_HI = 0.1
SCORE_FLOOR = 1e-300  # for log10 in ranking only (not an axis tick)

HIGHLIGHT_GENES = frozenset(pal.F4A_FRAME_STYLES.keys())


def _build_frame_styles(
    genes_arr: np.ndarray,
    *,
    min_mps: int,
    use_light: bool,
) -> tuple[frozenset[str], dict[str, dict[str, str]]]:
    """Genes with framed callouts and per-gene bbox styles (highlights + multi-MP)."""
    base = dict(pal.F4A_LIGHT_FRAME_STYLES if use_light else pal.F4A_FRAME_STYLES)
    extras = pal.F4A_EXTRA_FRAME_STYLES_LIGHT if use_light else pal.F4A_EXTRA_FRAME_STYLES_DARK
    counts: dict[str, int] = {}
    for g in genes_arr:
        s = str(g).strip()
        if s:
            counts[s] = counts.get(s, 0) + 1
    framed = {g for g, n in counts.items() if n >= min_mps} | set(HIGHLIGHT_GENES)
    styles = dict(base)
    extra_genes = sorted(g for g in framed if g not in styles)
    for i, g in enumerate(extra_genes):
        styles[g] = extras[i % len(extras)]
    return frozenset(framed), styles

# Log-space label offsets (decades) around each point; cycle by slot index.
_LABEL_OFFSETS_LOG: tuple[tuple[float, float], ...] = (
    (0.09, 0.07),
    (-0.07, 0.09),
    (0.10, -0.06),
    (-0.08, -0.07),
    (0.06, 0.10),
    (-0.10, 0.06),
    (0.08, 0.08),
    (-0.06, -0.09),
)


def _label_xy_log_near(
    x: float,
    y: float,
    slot: int,
    *,
    max_decade: float = 0.24,
) -> tuple[float, float]:
    """Place label text near (x, y) in log coordinates (never snapped to the 10⁻⁴ border)."""
    lx = float(np.log10(max(x, P_LO)))
    ly = float(np.log10(max(y, SCORE_FLOOR)))
    dx, dy = _LABEL_OFFSETS_LOG[slot % len(_LABEL_OFFSETS_LOG)]
    tx = 10.0 ** (lx + dx)
    ty = 10.0 ** (ly + dy)
    tx = float(np.clip(tx, 10.0 ** (lx - max_decade), 10.0 ** (lx + max_decade)))
    ty = float(np.clip(ty, 10.0 ** (ly - max_decade), 10.0 ** (ly + max_decade)))
    tx = float(np.clip(tx, P_LO * 1.5, P_HI * 0.95))
    ty = float(np.clip(ty, P_LO * 1.5, P_HI * 0.95))
    return tx, ty


def _label_clear_of_left_edge(
    x: float,
    y: float,
    slot: int,
    *,
    margin_decades: float = 1.35,
) -> bool:
    """False when the callout would be clipped by the left (TIS) axis margin."""
    left_limit = 10.0 ** (np.log10(P_LO) + margin_decades)
    if x < left_limit:
        return False
    tx, _ = _label_xy_log_near(x, y, slot)
    if tx < left_limit:
        return False
    if tx < x and tx * (10.0**0.55) < left_limit:
        return False
    return True


def _annotate_gene_callout(
    ax,
    gene: str,
    x: float,
    y: float,
    slot: int,
    *,
    lbl_c: str,
    label_accent: str,
    fontsize: float,
    frame_style: dict[str, str] | None = None,
    with_arrow: bool = False,
    zorder: int = 4,
) -> None:
    tx, ty = _label_xy_log_near(x, y, slot)
    ha = "left" if tx >= x else "right"
    va = "bottom" if ty >= y else "top"
    kw: dict = dict(
        xy=(x, y),
        xytext=(tx, ty),
        textcoords="data",
        fontsize=fontsize,
        color=lbl_c if frame_style else label_accent,
        ha=ha,
        va=va,
        zorder=zorder,
        clip_on=True,
    )
    if frame_style:
        kw["bbox"] = dict(
            boxstyle="round,pad=0.22",
            fc=frame_style["facecolor"],
            ec=frame_style["edgecolor"],
            lw=1.2,
        )
    if with_arrow:
        kw["arrowprops"] = dict(
            arrowstyle="-",
            color=label_accent,
            lw=0.65,
            alpha=0.75,
            shrinkA=2,
            shrinkB=3,
            connectionstyle="arc3,rad=0.08",
        )
    ax.annotate(str(gene), **kw)


def load_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    need = {"GeneSymbol", "TISPvalue", "RiboPvalue"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{path}: missing columns {sorted(missing)}")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Fig 4A TIS vs Ribo-seq p-values (501 significant_lnc_peptides.tsv by default)."
    )
    ap.add_argument(
        "--peptides-tsv",
        type=Path,
        default=DEFAULT_TSV,
        help=f"Peptide table (default: {DEFAULT_TSV.name}, ~501 MPs).",
    )
    ap.add_argument(
        "--cohort",
        choices=("tr_lncrna", "tcga_matrix"),
        default=None,
        help="Optional filter. tr_lncrna: symbol or Ensembl vs significant_lncs / canonical "
        "Tr genes. tcga_matrix: read smprot_filtered_tcga_expr_genes.tsv (all rows). "
        "Default: all rows in --peptides-tsv.",
    )
    ap.add_argument("--out-dir", type=Path, default=FIGURES)
    ap.add_argument("--out-name", type=str, default="fig4a_tr_lncrna_mp_tis_vs_riboseq_pvalues.png")
    ap.add_argument(
        "--top-extreme-labels",
        type=int,
        default=28,
        metavar="N",
        help="Label the N most extreme MPs outside the violet corner (by -log10 TIS - log10 Ribo).",
    )
    ap.add_argument(
        "--top-corner-labels",
        type=int,
        default=8,
        metavar="N",
        help="Label the N strongest MPs in the violet corner (after filters).",
    )
    ap.add_argument(
        "--extra-corner-gene-labels",
        type=int,
        default=8,
        metavar="M",
        help="Additional sparse lavender gene labels in the violet corner (0 = none). MPs with "
        "either raw p < axis floor 10⁻¹² are skipped except the red LINC00958 anchor.",
    )
    ap.add_argument(
        "--no-light",
        action="store_true",
        help="Use the legacy black background instead of the default white figure.",
    )
    ap.add_argument(
        "--frame-min-mps",
        type=int,
        default=17,
        metavar="N",
        help="Framed bbox labels for every gene with at least N MPs in the plot (default 17). "
        "PTPRG-AS1, LINC00326, LINC00958 are always framed when in a shaded region.",
    )
    add_publication_args(ap)
    args = ap.parse_args()
    use_light = not args.no_light

    if args.cohort == "tcga_matrix":
        args.peptides_tsv = TCGA_MATRIX_TSV
    elif args.peptides_tsv.resolve() == TCGA_MATRIX_TSV.resolve():
        raise SystemExit(
            f"Fig 4A manuscript panel uses {DEFAULT_TSV.name} (~501 significant Tr MPs). "
            f"Pass --cohort tcga_matrix only for the broader TCGA-matrix export (~2606 MPs)."
        )

    if not args.peptides_tsv.exists():
        raise SystemExit(f"Missing {args.peptides_tsv}")

    df = load_table(args.peptides_tsv)
    if args.cohort == "tr_lncrna":
        m = peptide_rows_match_significant_lnc(df)
        df = df.loc[m].copy()
        if df.empty:
            raise SystemExit(
                f"No rows in {args.peptides_tsv} match Tr lncRNAs via GeneSymbol or GeneID "
                f"(see pipeline/smprot_gene_match.py)."
            )
    df = df[(df["TISPvalue"] <= P_BOTH) & (df["RiboPvalue"] <= P_BOTH)].copy()
    df = df.dropna(subset=["TISPvalue", "RiboPvalue"]).copy()
    if df.empty:
        raise SystemExit("No rows after p ≤ 0.05 on both axes.")

    is_manuscript_default = (
        args.cohort is None and args.peptides_tsv.resolve() == DEFAULT_TSV.resolve()
    )
    if is_manuscript_default and len(df) > MANUSCRIPT_COHORT_MAX_MPS:
        raise SystemExit(
            f"Expected ≤{MANUSCRIPT_COHORT_MAX_MPS} MPs from {DEFAULT_TSV.name} "
            f"(got {len(df)}). Use --cohort tcga_matrix for the full TCGA-matrix cohort."
        )

    floor = 1e-300
    tis_raw = df["TISPvalue"].to_numpy(dtype=float)
    ribo_raw = df["RiboPvalue"].to_numpy(dtype=float)
    tis = np.clip(tis_raw, floor, 1.0)
    ribo = np.clip(ribo_raw, floor, 1.0)
    genes = df["GeneSymbol"].to_numpy()

    # Most TIS-extreme LINC00958 MP → red highlight (may sit on the left margin if TIS < 10⁻¹²).
    lin_mask = genes == "LINC00958"
    extreme_idx: int | None = None
    idx_lin = np.where(lin_mask)[0]
    if len(idx_lin):
        extreme_idx = int(idx_lin[int(np.argmin(tis[idx_lin]))])

    bulk_mask = np.ones(len(df), dtype=bool)
    if extreme_idx is not None:
        bulk_mask[extreme_idx] = False

    # Shaded "strong single evidence" bands plus both-strong corner (for labeling LINC00326-type MPs).
    in_green = (tis <= STRONG) & (ribo >= STRONG)
    in_blue = (tis >= STRONG) & (ribo <= STRONG)
    in_corner = (tis < STRONG) & (ribo < STRONG)
    in_any_shade = in_green | in_blue | in_corner

    def allow_auto_label(i: int) -> bool:
        """False for MPs pinned below the 10⁻¹² axis clip (unless red LINC00958 anchor)."""
        if extreme_idx is not None and i == extreme_idx:
            return True
        return bool(tis_raw[i] >= P_LO and ribo_raw[i] >= P_LO)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    fig_bg = pal.F4A_LIGHT_AX if use_light else "black"
    ax_bg = pal.F4A_LIGHT_AX if use_light else "black"
    fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=150, facecolor=fig_bg)
    ax.set_facecolor(ax_bg)

    framed_genes, frame_styles = _build_frame_styles(
        genes, min_mps=max(1, args.frame_min_mps), use_light=use_light
    )

    # Polygons from p = 10⁻¹² to 10⁻¹; outer boundary at P_HI = 0.1.
    shade_alpha = pal.F4A_LIGHT_SHADE_ALPHA if use_light else 0.24
    green = Polygon(
        [(P_LO, STRONG), (STRONG, STRONG), (STRONG, P_HI), (P_LO, P_HI)],
        closed=True,
        facecolor=pal.F4A_SHADE_TIS_STRONG,
        edgecolor="#7fbf9f",
        linewidth=1.0,
        alpha=shade_alpha,
        zorder=1,
    )
    blue = Polygon(
        [(STRONG, P_LO), (P_HI, P_LO), (P_HI, STRONG), (STRONG, STRONG)],
        closed=True,
        facecolor=pal.F4A_SHADE_RIBO_STRONG,
        edgecolor="#9ecae1",
        linewidth=1.0,
        alpha=shade_alpha,
        zorder=1,
    )
    corner = Polygon(
        [(P_LO, P_LO), (STRONG, P_LO), (STRONG, STRONG), (P_LO, STRONG)],
        closed=True,
        facecolor=pal.F4A_SHADE_BOTH_STRONG,
        edgecolor="#d4a5df",
        linewidth=0.9,
        alpha=shade_alpha * 0.75,
        zorder=1,
    )
    ax.add_patch(green)
    ax.add_patch(blue)
    ax.add_patch(corner)

    # Reference lines at p = 1e-4 from 10⁻¹² to 10⁻¹.
    grid_c = pal.F4A_LIGHT_GRID_MAJOR if use_light else pal.F4A_GRID
    grid_min_c = pal.F4A_LIGHT_GRID_MINOR if use_light else pal.F4A_GRID
    frame_c = pal.F4A_LIGHT_FRAME if use_light else pal.F4A_FRAME
    ax.plot([P_LO, P_HI], [STRONG, STRONG], color=grid_c, lw=1.05, alpha=0.55, zorder=2)
    ax.plot([STRONG, STRONG], [P_LO, P_HI], color=grid_c, lw=1.05, alpha=0.55, zorder=2)

    # Outer frame at p = 10⁻¹ (top and right) so the 10⁻¹ boundary reads clearly.
    ax.plot(
        [P_LO, P_HI],
        [P_HI, P_HI],
        color=frame_c,
        lw=1.35,
        alpha=0.9 if use_light else 0.75,
        solid_capstyle="round",
        zorder=2,
    )
    ax.plot(
        [P_HI, P_HI],
        [P_LO, P_HI],
        color=frame_c,
        lw=1.35,
        alpha=0.9 if use_light else 0.75,
        solid_capstyle="round",
        zorder=2,
    )

    bulk_c = pal.F4A_LIGHT_POINTS_BULK if use_light else pal.F4A_POINTS_BULK
    bulk_alpha = 0.45 if use_light else 0.35
    ax.scatter(
        tis[bulk_mask],
        ribo[bulk_mask],
        s=22,
        c=bulk_c,
        alpha=bulk_alpha,
        linewidths=0,
        zorder=3,
        rasterized=True,
    )
    if extreme_idx is not None:
        te = float(tis[extreme_idx])
        re = float(ribo[extreme_idx])
        # If TIS is below the axis floor, pin a marker on the left spine (true Ribo unchanged).
        xe = P_LO if te < P_LO else te
        ax.scatter(
            [xe],
            [re],
            s=110,
            c=pal.F4A_POINT_EXTREME,
            alpha=0.98,
            linewidths=0.75,
            edgecolors=frame_c,
            marker="o",
            zorder=6,
            clip_on=False,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(P_LO, P_HI)
    ax.set_ylim(P_LO, P_HI)
    # Explicit decade ticks 10^-12 … 10^-1 so $10^{-1}$ appears on both axes.
    decade_ticks = [10**k for k in range(-12, 0)]
    ax.set_xticks(decade_ticks)
    ax.set_yticks(decade_ticks)
    ax.xaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10))
    ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10))
    lbl_c = pal.F4A_LIGHT_TEXT if use_light else "0.9"
    tick_c = pal.F4A_LIGHT_TICK if use_light else "0.85"
    label_accent = pal.OI_VERMILLION if use_light else pal.OI_YELLOW
    ax.set_xlabel("P-value of TIS (translation-initiation site)", color=lbl_c, fontsize=11)
    ax.set_ylabel("P-value of Ribo-seq", color=lbl_c, fontsize=11)
    ax.tick_params(colors=tick_c, which="both", labelsize=8)
    spine_c = pal.F4A_LIGHT_TICK if use_light else "0.5"
    for spine in ax.spines.values():
        spine.set_color(spine_c)
    ax.grid(True, which="major", linestyle="-", alpha=0.35 if use_light else 0.25, color=grid_c)
    ax.grid(True, which="minor", linestyle=":", alpha=0.22 if use_light else 0.12, color=grid_min_c)

    # Track row indices that already carry a callout (same MP not labeled twice).
    labeled: set[int] = set()
    label_slot = 0

    # Framed labels: highlight trio + all genes with ≥ frame-min-mps (each MP in a shaded region).
    for i, g in enumerate(genes):
        gs = str(g)
        if gs not in framed_genes:
            continue
        if not bool(in_any_shade[i]):
            continue
        if gs == "LINC00958" and i == extreme_idx:
            continue  # extreme MP: highlight dot + arrow label
        if i in labeled:
            continue
        sty = frame_styles[gs]
        slot_i = label_slot
        if not _label_clear_of_left_edge(float(tis[i]), float(ribo[i]), slot_i):
            continue
        _annotate_gene_callout(
            ax,
            gs,
            float(tis[i]),
            float(ribo[i]),
            slot_i,
            lbl_c=lbl_c,
            label_accent=label_accent,
            fontsize=7 if gs in HIGHLIGHT_GENES else 6.5,
            frame_style=sty,
            with_arrow=bool(in_corner[i]),
            zorder=5,
        )
        labeled.add(i)
        label_slot += 1

    if extreme_idx is not None and bool(in_any_shade[extreme_idx]):
        te = float(tis[extreme_idx])
        rib_e = float(ribo[extreme_idx])
        xe = P_LO if te < P_LO else te
        # Tight callout: small multiplicative offsets in log space (label sits just NE of the point).
        lx = np.log10(max(float(xe), P_LO))
        ly = np.log10(max(rib_e, SCORE_FLOOR))
        tx = min(10 ** (lx + 0.18), P_HI * 0.92)
        ty = min(10 ** (ly + 0.10), P_HI * 0.92)
        ax.annotate(
            "LINC00958",
            xy=(xe, rib_e),
            xytext=(tx, ty),
            textcoords="data",
            fontsize=8.5,
            color=lbl_c,
            ha="left",
            va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.22",
                fc=frame_styles["LINC00958"]["facecolor"],
                ec=frame_styles["LINC00958"]["edgecolor"],
                lw=1.5,
            ),
            arrowprops=dict(
                arrowstyle="-|>",
                color=pal.OI_SKY_BLUE,
                lw=1.15,
                shrinkA=1,
                shrinkB=2,
                mutation_scale=12,
                connectionstyle="arc3,rad=0.05",
            ),
            zorder=7,
            clip_on=True,
        )
        labeled.add(extreme_idx)

    log10_t = np.log10(np.maximum(tis, SCORE_FLOOR))
    log10_r = np.log10(np.maximum(ribo, SCORE_FLOOR))
    score = -(log10_t + log10_r)

    def _label_corner_candidate(j: int, g: str, slot: int) -> bool:
        if j in labeled:
            return False
        if extreme_idx is not None and j == extreme_idx:
            return False
        if g in framed_genes:
            return False
        if not allow_auto_label(j):
            return False
        if not _label_clear_of_left_edge(float(tis[j]), float(ribo[j]), slot):
            return False
        return True

    # Violet corner: strongest MPs first (multiple MPs per gene allowed); log-space anchor + arrow.
    n_corner = max(0, args.top_corner_labels) + max(0, args.extra_corner_gene_labels)
    if n_corner > 0:
        corner_idx = np.where(in_corner)[0]
        if len(corner_idx):
            sub_order = corner_idx[np.argsort(-score[corner_idx])]
            ck = 0
            for j in sub_order:
                if ck >= n_corner:
                    break
                g = str(genes[j])
                if not _label_corner_candidate(j, g, label_slot + ck):
                    continue
                _annotate_gene_callout(
                    ax,
                    g,
                    float(tis[j]),
                    float(ribo[j]),
                    label_slot + ck,
                    lbl_c=lbl_c,
                    label_accent=label_accent,
                    fontsize=6.0 if ck < args.top_corner_labels else 5.8,
                    with_arrow=True,
                    zorder=4,
                )
                labeled.add(j)
                ck += 1

    # Plain-text labels for extreme MPs outside the violet corner (multiple MPs per gene allowed).
    if args.top_extreme_labels > 0:
        order = np.argsort(-score)
        k = 0
        for j in order:
            if k >= args.top_extreme_labels:
                break
            if bool(in_corner[j]):
                continue
            if j in labeled:
                continue
            g = str(genes[j])
            if g in framed_genes and bool(in_any_shade[j]):
                continue
            if extreme_idx is not None and j == extreme_idx:
                continue
            if not allow_auto_label(j):
                continue
            slot_k = label_slot + k
            if not _label_clear_of_left_edge(float(tis[j]), float(ribo[j]), slot_k):
                continue
            _annotate_gene_callout(
                ax,
                g,
                float(tis[j]),
                float(ribo[j]),
                slot_k,
                lbl_c=lbl_c,
                label_accent=label_accent,
                fontsize=5.8,
                with_arrow=False,
                zorder=4,
            )
            labeled.add(j)
            k += 1

    fig.tight_layout()
    out = args.out_dir / args.out_name
    save_figure_bundle(
        fig,
        out,
        png_dpi=150,
        publication_dir=args.publication_dir,
        publication_tiff_kind=args.publication_tiff_kind,
        figures_root=FIGURES,
        bbox_inches="tight",
        facecolor=fig_bg,
    )
    plt.close(fig)
    if args.cohort is None and is_manuscript_default:
        cohort_note = (
            f"{args.peptides_tsv.name} (significant Tr MPs; TIS<={P_BOTH:g} & Ribo<={P_BOTH:g})"
        )
    elif args.cohort is None:
        cohort_note = f"{args.peptides_tsv.name} (TIS<={P_BOTH:g} & Ribo<={P_BOTH:g})"
    elif args.cohort == "tr_lncrna":
        cohort_note = "Tr lncRNA (symbol+Ensembl)"
    else:
        cohort_note = "all TCGA-matrix"
    print(
        f"Wrote {out} ({len(df)} MPs, filter={cohort_note}, "
        f"framed genes={len(framed_genes)})"
    )


if __name__ == "__main__":
    main()
