"""
Figure 4A: TIS vs Ribo-seq p-values for TCGA-matrix filtered lncRNA micropeptides (MPs).

Each point is one row from ``data/smprot_filtered_tcga_expr_genes.tsv`` (SmProt-filtered
lncRNA peptides whose ``GeneSymbol`` is a column in the TCGA expression matrices), with
**both** ``TISPvalue`` and ``RiboPvalue`` ≤ 0.05. This is **not** restricted to the
canonical Tr (limma∩z) gene list.

Axes are **p-values** on log scales from **10⁻¹²** to **10⁻¹** on both dimensions (nothing
below 10⁻¹² is shown on the axes; ultra–small p-values still participate in ranking).
Shaded
regions: **bluish-green** (TIS ≤ 1e-4 and Ribo ≥ 1e-4), **sky blue** (Ribo ≤ 1e-4 and TIS ≥ 1e-4),
and a **reddish-purple** patch where **both** are ≤ 1e-4 (so LINC00326-type MPs are inside shaded
area). Outer extent and reference lines use **10⁻¹** on both axes.

Framed labels for **PTPRG-AS1**, **LINC00326**, and **LINC00958** only at MPs lying in
one of those shaded regions (the extreme LINC00958 uses the highlight dot + arrow). Additional
**plain-text** labels mark the most extreme MPs by combined significance (see
``--top-extreme-labels``). **Violet-corner** lavender labels are **gene symbols only**
(sparse: top ``--top-corner-labels`` in the purple patch, then up to
``--extra-corner-gene-labels`` more among still-unlabeled corner MPs). MPs with **either**
raw p **&lt;** ``P_LO`` (10⁻¹²) are **skipped** for those automatic labels (they sit on
the axis clip) **except** the special **LINC00958** anchor, which may be below the
floor on TIS.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES

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
DEFAULT_TSV = DATA / "smprot_filtered_tcga_expr_genes.tsv"

P_BOTH = 0.05
STRONG = 1e-4
P_LO = 1e-12  # lower axis limit (both axes); do not display p < this.
P_HI = 0.1
SCORE_FLOOR = 1e-300  # for log10 in ranking only (not an axis tick)

THREE_GENES = frozenset(pal.F4A_FRAME_STYLES.keys())


def load_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    need = {"GeneSymbol", "TISPvalue", "RiboPvalue"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{path}: missing columns {sorted(missing)}")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Fig 4A TIS vs Ribo-seq p-values (TCGA-matrix filtered lncRNA MPs)."
    )
    ap.add_argument("--peptides-tsv", type=Path, default=DEFAULT_TSV)
    ap.add_argument("--out-dir", type=Path, default=FIGURES)
    ap.add_argument("--out-name", type=str, default="fig4a_tr_lncrna_mp_tis_vs_riboseq_pvalues.png")
    ap.add_argument(
        "--top-extreme-labels",
        type=int,
        default=28,
        metavar="N",
        help="Label the N most extreme MPs (by -log10 TIS - log10 Ribo), plain text.",
    )
    ap.add_argument(
        "--top-corner-labels",
        type=int,
        default=8,
        metavar="N",
        help="Lavender gene-only: label the N strongest MPs in the violet corner (after filters).",
    )
    ap.add_argument(
        "--extra-corner-gene-labels",
        type=int,
        default=8,
        metavar="M",
        help="Additional sparse lavender gene labels in the violet corner (0 = none). MPs with "
        "either raw p < axis floor 10⁻¹² are skipped except the red LINC00958 anchor.",
    )
    args = ap.parse_args()

    if not args.peptides_tsv.exists():
        raise SystemExit(f"Missing {args.peptides_tsv}")

    df = load_table(args.peptides_tsv)
    df = df[(df["TISPvalue"] <= P_BOTH) & (df["RiboPvalue"] <= P_BOTH)].copy()
    df = df.dropna(subset=["TISPvalue", "RiboPvalue"]).copy()
    if df.empty:
        raise SystemExit("No rows after p ≤ 0.05 on both axes.")

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

    fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=150, facecolor="black")
    ax.set_facecolor("black")

    # Polygons from p = 10⁻¹² to 10⁻¹; outer boundary at P_HI = 0.1.
    green = Polygon(
        [(P_LO, STRONG), (STRONG, STRONG), (STRONG, P_HI), (P_LO, P_HI)],
        closed=True,
        facecolor=pal.F4A_SHADE_TIS_STRONG,
        edgecolor="#7fbf9f",
        linewidth=1.0,
        alpha=0.24,
        zorder=1,
    )
    blue = Polygon(
        [(STRONG, P_LO), (P_HI, P_LO), (P_HI, STRONG), (STRONG, STRONG)],
        closed=True,
        facecolor=pal.F4A_SHADE_RIBO_STRONG,
        edgecolor="#9ecae1",
        linewidth=1.0,
        alpha=0.24,
        zorder=1,
    )
    corner = Polygon(
        [(P_LO, P_LO), (STRONG, P_LO), (STRONG, STRONG), (P_LO, STRONG)],
        closed=True,
        facecolor=pal.F4A_SHADE_BOTH_STRONG,
        edgecolor="#d4a5df",
        linewidth=0.9,
        alpha=0.18,
        zorder=1,
    )
    ax.add_patch(green)
    ax.add_patch(blue)
    ax.add_patch(corner)

    # Reference lines at p = 1e-4 from 10⁻¹² to 10⁻¹.
    ax.plot([P_LO, P_HI], [STRONG, STRONG], color=pal.F4A_GRID, lw=1.05, alpha=0.55, zorder=2)
    ax.plot([STRONG, STRONG], [P_LO, P_HI], color=pal.F4A_GRID, lw=1.05, alpha=0.55, zorder=2)

    # Outer frame at p = 10⁻¹ (top and right) so the 10⁻¹ boundary reads clearly.
    ax.plot(
        [P_LO, P_HI],
        [P_HI, P_HI],
        color=pal.F4A_FRAME,
        lw=1.35,
        alpha=0.75,
        solid_capstyle="round",
        zorder=2,
    )
    ax.plot(
        [P_HI, P_HI],
        [P_LO, P_HI],
        color=pal.F4A_FRAME,
        lw=1.35,
        alpha=0.75,
        solid_capstyle="round",
        zorder=2,
    )

    ax.scatter(
        tis[bulk_mask],
        ribo[bulk_mask],
        s=22,
        c=pal.F4A_POINTS_BULK,
        alpha=0.35,
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
            edgecolors=pal.F4A_FRAME,
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
    ax.set_xlabel("P-value of TIS (translation-initiation site)", color="0.9", fontsize=11)
    ax.set_ylabel("P-value of Ribo-seq", color="0.9", fontsize=11)
    ax.tick_params(colors="0.85", which="both", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("0.5")
    ax.grid(True, which="major", linestyle="-", alpha=0.25, color="0.65")
    ax.grid(True, which="minor", linestyle=":", alpha=0.12, color="0.5")

    # Track row indices that already carry a text/marker callout (avoid duplicate labels).
    labeled: set[int] = set()

    # Framed labels for the three highlighted genes only inside shaded regions (not bulk quadrant).
    for i, g in enumerate(genes):
        if g not in THREE_GENES:
            continue
        if not bool(in_any_shade[i]):
            continue
        if g == "LINC00958" and i == extreme_idx:
            continue  # extreme MP: highlight dot + arrow label
        sty = pal.F4A_FRAME_STYLES[g]
        ax.annotate(
            str(g),
            xy=(tis[i], ribo[i]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
            color="0.95",
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.22", fc=sty["facecolor"], ec=sty["edgecolor"], lw=1.2),
            zorder=5,
            # LINC00326 sits on the right: allow full bbox past the axis edge; other framed genes stay clipped.
            clip_on=g != "LINC00326",
        )
        labeled.add(i)

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
            color="0.95",
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.22", fc="#050a28", ec=pal.OI_SKY_BLUE, lw=1.5),
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

    # Plain-text labels for the most extreme MPs (gene symbol only; no peptide IDs).
    if args.top_extreme_labels > 0:
        log10_t = np.log10(np.maximum(tis, SCORE_FLOOR))
        log10_r = np.log10(np.maximum(ribo, SCORE_FLOOR))
        score = -(log10_t + log10_r)
        order = np.argsort(-score)
        k = 0
        for j in order:
            if k >= args.top_extreme_labels:
                break
            g = str(genes[j])
            if g in THREE_GENES and bool(in_any_shade[j]):
                continue  # already framed (or red marker + callout for extreme LINC00958)
            if extreme_idx is not None and j == extreme_idx:
                continue
            if not allow_auto_label(j):
                continue
            ax.annotate(
                g,
                xy=(tis[j], ribo[j]),
                xytext=(5 + (k % 3) * 3, 4 + (k % 6) * 2),
                textcoords="offset points",
                fontsize=5.8,
                color=pal.OI_YELLOW,
                ha="left",
                va="bottom",
                alpha=0.92,
                zorder=4,
                clip_on=True,
            )
            labeled.add(j)
            k += 1

    # Lavender gene-only labels in the violet corner (sparse).
    if args.top_corner_labels > 0:
        log10_t = np.log10(np.maximum(tis, SCORE_FLOOR))
        log10_r = np.log10(np.maximum(ribo, SCORE_FLOOR))
        score = -(log10_t + log10_r)
        corner_idx = np.where(in_corner)[0]
        if len(corner_idx):
            sub_order = corner_idx[np.argsort(-score[corner_idx])]
            ck = 0
            for j in sub_order:
                if ck >= args.top_corner_labels:
                    break
                g = str(genes[j])
                if extreme_idx is not None and j == extreme_idx:
                    continue
                if g in THREE_GENES:
                    continue  # framed / callout; no duplicate lavender gene
                if not allow_auto_label(j):
                    continue
                ax.annotate(
                    g,
                    xy=(tis[j], ribo[j]),
                    xytext=(-10 - (ck % 3) * 16, -8 - (ck % 4) * 14),
                    textcoords="offset points",
                    fontsize=6.0,
                    color=pal.OI_YELLOW,
                    ha="right",
                    va="top",
                    alpha=0.95,
                    zorder=4,
                    clip_on=True,
                )
                labeled.add(j)
                ck += 1

    # A few more corner MPs (other lncRNAs), still sparse and still ≥10⁻¹² on both raw p.
    if args.extra_corner_gene_labels > 0:
        log10_t_all = np.log10(np.maximum(tis, SCORE_FLOOR))
        log10_r_all = np.log10(np.maximum(ribo, SCORE_FLOOR))
        score_all = -(log10_t_all + log10_r_all)
        rest: list[int] = []
        for i in np.where(in_corner)[0]:
            g = str(genes[i])
            if g in THREE_GENES:
                continue
            if extreme_idx is not None and i == extreme_idx:
                continue
            if i in labeled:
                continue
            if not allow_auto_label(i):
                continue
            rest.append(i)
        rest.sort(key=lambda i: float(score_all[i]), reverse=True)
        for uk, i in enumerate(rest[: args.extra_corner_gene_labels]):
            g = str(genes[i])
            ax.annotate(
                g,
                xy=(tis[i], ribo[i]),
                xytext=(-12 - (uk % 3) * 18, -10 - (uk % 4) * 16),
                textcoords="offset points",
                fontsize=5.8,
                color=pal.OI_YELLOW,
                ha="right",
                va="top",
                alpha=0.9,
                zorder=4,
                clip_on=True,
            )
            labeled.add(i)


    fig.tight_layout()
    out = args.out_dir / args.out_name
    fig.savefig(out, bbox_inches="tight", facecolor="black")
    plt.close(fig)
    print(f"Wrote {out} ({len(df)} MPs)")


if __name__ == "__main__":
    main()
