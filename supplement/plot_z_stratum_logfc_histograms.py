"""
Histograms of log2FC (mean late − mean early on log2(expected count+1) scale) for every
lncRNA that **enters the z-score calculation** in each stratum: same filters as
``tr_lncrna_de_analysis.py`` (≥40% of union samples with log2(expr+1)≥1), per
cancer type × transition × analysis (stage / metastasis).

Outputs under ``figures/z_stratum_logfc_histograms/``:
  - One PNG per (analysis, transition) with a grid of per-cancer histograms.
  - ``stratum_logfc_summary.csv`` — counts and basic moments per stratum.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_REPO / "pipeline")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES

ROOT = REPO_ROOT


import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure_palettes as pal
from tr_lncrna_de_analysis import (
    DATA,
    META_META,
    META_STAGE,
    M_TRANSITIONS,
    STAGE_TRANSITIONS,
    cancer_types_over_n,
    genes_passing_filter,
    log2fc_and_z,
    subset_two_groups,
)

OUT_DIR = ROOT / "figures" / "z_stratum_logfc_histograms"


def collect_stratum_logfc(
    df: pd.DataFrame,
    meta_cols: list[str],
    stage_col: str,
    cancer_type: str,
    transition: tuple,
) -> np.ndarray:
    tname, early, late = transition
    gene_cols = [c for c in df.columns if c not in meta_cols]
    dfc = df.loc[df["cancer_type"].astype(str) == str(cancer_type)].copy()
    Xe, Xl, _ = subset_two_groups(dfc, gene_cols, stage_col, early, late)
    if Xe.shape[0] < 1 or Xl.shape[0] < 1:
        return np.array([])
    mpass = genes_passing_filter(Xe, Xl)
    logfc, _, _ = log2fc_and_z(Xe, Xl, mpass)
    if logfc.size == 0:
        return np.array([])
    return logfc[np.isfinite(logfc)]


def plot_transition_grid(
    df: pd.DataFrame,
    meta_cols: list[str],
    stage_col: str,
    analysis: str,
    transition: tuple,
    cancers: list[str],
    out_path: Path,
) -> list[dict]:
    tname, _, _ = transition
    n_c = len(cancers)
    ncols = 4
    nrows = max(1, math.ceil(n_c / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.4 * nrows), sharey=False)
    axes_flat = np.atleast_1d(axes).ravel()
    summary: list[dict] = []

    for ax in axes_flat[n_c:]:
        ax.set_visible(False)

    for idx, ct in enumerate(cancers):
        ax = axes_flat[idx]
        vals = collect_stratum_logfc(df, meta_cols, stage_col, ct, transition)
        summary.append(
            {
                "analysis": analysis,
                "transition": tname,
                "cancer_type": ct,
                "n_genes_z_stratum": int(vals.size),
                "log2fc_mean": float(vals.mean()) if vals.size else float("nan"),
                "log2fc_std": float(vals.std(ddof=1)) if vals.size > 1 else float("nan"),
            }
        )
        if vals.size == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes, fontsize=9)
        else:
            ax.hist(vals, bins=min(40, max(10, int(np.sqrt(vals.size)))), color=pal.OI_SKY_BLUE, edgecolor="white", linewidth=0.4)
        ax.set_title(str(ct), fontsize=9)
        ax.set_xlabel("log2FC", fontsize=8)
        ax.tick_params(axis="both", labelsize=7)

    fig.suptitle(
        f"{analysis}: log2FC among genes used for z-scores — {tname}\n"
        f"(expr filter: ≥40% union samples with log2(count+1)≥1; transition arms need ≥1 sample each)",
        fontsize=11,
        y=1.02,
    )
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df_s = pd.read_csv(DATA / "primary_exp_stage_lncRNA.csv")
    df_m = pd.read_csv(DATA / "primary_exp_metastasis_lncRNA.csv")
    df_m = df_m.loc[df_m["M_stage"].astype(str) != "M1_s"].copy()

    cancers_s = sorted(cancer_types_over_n(df_s[META_STAGE]))
    cancers_m = sorted(cancer_types_over_n(df_m[META_META]))

    all_summary: list[dict] = []

    for tr in STAGE_TRANSITIONS:
        rows = plot_transition_grid(
            df_s, META_STAGE, "stage", "stage", tr, cancers_s, OUT_DIR / f"stage_{tr[0]}_log2fc_by_cancer.png"
        )
        all_summary.extend(rows)

    for tr in M_TRANSITIONS:
        rows = plot_transition_grid(
            df_m, META_META, "M_stage", "metastasis", tr, cancers_m, OUT_DIR / f"metastasis_{tr[0]}_log2fc_by_cancer.png"
        )
        all_summary.extend(rows)

    sum_df = pd.DataFrame(all_summary)
    sum_df.to_csv(OUT_DIR / "stratum_logfc_summary.csv", index=False)
    print(f"Wrote {len(all_summary)} stratum summaries -> {OUT_DIR / 'stratum_logfc_summary.csv'}")
    print(f"PNGs under {OUT_DIR}")


if __name__ == "__main__":
    main()
