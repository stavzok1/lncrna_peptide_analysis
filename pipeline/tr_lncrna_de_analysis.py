"""
Tr-lncRNA identification: per-cancer-type stage and M_stage transitions,
log2FC from group means, gene filtering (expected count >= 1 in >=40% union samples),
z-scores across genes, |z| >= 3.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA

ROOT = REPO_ROOT


import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure_palettes as pal

OUT = ROOT / "tr_lncrna_output"
OUT.mkdir(parents=True, exist_ok=True)

META_STAGE = ["sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage"]
META_META = ["sample_id", "cancer_type", "stage", "M_stage"]

# log2(expected_count + 1) >= 1  <=>  expected_count >= 1
LOG2_THRESH = 1.0
EXPR_FRAC_MIN = 0.4
MIN_SAMPLES_CANCER = 100
Z_ABS_MIN = 3.0  # |z| >= 3 per user (z <= -3 or z >= 3)

STAGE_TRANSITIONS = [
    ("I_II", "I", "II"),
    ("II_III", "II", "III"),
    ("III_IV", "III", "IV"),
    ("E_L", ("I", "II"), ("III", "IV")),
]

# M1_s excluded from all analyses; use exact labels from data
M_TRANSITIONS = [
    ("M0s_M0l", "M0_s", "M0_l"),
    ("M0l_M1L", "M0_l", "M1_L"),
    ("ME_M1L", ("M0_s", "M0_l"), "M1_L"),
]


def cancer_types_over_n(meta: pd.DataFrame, n: int = MIN_SAMPLES_CANCER) -> list[str]:
    counts = meta.groupby("cancer_type").size()
    return counts[counts > n].index.tolist()


def subset_two_groups(
    df: pd.DataFrame,
    gene_cols: list[str],
    col: str,
    earlier,
    later,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (X_early, X_late, sample_ids_union) as float matrices n_samples x n_genes."""
    if isinstance(earlier, str):
        mask_e = df[col].astype(str) == earlier
    else:
        mask_e = df[col].astype(str).isin(list(earlier))
    if isinstance(later, str):
        mask_l = df[col].astype(str) == later
    else:
        mask_l = df[col].astype(str).isin(list(later))
    sub = df.loc[mask_e | mask_l].copy()
    if sub.empty:
        return np.empty((0, len(gene_cols))), np.empty((0, len(gene_cols))), []
    X = sub[gene_cols].to_numpy(dtype=np.float64)
    me = sub[col].astype(str).values
    if isinstance(earlier, str):
        idx_e = me == earlier
    else:
        set_e = set(earlier)
        idx_e = np.array([m in set_e for m in me])
    if isinstance(later, str):
        idx_l = me == later
    else:
        set_l = set(later)
        idx_l = np.array([m in set_l for m in me])
    return X[idx_e], X[idx_l], sub["sample_id"].astype(str).tolist()


def genes_passing_filter(Xe: np.ndarray, Xl: np.ndarray, frac: float = EXPR_FRAC_MIN) -> np.ndarray:
    """Boolean mask length n_genes: >= frac samples in union have expr >= LOG2_THRESH."""
    if Xe.size == 0 and Xl.size == 0:
        return np.zeros(Xe.shape[1] if Xe.size else 0, dtype=bool)
    Xu = np.vstack([Xe, Xl]) if Xe.size and Xl.size else (Xe if Xe.size else Xl)
    n = Xu.shape[0]
    if n == 0:
        return np.zeros(Xu.shape[1], dtype=bool)
    ok = (Xu >= LOG2_THRESH).sum(axis=0) / n >= frac
    return ok


def log2fc_and_z(Xe: np.ndarray, Xl: np.ndarray, mask_genes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """log2FC = mean_late - mean_early for masked genes; z across genes."""
    if not mask_genes.any():
        return np.array([]), np.array([]), np.array([])
    Ge = Xe[:, mask_genes] if Xe.size else np.empty((0, mask_genes.sum()))
    Gl = Xl[:, mask_genes] if Xl.size else np.empty((0, mask_genes.sum()))
    mean_e = Ge.mean(axis=0) if Ge.shape[0] else np.full(mask_genes.sum(), np.nan)
    mean_l = Gl.mean(axis=0) if Gl.shape[0] else np.full(mask_genes.sum(), np.nan)
    logfc = mean_l - mean_e
    valid = np.isfinite(logfc)
    if valid.sum() < 2:
        z = np.full_like(logfc, np.nan)
        return logfc, z, np.where(mask_genes)[0]
    mu = np.nanmean(logfc[valid])
    sd = np.nanstd(logfc[valid], ddof=1)
    if sd == 0 or not np.isfinite(sd):
        z = np.where(valid, 0.0, np.nan)
    else:
        z = (logfc - mu) / sd
    return logfc, z, np.where(mask_genes)[0]


def run_transitions(
    df: pd.DataFrame,
    meta_cols: list[str],
    stage_col: str,
    transitions: list,
    analysis_label: str,
) -> tuple[pd.DataFrame, set[str]]:
    gene_cols = [c for c in df.columns if c not in meta_cols]
    meta_only = df[meta_cols].copy()
    cancers = cancer_types_over_n(meta_only)
    rows = []
    tr_genes: set[str] = set()

    for ctype in sorted(cancers):
        dfc = df[df["cancer_type"] == ctype].copy()
        for tname, early, late in transitions:
            Xe, Xl, _ = subset_two_groups(dfc, gene_cols, stage_col, early, late)
            if Xe.shape[0] < 1 or Xl.shape[0] < 1:
                continue
            mpass = genes_passing_filter(Xe, Xl)
            logfc, z, idx_pass = log2fc_and_z(Xe, Xl, mpass)
            if logfc.size == 0:
                continue
            genes_pass = [gene_cols[i] for i in range(len(gene_cols)) if mpass[i]]
            # idx_pass maps to positions within masked subset
            for j, g in enumerate(genes_pass):
                zz = z[j]
                if np.isfinite(zz) and abs(zz) >= Z_ABS_MIN:
                    tr_genes.add(g)
                    rows.append(
                        {
                            "analysis": analysis_label,
                            "cancer_type": ctype,
                            "transition": tname,
                            "gene": g,
                            "log2FC": float(logfc[j]),
                            "z": float(zz),
                            "n_early": int(Xe.shape[0]),
                            "n_late": int(Xl.shape[0]),
                        }
                    )

    return pd.DataFrame(rows), tr_genes


def plot_counts(df_stage: pd.DataFrame, df_meta: pd.DataFrame, union_genes: set[str]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    def bar_unique_genes(ax, df, title):
        if df.empty:
            ax.text(0.5, 0.5, "No results", ha="center", va="center")
            ax.set_title(title)
            return
        u = df.groupby("transition")["gene"].nunique().sort_index()
        ax.bar(u.index.astype(str), u.values, color=pal.OI_SKY_BLUE)
        ax.set_ylabel("Unique Tr-lncRNAs")
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=45)

    bar_unique_genes(axes[0], df_stage, "Stage transitions (unique genes per transition)")
    bar_unique_genes(axes[1], df_meta, "M_stage transitions (unique genes per transition)")
    axes[2].bar(
        ["Stage (all)", "Metastasis (all)", "Union"],
        [
            df_stage["gene"].nunique() if not df_stage.empty else 0,
            df_meta["gene"].nunique() if not df_meta.empty else 0,
            len(union_genes),
        ],
        color=[pal.OI_BLUISH_GREEN, pal.OI_ORANGE, pal.OI_REDDISH_PURPLE],
    )
    axes[2].set_ylabel("Unique Tr-lncRNAs")
    axes[2].set_title("Overall unique genes")
    plt.tight_layout()
    fig.savefig(OUT / "tr_lncrna_counts.png", dpi=150)
    plt.close()


def plot_per_cancer_tr_counts(
    df_stage: pd.DataFrame,
    df_meta: pd.DataFrame,
    cancers_stage_order: list[str],
    cancers_meta_order: list[str],
) -> None:
    """Bar chart of unique Tr-lncRNAs per cancer (all transitions pooled)."""
    n_s = len(cancers_stage_order)
    n_m = len(cancers_meta_order)
    fig, axes = plt.subplots(1, 2, figsize=(max(12, 0.45 * n_s), 5))

    def panel(ax: plt.Axes, df: pd.DataFrame, cancers: list[str], title_core: str, n_cancers: int) -> None:
        if not cancers:
            ax.text(0.5, 0.5, "No cancers pass threshold", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"{title_core} (N cancers included = {n_cancers})")
            return
        if df.empty:
            counts = pd.Series(0, index=cancers)
        else:
            c = df.groupby("cancer_type")["gene"].nunique()
            counts = c.reindex(cancers, fill_value=0)
        x = np.arange(len(cancers))
        ax.bar(x, counts.values, color=pal.OI_SKY_BLUE)
        ax.set_xticks(x)
        ax.set_xticklabels(cancers, rotation=45, ha="right")
        ax.set_ylabel("Unique Tr-lncRNAs (all transitions)")
        ax.set_title(f"{title_core} (N cancers included = {n_cancers})")

    panel(axes[0], df_stage, cancers_stage_order, "Stage analysis", n_s)
    panel(axes[1], df_meta, cancers_meta_order, "Metastasis (M_stage) analysis", n_m)
    plt.tight_layout()
    fig.savefig(OUT / "tr_lncrna_per_cancer_counts.png", dpi=150)
    plt.close()


def plot_overlap(set_s: set[str], set_m: set[str]) -> None:
    only_s = len(set_s - set_m)
    only_m = len(set_m - set_s)
    both = len(set_s & set_m)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(
        ["Stage only", "Both", "Metastasis only"],
        [only_s, both, only_m],
        color=[pal.OI_BLUISH_GREEN, pal.VOLCANO_NS, pal.OI_ORANGE],
    )
    ax.set_ylabel("Unique Tr-lncRNAs")
    ax.set_title("Overlap between stage-derived and metastasis-derived Tr sets")
    for i, v in enumerate([only_s, both, only_m]):
        ax.text(i, v + max(1, 0.02 * (only_s + both + only_m)), str(v), ha="center", fontsize=10)
    plt.tight_layout()
    fig.savefig(OUT / "tr_lncrna_overlap.png", dpi=150)
    plt.close()


def main() -> None:
    print("Loading stage matrix...")
    df_s = pd.read_csv(DATA / "primary_exp_stage_lncRNA.csv")
    print("Loading metastasis matrix...")
    df_m = pd.read_csv(DATA / "primary_exp_metastasis_lncRNA.csv")

    # Drop M1_s from metastasis df entirely
    df_m = df_m[df_m["M_stage"].astype(str) != "M1_s"].copy()

    print("Running stage transitions...")
    res_s, set_s = run_transitions(df_s, META_STAGE, "stage", STAGE_TRANSITIONS, "stage")
    print("Running M_stage transitions...")
    res_m, set_m = run_transitions(df_m, META_META, "M_stage", M_TRANSITIONS, "metastasis")

    union = set_s | set_m
    cancers_stage = cancer_types_over_n(df_s[META_STAGE])
    cancers_meta = cancer_types_over_n(df_m[META_META])

    res_s.to_csv(OUT / "tr_lncrnas_stage_detail.csv", index=False)
    res_m.to_csv(OUT / "tr_lncrnas_metastasis_detail.csv", index=False)

    summary = {
        "min_samples_per_cancer": MIN_SAMPLES_CANCER,
        "expr_frac_min_union": EXPR_FRAC_MIN,
        "log2_threshold_for_count_ge_1": LOG2_THRESH,
        "z_abs_min": Z_ABS_MIN,
        "unique_tr_genes_stage": sorted(set_s),
        "unique_tr_genes_metastasis": sorted(set_m),
        "unique_tr_genes_union": sorted(union),
        "n_unique_stage": len(set_s),
        "n_unique_metastasis": len(set_m),
        "n_unique_union": len(union),
        "n_rows_stage_detail": int(len(res_s)),
        "n_rows_metastasis_detail": int(len(res_m)),
        "stage_by_transition": res_s.groupby("transition")["gene"].nunique().to_dict() if not res_s.empty else {},
        "meta_by_transition": res_m.groupby("transition")["gene"].nunique().to_dict() if not res_m.empty else {},
        "cancer_types_included_stage": sorted(cancers_stage),
        "cancer_types_included_metastasis": sorted(cancers_meta),
        "n_overlap_stage_and_metastasis": len(set_s & set_m),
        "n_stage_only": len(set_s - set_m),
        "n_metastasis_only": len(set_m - set_s),
    }
    with open(OUT / "tr_lncrna_summary.json", "w", encoding="utf-8") as f:
        json.dump({k: v for k, v in summary.items() if k not in (
            "unique_tr_genes_stage", "unique_tr_genes_metastasis", "unique_tr_genes_union"
        )}, f, indent=2)

    with open(OUT / "tr_genes_stage_unique.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(sorted(set_s)))
    with open(OUT / "tr_genes_metastasis_unique.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(sorted(set_m)))
    with open(OUT / "tr_genes_union.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(sorted(union)))

    plot_counts(res_s, res_m, union)
    plot_overlap(set_s, set_m)
    plot_per_cancer_tr_counts(
        res_s,
        res_m,
        sorted(cancers_stage),
        sorted(cancers_meta),
    )

    print("Done. Output:", OUT)
    print("Unique Tr genes stage:", len(set_s))
    print("Unique Tr genes metastasis:", len(set_m))
    print("Unique Tr genes union:", len(union))


if __name__ == "__main__":
    main()
