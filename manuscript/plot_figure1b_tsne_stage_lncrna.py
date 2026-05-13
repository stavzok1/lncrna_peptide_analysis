"""
Figure 1B — t-SNE on TCGA primary tumor samples with AJCC stage (lncRNA expression matrix).

Loads ``data/primary_exp_stage_lncRNA.csv`` (same samples as ``tr_lncrna_de_analysis``),
standardizes genes across samples, optionally PCA-prelimits dimensionality, then fits a
**4-dimensional** t-SNE (OpenTSNE; sklearn Barnes-Hut does not support ``n_components > 3``).

Writes four PNGs under ``figures/``:

  - components **1 vs 2**, colored by **cancer_type**
  - components **1 vs 2**, colored by **stage**
  - components **3 vs 4**, colored by **cancer_type**
  - components **3 vs 4**, colored by **stage**
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import DATA, FIGURES, REPO_ROOT

ROOT = REPO_ROOT


import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    from openTSNE import TSNE as OpenTSNE
except ImportError as e:  # pragma: no cover
    raise SystemExit(
        "openTSNE is required for 4D t-SNE at this sample size. Install: pip install opentsne"
    ) from e

DEFAULT_CSV = DATA / "primary_exp_stage_lncRNA.csv"

META_COLS = ("sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage")


def _categorical_colors(labels: np.ndarray) -> tuple[dict[str, tuple], list[str]]:
    """Map each unique string label to an RGBA color; return (color_by_label, uniques_sorted)."""
    uniques = sorted({str(x) for x in labels})
    # Enough distinct colors for ~28 cancer types (tab20 + tab20b)
    c20 = plt.colormaps["tab20"]
    c20b = plt.colormaps["tab20b"]
    colors: list[tuple] = []
    for i in range(len(uniques)):
        if i < 20:
            colors.append(tuple(float(x) for x in c20(i / 19.0)))
        else:
            colors.append(tuple(float(x) for x in c20b((i - 20) / 19.0)))
    return {u: colors[j] for j, u in enumerate(uniques)}, uniques


def _stage_colors() -> dict[str, tuple]:
    return {
        "I": tuple(float(x) for x in plt.colormaps["viridis"](0.15)),
        "II": tuple(float(x) for x in plt.colormaps["viridis"](0.40)),
        "III": tuple(float(x) for x in plt.colormaps["viridis"](0.65)),
        "IV": tuple(float(x) for x in plt.colormaps["viridis"](0.90)),
    }


def _scatter_panel(
    ax: plt.Axes,
    emb: np.ndarray,
    i: int,
    j: int,
    hue: np.ndarray,
    title: str,
    legend_title: str,
    stage_mode: bool,
) -> None:
    hue = hue.astype(str)
    if stage_mode:
        color_map = _stage_colors()
        order = ["I", "II", "III", "IV"]
        rank = {s: k for k, s in enumerate(order)}
        uniques = sorted(set(hue.tolist()), key=lambda s: rank.get(s, 999))
        for lab in uniques:
            m = hue == lab
            c = color_map.get(lab, (0.5, 0.5, 0.5, 1.0))
            ax.scatter(emb[m, i], emb[m, j], s=4, c=[c], alpha=0.75, linewidths=0, label=str(lab))
    else:
        cmap, uniques = _categorical_colors(hue)
        for lab in uniques:
            m = hue == lab
            c = cmap[lab]
            ax.scatter(emb[m, i], emb[m, j], s=4, c=[c], alpha=0.75, linewidths=0, label=str(lab))
    ax.set_xlabel(f"t-SNE {i + 1}")
    ax.set_ylabel(f"t-SNE {j + 1}")
    ax.set_title(title)
    leg = ax.legend(
        title=legend_title,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0.0,
        fontsize=7,
        title_fontsize=8,
        frameon=False,
        markerscale=2.0,
    )
    for lh in leg.legend_handles:
        lh.set_alpha(1.0)


def main() -> None:
    ap = argparse.ArgumentParser(description="Figure 1B: 4D t-SNE on stage-sample lncRNA matrix → four PNGs.")
    ap.add_argument("--matrix-csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--out-dir", type=Path, default=FIGURES)
    ap.add_argument(
        "--n-pca",
        type=int,
        default=50,
        metavar="P",
        help="PCA components before t-SNE (0 = skip PCA; slow / unstable on full gene count).",
    )
    ap.add_argument("--perplexity", type=float, default=30.0)
    ap.add_argument("--random-state", type=int, default=0)
    args = ap.parse_args()

    if not args.matrix_csv.exists():
        raise SystemExit(f"Missing {args.matrix_csv}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading matrix (this can take a minute)...")
    df = pd.read_csv(args.matrix_csv)
    miss = [c for c in META_COLS if c not in df.columns]
    if miss:
        raise SystemExit(f"{args.matrix_csv}: missing columns {miss}")

    gene_cols = [c for c in df.columns if c not in META_COLS]
    X = df[gene_cols].to_numpy(dtype=np.float64)
    if not np.isfinite(X).all():
        raise SystemExit("Expression matrix contains non-finite values.")

    print(f"Samples: {X.shape[0]}, genes: {X.shape[1]}")
    print("Standardizing genes across samples...")
    Xs = StandardScaler().fit_transform(X).astype(np.float32, copy=False)

    if args.n_pca and args.n_pca > 0:
        n_use = min(args.n_pca, Xs.shape[1], max(1, Xs.shape[0] - 1))
        print(f"PCA to {n_use} components...")
        X_in = PCA(n_components=n_use, random_state=args.random_state, svd_solver="randomized").fit_transform(Xs)
        X_in = X_in.astype(np.float32, copy=False)
    else:
        X_in = Xs

    print("Fitting 4D t-SNE (OpenTSNE)...")
    tsne = OpenTSNE(
        n_components=4,
        perplexity=args.perplexity,
        random_state=args.random_state,
        n_jobs=-1,
    )
    emb = np.asarray(tsne.fit(X_in), dtype=np.float64)
    print("Embedding shape:", emb.shape)

    cancer = df["cancer_type"].astype(str).to_numpy()
    stage = df["stage"].astype(str).to_numpy()

    base = "fig1b_tsne_stage_lncrna_samples"

    def save_pair(dim_i: int, dim_j: int, suffix: str) -> Path:
        fig, ax = plt.subplots(figsize=(7.0, 6.0))
        if suffix.endswith("cancer_type"):
            _scatter_panel(
                ax,
                emb,
                dim_i,
                dim_j,
                cancer,
                f"Fig 1B — t-SNE (dims {dim_i + 1} vs {dim_j + 1}) — cancer type",
                "cancer_type",
                stage_mode=False,
            )
        else:
            _scatter_panel(
                ax,
                emb,
                dim_i,
                dim_j,
                stage,
                f"Fig 1B — t-SNE (dims {dim_i + 1} vs {dim_j + 1}) — AJCC stage",
                "stage",
                stage_mode=True,
            )
        fig.patch.set_facecolor("white")
        ax.set_facecolor("#fafafa")
        for spine in ax.spines.values():
            spine.set_color("#bbbbbb")
        ax.grid(False)
        path = args.out_dir / f"{base}_{suffix}.png"
        fig.tight_layout()
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return path

    p1 = save_pair(0, 1, "dims12_cancer_type")
    p2 = save_pair(0, 1, "dims12_stage")
    p3 = save_pair(2, 3, "dims34_cancer_type")
    p4 = save_pair(2, 3, "dims34_stage")
    print("Wrote:")
    for p in (p1, p2, p3, p4):
        print(" ", p)


if __name__ == "__main__":
    main()
