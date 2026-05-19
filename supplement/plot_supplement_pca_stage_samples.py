"""
Supplementary sample scatter: **PCA** of standardized lncRNA expression (same matrix as Fig 1B).

Four panels (same layout as ``manuscript/plot_figure1b_tsne_stage_lncrna.py``):

- **PC1 vs PC2**, coloured by **cancer_type**
- **PC1 vs PC2**, coloured by **AJCC stage**
- **PC3 vs PC4**, coloured by **cancer_type**
- **PC3 vs PC4**, coloured by **stage**

PCA uses ``sklearn.decomposition.PCA(..., svd_solver='randomized')`` on the full
sample × gene matrix after per-gene standardization (no gene-level PCA truncation
before PCA, unlike the default Fig 1B t-SNE pipeline unless you match it with ``--n-pca``).

Default cohort: ``primary_exp_stage_lncRNA.csv`` with cancer types that have **>100**
samples (same rule as Fig 1B / DE; ``--min-samples-per-cancer 0`` for all types).

Default PNG output directory: ``figures/supplementary/pca/`` (see ``repo_paths.FIGURES_SUPPLEMENTARY_PCA``).
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = _REPO / "manuscript"
for _p in (str(_REPO), str(_REPO / "scripts"), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from repo_paths import DATA, FIGURES, FIGURES_SUPPLEMENTARY_PCA, REPO_ROOT
from figure_export import add_publication_args, save_figure_bundle
from tcga_cohort_filters import MIN_SAMPLES_CANCER, filter_samples_min_cancer_count

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

DEFAULT_CSV = DATA / "primary_exp_stage_lncRNA.csv"
META_COLS = ("sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage")


def _categorical_colors(labels: np.ndarray) -> tuple[dict[str, tuple], list[str]]:
    uniques = sorted({str(x) for x in labels})
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
    x: np.ndarray,
    y: np.ndarray,
    hue: np.ndarray,
    title: str,
    legend_title: str,
    *,
    xlab: str,
    ylab: str,
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
            ax.scatter(x[m], y[m], s=4, c=[c], alpha=0.75, linewidths=0, label=str(lab))
    else:
        cmap, uniques = _categorical_colors(hue)
        for lab in uniques:
            m = hue == lab
            c = cmap[lab]
            ax.scatter(x[m], y[m], s=4, c=[c], alpha=0.75, linewidths=0, label=str(lab))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
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
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=FIGURES_SUPPLEMENTARY_PCA,
        help=f"PNG output directory (default: {FIGURES_SUPPLEMENTARY_PCA.relative_to(REPO_ROOT)}).",
    )
    ap.add_argument(
        "--n-pca",
        type=int,
        default=0,
        metavar="P",
        help="If >0, run gene-level PCA to P components before sample PCA (matches Fig 1B truncation). "
        "Default 0 = sample PCA on all standardized genes (randomized SVD).",
    )
    ap.add_argument("--random-state", type=int, default=0)
    ap.add_argument(
        "--out-prefix",
        type=str,
        default="figS_pca_stage_lncrna_samples",
        help="Output basename prefix for four PNGs.",
    )
    ap.add_argument(
        "--min-samples-per-cancer",
        type=int,
        default=MIN_SAMPLES_CANCER,
        metavar="N",
        help=f"Keep only cancer types with >N samples (default {MIN_SAMPLES_CANCER}; 0 = all types).",
    )
    add_publication_args(ap)
    args = ap.parse_args()

    if not args.matrix_csv.exists():
        raise SystemExit(f"Missing {args.matrix_csv}")

    df = pd.read_csv(args.matrix_csv)
    miss = [c for c in META_COLS if c not in df.columns]
    if miss:
        raise SystemExit(f"{args.matrix_csv}: missing columns {miss}")

    n_before = len(df)
    df, kept_cancers = filter_samples_min_cancer_count(
        df, min_samples_per_cancer=args.min_samples_per_cancer
    )
    if args.min_samples_per_cancer > 0:
        print(
            f"Cohort filter (> {args.min_samples_per_cancer} samples per cancer type): "
            f"{n_before} -> {len(df)} samples; {len(kept_cancers)} cancer types"
        )
    else:
        print(f"No per-cancer sample filter ({len(df)} samples)")

    gene_cols = [c for c in df.columns if c not in META_COLS]
    X = df[gene_cols].to_numpy(dtype=np.float64)
    if not np.isfinite(X).all():
        raise SystemExit("Expression matrix contains non-finite values.")

    print(f"Samples: {X.shape[0]}, genes: {X.shape[1]}")
    Xs = StandardScaler().fit_transform(X).astype(np.float32, copy=False)

    if args.n_pca and args.n_pca > 0:
        n_use = min(args.n_pca, Xs.shape[1], max(1, Xs.shape[0] - 1))
        print(f"Gene-level PCA to {n_use} components (pre-sample-PCA)...")
        X_work = PCA(n_components=n_use, random_state=args.random_state, svd_solver="randomized").fit_transform(
            Xs
        ).astype(np.float32, copy=False)
    else:
        X_work = Xs

    n_pc = min(4, X_work.shape[0] - 1, X_work.shape[1])
    if n_pc < 4:
        raise SystemExit(f"Need at least 4 PCA components; got n_pc={n_pc} (samples={X_work.shape[0]}).")
    print(f"Sample PCA (n_components=4) on matrix shape {X_work.shape}...")
    pca4 = PCA(n_components=4, random_state=args.random_state, svd_solver="randomized")
    Z = pca4.fit_transform(X_work)
    ev = pca4.explained_variance_ratio_
    print(
        "Explained variance ratio (PC1–PC4): "
        + ", ".join(f"{100 * float(ev[i]):.2f}%" for i in range(4))
    )

    cancer = df["cancer_type"].astype(str).to_numpy()
    stage = df["stage"].astype(str).to_numpy()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    base = args.out_prefix

    def save_pair(pc_i: int, pc_j: int, suffix: str) -> Path:
        xi, xj = pc_i - 1, pc_j - 1
        fig, ax = plt.subplots(figsize=(7.0, 6.0))
        if suffix.endswith("cancer_type"):
            _scatter_panel(
                ax,
                Z[:, xi],
                Z[:, xj],
                cancer,
                f"Supplement — PCA (PC{pc_i} vs PC{pc_j}) — cancer type",
                "cancer_type",
                xlab=f"PC{pc_i} ({100 * float(ev[xi]):.1f}% var.)",
                ylab=f"PC{pc_j} ({100 * float(ev[xj]):.1f}% var.)",
                stage_mode=False,
            )
        else:
            _scatter_panel(
                ax,
                Z[:, xi],
                Z[:, xj],
                stage,
                f"Supplement — PCA (PC{pc_i} vs PC{pc_j}) — AJCC stage",
                "stage",
                xlab=f"PC{pc_i} ({100 * float(ev[xi]):.1f}% var.)",
                ylab=f"PC{pc_j} ({100 * float(ev[xj]):.1f}% var.)",
                stage_mode=True,
            )
        fig.patch.set_facecolor("white")
        ax.set_facecolor("#fafafa")
        for spine in ax.spines.values():
            spine.set_color("#bbbbbb")
        ax.grid(False)
        path = args.out_dir / f"{base}_pc{pc_i}_pc{pc_j}_{suffix}.png"
        fig.tight_layout()
        save_figure_bundle(
            fig,
            path,
            png_dpi=200,
            publication_dir=args.publication_dir,
            publication_tiff_kind=args.publication_tiff_kind,
            figures_root=FIGURES,
            bbox_inches="tight",
        )
        plt.close(fig)
        return path

    p1 = save_pair(1, 2, "cancer_type")
    p2 = save_pair(1, 2, "stage")
    p3 = save_pair(3, 4, "cancer_type")
    p4 = save_pair(3, 4, "stage")
    print("Wrote:")
    for p in (p1, p2, p3, p4):
        print(" ", p)


if __name__ == "__main__":
    main()
