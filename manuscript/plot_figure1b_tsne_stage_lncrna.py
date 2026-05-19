"""
Figure 1B — sample embedding on TCGA primary tumor lncRNA expression (same matrix as DE pipeline).

Default (**``--embedding sklearn2_pca34``**): **sklearn t-SNE** in **2D** only (Barnes–Hut on ``X_in``;
``X_in`` may be gene-level PCA–truncated first via ``--n-pca``). Writes separate PNGs coloured by
**cancer_type** and **AJCC stage**, plus a **two-panel** figure (**A** = cancer type, **B** = stage).
There are **no** PC3/PC4 or ``_dims34_`` files in this mode.

Optional **``--embedding opentsne4``**: **4D OpenTSNE** — separate and combined panels for dims 1–2 and 3–4.
OpenTSNE may emit a **FutureWarning** about Barnes–Hut in >3 dimensions; that path filters it.

Loads ``data/primary_exp_stage_lncRNA.csv`` (stage primary-tumor matrix only; not the
metastasis matrix). By default keeps samples from cancer types with **>100** samples in that
matrix (same rule as ``tr_lncrna_de_analysis`` / Fig 2 bars; ``--min-samples-per-cancer 0`` for all types).
Standardizes genes across samples, optionally PCA-prelimits dimensionality, then fits the chosen embedding.

Basenames use ``--filename-prefix`` (default ``fig1b_tsne_stage_lncrna_samples``) and ``--out-dir``
(default ``figures/``) so alternate embeddings can live under ``figures/supplementary/embedding/``
without clobbering canonical files.
"""
from __future__ import annotations

from pathlib import Path
import sys
import warnings

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import DATA, FIGURES, REPO_ROOT
from figure_export import add_publication_args, save_figure_bundle
from tcga_cohort_filters import MIN_SAMPLES_CANCER, filter_samples_min_cancer_count

ROOT = REPO_ROOT


import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE as SklearnTSNE
from sklearn.preprocessing import StandardScaler

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


def _style_embedding_figure(fig: plt.Figure, axes: list[plt.Axes]) -> None:
    fig.patch.set_facecolor("white")
    for ax in axes:
        ax.set_facecolor("#fafafa")
        for spine in ax.spines.values():
            spine.set_color("#bbbbbb")
        ax.grid(False)


def _scatter_panel(
    ax: plt.Axes,
    emb: np.ndarray,
    i: int,
    j: int,
    hue: np.ndarray,
    title: str,
    legend_title: str,
    stage_mode: bool,
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    panel_label: str | None = None,
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
    ax.set_xlabel(xlabel if xlabel is not None else f"t-SNE {i + 1}")
    ax.set_ylabel(ylabel if ylabel is not None else f"t-SNE {j + 1}")
    if title:
        ax.set_title(title)
    if panel_label:
        # Outside and left of the y-axis label (t-SNE 2), not over the scatter.
        ax.text(
            -0.11,
            1.0,
            panel_label,
            transform=ax.transAxes,
            fontsize=14,
            fontweight=400,
            va="bottom",
            ha="right",
            color="black",
            clip_on=False,
        )
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
    ap = argparse.ArgumentParser(
        description="Figure 1B: sample embedding on stage lncRNA matrix → two PNGs (sklearn default) or four (opentsne4)."
    )
    ap.add_argument("--matrix-csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--out-dir", type=Path, default=FIGURES)
    ap.add_argument(
        "--filename-prefix",
        type=str,
        default="fig1b_tsne_stage_lncrna_samples",
        metavar="STEM",
        help="Basename prefix for PNGs (before _dims12_* and, with opentsne4, _dims34_*).",
    )
    ap.add_argument(
        "--n-pca",
        type=int,
        default=50,
        metavar="P",
        help="PCA components before t-SNE (0 = skip PCA; slow / unstable on full gene count).",
    )
    ap.add_argument("--perplexity", type=float, default=30.0)
    ap.add_argument("--random-state", type=int, default=0)
    ap.add_argument(
        "--embedding",
        choices=("sklearn2_pca34", "opentsne4"),
        default="sklearn2_pca34",
        help="sklearn2_pca34: 2D sklearn t-SNE only (two PNGs: dims 1–2). "
        "opentsne4: 4D OpenTSNE (four PNGs: t-SNE dims 1–2 and 3–4).",
    )
    ap.add_argument(
        "--min-samples-per-cancer",
        type=int,
        default=MIN_SAMPLES_CANCER,
        metavar="N",
        help=f"Keep only cancer types with >N samples in the matrix (default {MIN_SAMPLES_CANCER}, "
        "aligned with DE/peptide-fraction cohort). Use 0 for all cancer types.",
    )
    ap.add_argument(
        "--no-combined-ab-panel",
        action="store_true",
        help="Skip the two-panel figure (A: cancer type, B: AJCC stage) in one PNG.",
    )
    add_publication_args(ap)
    args = ap.parse_args()

    if not args.matrix_csv.exists():
        raise SystemExit(f"Missing {args.matrix_csv}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading matrix (this can take a minute)...")
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
            f"{n_before} -> {len(df)} samples; {len(kept_cancers)} cancer types: "
            f"{', '.join(kept_cancers)}"
        )
    else:
        print(f"No per-cancer sample filter ({len(df)} samples, all cancer types in matrix)")

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

    n_s = X_in.shape[0]
    if args.embedding == "sklearn2_pca34":
        perplexity = float(min(args.perplexity, max(5.0, (n_s - 1) / 3.0)))
        print(f"Fitting 2D t-SNE (sklearn, perplexity={perplexity:.3g}, Barnes–Hut)...")
        tsne2 = SklearnTSNE(
            n_components=2,
            perplexity=perplexity,
            random_state=args.random_state,
            max_iter=1000,
            init="pca",
            method="barnes_hut",
            learning_rate="auto",
        )
        emb12 = np.asarray(tsne2.fit_transform(X_in), dtype=np.float64)
        emb = emb12.astype(np.float64, copy=False)
    else:
        try:
            from openTSNE import TSNE as OpenTSNE
        except ImportError as e:  # pragma: no cover
            raise SystemExit("Install openTSNE for --embedding opentsne4: pip install opentsne") from e
        print("Fitting 4D t-SNE (OpenTSNE)...")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*BH t-SNE for >3 dimensions.*",
                category=FutureWarning,
            )
            tsne = OpenTSNE(
                n_components=4,
                perplexity=args.perplexity,
                random_state=args.random_state,
                n_jobs=-1,
            )
            emb = np.asarray(tsne.fit(X_in), dtype=np.float64)
        print("Embedding shape:", emb.shape)

    if args.embedding == "sklearn2_pca34":
        print("Embedding shape:", emb.shape, "(2D t-SNE only; no dims 3–4 PNGs — use --embedding opentsne4 for four panels)")

    cancer = df["cancer_type"].astype(str).to_numpy()
    stage = df["stage"].astype(str).to_numpy()

    base = args.filename_prefix.strip() or "fig1b_tsne_stage_lncrna_samples"

    def save_pair(
        dim_i: int,
        dim_j: int,
        suffix: str,
        *,
        dim_tag: str,
        xlabel: str | None = None,
        ylabel: str | None = None,
    ) -> Path:
        fig, ax = plt.subplots(figsize=(7.0, 6.0))
        if suffix.endswith("cancer_type"):
            _scatter_panel(
                ax,
                emb,
                dim_i,
                dim_j,
                cancer,
                f"t-SNE ({dim_tag}) — cancer type",
                "cancer_type",
                stage_mode=False,
                xlabel=xlabel,
                ylabel=ylabel,
            )
        else:
            _scatter_panel(
                ax,
                emb,
                dim_i,
                dim_j,
                stage,
                f"t-SNE ({dim_tag}) — AJCC stage",
                "stage",
                stage_mode=True,
                xlabel=xlabel,
                ylabel=ylabel,
            )
        _style_embedding_figure(fig, [ax])
        path = args.out_dir / f"{base}_{suffix}.png"
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

    def save_combined_ab(
        dim_i: int,
        dim_j: int,
        suffix: str,
        *,
        dim_tag: str,
    ) -> Path:
        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.2, 5.4))
        _scatter_panel(
            ax_a,
            emb,
            dim_i,
            dim_j,
            cancer,
            "",
            "cancer_type",
            stage_mode=False,
            panel_label="A",
        )
        _scatter_panel(
            ax_b,
            emb,
            dim_i,
            dim_j,
            stage,
            "",
            "stage",
            stage_mode=True,
            panel_label="B",
        )
        ax_b.set_ylabel("")
        _style_embedding_figure(fig, [ax_a, ax_b])
        path = args.out_dir / f"{base}_{suffix}.png"
        fig.tight_layout(w_pad=0.15, pad=0.4)
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

    paths_out: list[Path] = []
    if args.embedding == "sklearn2_pca34":
        paths_out.append(save_pair(0, 1, "dims12_cancer_type", dim_tag="dims 1 vs 2"))
        paths_out.append(save_pair(0, 1, "dims12_stage", dim_tag="dims 1 vs 2"))
        if not args.no_combined_ab_panel:
            paths_out.append(save_combined_ab(0, 1, "dims12_panels_AB", dim_tag="dims 1 vs 2"))
    else:
        paths_out.append(save_pair(0, 1, "dims12_cancer_type", dim_tag="dims 1 vs 2"))
        paths_out.append(save_pair(0, 1, "dims12_stage", dim_tag="dims 1 vs 2"))
        if not args.no_combined_ab_panel:
            paths_out.append(save_combined_ab(0, 1, "dims12_panels_AB", dim_tag="dims 1 vs 2"))
        paths_out.append(save_pair(2, 3, "dims34_cancer_type", dim_tag="dims 3 vs 4"))
        paths_out.append(save_pair(2, 3, "dims34_stage", dim_tag="dims 3 vs 4"))
        if not args.no_combined_ab_panel:
            paths_out.append(save_combined_ab(2, 3, "dims34_panels_AB", dim_tag="dims 3 vs 4"))
    print("Wrote:")
    for p in paths_out:
        print(" ", p)


if __name__ == "__main__":
    main()
