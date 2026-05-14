"""
Volcano plot: overlapping 2-mer (dipeptide) counts in lncRNA-derived peptides vs
``data/known_proteins.fasta``.

- **log2 fold change (x):** same smoothing as the heatmaps in ``plot_dipeptide_mp_figure2.py``:
  Dirichlet-style ``(count + 0.5)`` over all 400 cells, then ``log2(f_lnc / f_proteome)``.
- **Significance (y):** two-sided Fisher exact test per dipeptide on a 2×2 table
  (this dipeptide vs all other dipeptides × lnc vs proteome); **Benjamini–Hochberg FDR**
  across **400** tests. **q-value** = FDR-adjusted *p*.

**Point colors (default thresholds; color-blind-friendly palette):**
- **Vermillion / orange:** ``log2FC > 1`` and ``q < 0.05`` (enriched in lnc peptides).
- **Blue:** ``log2FC < -1`` and ``q < 0.05`` (depleted).
- **Yellow:** ``q < 0.05`` and ``-1 <= log2FC <= 1`` (significant, small effect).
- **Gray:** ``q >= 0.05``.

Writes PNG + CSV + short report under ``figures/<peptide_set>/`` by default. For the
default TCGA-matrix run (no custom FASTA / ``--out-dir``, not ``--canonical-only``), a copy
is also written to ``figures/fig3b.png``.

**Manuscript Figure 3B:** default is TCGA-matrix FASTA (``smprot_tcga_filtered_peptides.faa``),
same default peptide universe as **Figure 3A** in ``plot_aa_frequency_tcga_vs_proteome.py``.
Use ``--peptide-set all_smprot_filtered`` on **both** 3A and 3B scripts for the full SmProt-filtered alternate.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_MS), str(_REPO / "supplement")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES
from figure_export import add_publication_args, save_figure_bundle

ROOT = REPO_ROOT


import argparse
import shutil

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import stats

import plot_dipeptide_mp_figure2 as dp

import figure_palettes as pal

DATA = ROOT / "data"
ALL_SMPROT_FILTERED_FA = DATA / "smprot_all_filtered_peptides.faa"
FIGURES = ROOT / "figures"


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, dtype=np.float64)
    ranked[order] = np.arange(1, n + 1)
    q = p * n / ranked
    adj_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    adj = np.empty(n, dtype=np.float64)
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    return adj


def classify(log2fc: np.ndarray, q: np.ndarray, fc_hi: float, fc_lo: float, q_cut: float) -> np.ndarray:
    """
    Return string category per point: enriched / depleted / sig_small / ns.
    Boundaries: red/blue require |fc| > threshold; yellow uses inclusive |fc| <= 1.
    """
    out = np.full(log2fc.shape, "ns", dtype=object)
    sig = q < q_cut
    out[sig & (log2fc > fc_hi)] = "enriched"
    out[sig & (log2fc < fc_lo)] = "depleted"
    out[sig & (log2fc >= fc_lo) & (log2fc <= fc_hi)] = "sig_small"
    return out


def dipeptide_labels() -> list[str]:
    return [f"{a}{b}" for a in dp.AA20 for b in dp.AA20]


def main() -> None:
    ap = argparse.ArgumentParser(description="Dipeptide volcano: lnc MPs vs proteome (Fisher + FDR).")
    ap.add_argument(
        "--peptide-set",
        choices=("tcga_matrix", "all_smprot_filtered"),
        default="tcga_matrix",
        help=(
            "tcga_matrix: data/smprot_tcga_filtered_peptides.faa (default). "
            "all_smprot_filtered: data/smprot_all_filtered_peptides.faa (export from smprot_filtered.tsv)."
        ),
    )
    ap.add_argument(
        "--peptide-fa",
        type=Path,
        default=None,
        help="Override peptide FASTA (ignores --peptide-set).",
    )
    ap.add_argument(
        "--proteome-fa",
        type=Path,
        default=dp.PROTEOME_FA,
        help="Reference proteome FASTA (default: data/known_proteins.fasta).",
    )
    ap.add_argument(
        "--canonical-only",
        action="store_true",
        help="Restrict peptides to GeneSymbol in canonical Tr list (heatmap panel B style).",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: figures/<peptide_set>/ when --peptide-fa is unset).",
    )
    ap.add_argument(
        "--fc-threshold",
        type=float,
        default=1.0,
        help="Absolute log2 fold change threshold for enriched/depleted (default 1.0).",
    )
    ap.add_argument(
        "--q-cutoff",
        type=float,
        default=0.05,
        help="FDR q-value cutoff for significance (default 0.05).",
    )
    ap.add_argument(
        "--annotate-n",
        type=int,
        default=22,
        help="Label up to this many extra points (top by -log10 q) beyond mandatory extremes.",
    )
    add_publication_args(ap)
    args = ap.parse_args()

    if args.peptide_fa is not None:
        peptide_fa = args.peptide_fa
    elif args.peptide_set == "tcga_matrix":
        peptide_fa = dp.MP_FAA
    else:
        peptide_fa = ALL_SMPROT_FILTERED_FA

    if not peptide_fa.exists():
        hint = ""
        if peptide_fa == ALL_SMPROT_FILTERED_FA:
            hint = (
                "\nCreate it with:\n  python export_tcga_filtered_peptides_fasta.py "
                "--peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa"
            )
        raise SystemExit(f"Missing peptide FASTA: {peptide_fa}{hint}")
    if not args.proteome_fa.exists():
        raise SystemExit(f"Missing proteome FASTA: {args.proteome_fa}")

    if args.out_dir is not None:
        out_dir = args.out_dir
    elif args.peptide_fa is None:
        out_dir = FIGURES / args.peptide_set
    else:
        out_dir = FIGURES / args.peptide_fa.stem

    tr_genes = dp.load_tr_genes() if args.canonical_only else None
    mp_recs = dp.parse_mp_fasta(peptide_fa)

    print("Counting proteome dipeptides (may take a minute)...")
    c_ref, _, di_ref = dp.count_proteome(args.proteome_fa)
    c_mp, _, di_mp = dp.count_mp_subset(mp_recs, gene_filter=tr_genes)

    if di_mp <= 0 or di_ref <= 0:
        raise SystemExit("No overlapping dipeptides counted in one or both sets.")

    f_mp = dp.freq_matrix(c_mp)
    f_ref = dp.freq_matrix(c_ref)
    log2fc = dp.log2_ratio_mat(f_mp, f_ref).ravel()

    p_raw = np.zeros(400, dtype=np.float64)
    for k in range(400):
        i, j = k // 20, k % 20
        a = int(c_mp[i, j])
        b = int(di_mp - c_mp[i, j])
        c_ = int(c_ref[i, j])
        d = int(di_ref - c_ref[i, j])
        _, p_two = stats.fisher_exact([[a, b], [c_, d]], alternative="two-sided")
        p_raw[k] = p_two

    q = benjamini_hochberg(p_raw)
    fc_hi = args.fc_threshold
    fc_lo = -args.fc_threshold
    cats = classify(log2fc, q, fc_hi=fc_hi, fc_lo=fc_lo, q_cut=args.q_cutoff)

    labels = dipeptide_labels()
    df = pd.DataFrame(
        {
            "dipeptide": labels,
            "count_lnc": c_mp.ravel().astype(np.int64),
            "count_proteome": c_ref.ravel().astype(np.int64),
            "log2_fold_change_lnc_over_proteome": log2fc,
            "fisher_p_two_sided": p_raw,
            "q_value_fdr_bh": q,
            "category": cats,
            "peptide_fasta": str(peptide_fa),
            "peptide_set": args.peptide_set if args.peptide_fa is None else "custom_fasta",
        }
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    if args.canonical_only:
        suffix = (
            "_all_smprot_filtered_tr_canonical"
            if args.peptide_set == "all_smprot_filtered" and args.peptide_fa is None
            else "_tr_canonical"
        )
    else:
        suffix = (
            "_all_smprot_filtered"
            if args.peptide_set == "all_smprot_filtered" and args.peptide_fa is None
            else "_tcga_filtered"
        )
    csv_path = out_dir / f"dipeptide_volcano_lnc_vs_proteome_stats{suffix}.csv"
    df.to_csv(csv_path, index=False)

    # y-axis: cap display for q=0
    y = -np.log10(np.maximum(q, 1e-300))

    colors = np.full(400, pal.VOLCANO_NS, dtype=object)
    colors[cats == "enriched"] = pal.VOLCANO_ENRICHED
    colors[cats == "depleted"] = pal.VOLCANO_DEPLETED
    colors[cats == "sig_small"] = pal.VOLCANO_SMALL_EFFECT

    fig, ax = plt.subplots(figsize=(8.5, 6.2), dpi=150)
    ax.scatter(log2fc, y, c=list(colors), s=14, alpha=0.85, edgecolors="none")

    q_line = -np.log10(args.q_cutoff)
    ax.axhline(q_line, color="0.45", linestyle="--", linewidth=1.0, zorder=0)
    ax.axvline(fc_hi, color=pal.VOLCANO_ENRICHED, linestyle="--", linewidth=0.9, alpha=0.7, zorder=0)
    ax.axvline(fc_lo, color=pal.VOLCANO_DEPLETED, linestyle="--", linewidth=0.9, alpha=0.7, zorder=0)

    n_enr = int((cats == "enriched").sum())
    n_dep = int((cats == "depleted").sum())
    n_yel = int((cats == "sig_small").sum())
    n_gray = int((cats == "ns").sum())

    patch_kw = dict(edgecolor="0.35", linewidth=0.4)
    legend_handles = [
        Patch(facecolor=pal.VOLCANO_ENRICHED, label=f"Enriched (n={n_enr})", **patch_kw),
        Patch(facecolor=pal.VOLCANO_DEPLETED, label=f"Depleted (n={n_dep})", **patch_kw),
        Patch(facecolor=pal.VOLCANO_SMALL_EFFECT, label=f"Sig. small effect (n={n_yel})", **patch_kw),
        Patch(facecolor=pal.VOLCANO_NS, label=f"Not significant (n={n_gray})", **patch_kw),
    ]
    # Legend to the right of the axes (compact); full-width axes via tight_layout + bbox_extra_artists on save.
    leg = ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
        fontsize=8,
        handlelength=0.85,
        handleheight=0.55,
        borderpad=0.28,
        labelspacing=0.18,
    )

    ax.text(
        0.02,
        0.98,
        f"Thresholds\nq < {args.q_cutoff} (FDR)\n|log2FC| > {fc_hi}",
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.88, edgecolor="0.7"),
    )

    # Labels: top significant depleted / enriched by q, plus a few global extremes
    annotate: set[int] = set()
    for mask in (cats == "depleted", cats == "enriched"):
        idxs = np.where(mask)[0]
        if idxs.size:
            order = idxs[np.argsort(-y[idxs])]
            annotate.update(order[:10].tolist())
    glob_order = np.argsort(-y)
    for k in glob_order[: max(12, args.annotate_n)]:
        annotate.add(int(k))

    for k in sorted(annotate):
        ax.annotate(
            labels[k],
            (log2fc[k], y[k]),
            fontsize=7,
            alpha=0.92,
            xytext=(4, 3),
            textcoords="offset points",
            ha="left",
        )

    ax.set_xlabel("log2 fold change (lnc / proteome)")
    ax.set_ylabel("-log10(q-value)")
    if args.canonical_only:
        title_set = "Tr-lncRNA MPs (canonical)"
    elif args.peptide_set == "all_smprot_filtered" and args.peptide_fa is None:
        title_set = "lncRNA peptides (all SmProt-filtered)"
    else:
        title_set = "lncRNA peptides (TCGA-matrix)"
    ax.set_title(f"Volcano plot: 2-mer enrichment ({title_set} vs proteome)")
    ax.set_xlim(left=min(-3.2, float(np.nanmin(log2fc)) - 0.15), right=max(2.2, float(np.nanmax(log2fc)) + 0.15))
    ymax = float(np.nanmax(y[np.isfinite(y)]))
    ax.set_ylim(0, ymax * 1.05 + 0.05)
    fig.tight_layout()

    png_path = out_dir / f"dipeptide_volcano_lnc_vs_proteome{suffix}.png"
    save_figure_bundle(
        fig,
        png_path,
        png_dpi=150,
        publication_dir=args.publication_dir,
        publication_tiff_kind=args.publication_tiff_kind,
        figures_root=FIGURES,
        bbox_inches="tight",
        bbox_extra_artists=(leg,),
        pad_inches=0.12,
    )
    plt.close(fig)

    if args.peptide_fa is None and args.peptide_set == "tcga_matrix" and args.out_dir is None and not args.canonical_only:
        fig3b_root = FIGURES / "fig3b.png"
        shutil.copy2(png_path, fig3b_root)
        print(f"Copied to {fig3b_root} (Fig 3B)")

    report = out_dir / f"dipeptide_volcano_lnc_vs_proteome_report{suffix}.txt"
    report.write_text(
        "\n".join(
            [
                f"peptide_fasta: {peptide_fa}",
                f"peptide_set_arg: {args.peptide_set}",
                f"proteome_fasta: {args.proteome_fa}",
                f"canonical_only: {args.canonical_only}",
                f"n_dipeptides_lnc: {int(di_mp):,}",
                f"n_dipeptides_proteome: {int(di_ref):,}",
                f"log2FC: freq_matrix(+0.5) per corpus then log2(lnc/proteome), same as plot_dipeptide_mp_figure2.py",
                f"test: Fisher two-sided 2x2 per dipeptide; FDR-BH across 400",
                f"q_cutoff: {args.q_cutoff}; fc_threshold: {fc_hi}",
                f"csv: {csv_path}",
                f"png: {png_path}",
            ]
            + (
                [f"fig3b_root_copy: {FIGURES / 'fig3b.png'}"]
                if args.peptide_fa is None
                and args.peptide_set == "tcga_matrix"
                and args.out_dir is None
                and not args.canonical_only
                else []
            )
        )
        + "\n",
        encoding="utf-8",
    )
    print(report.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
