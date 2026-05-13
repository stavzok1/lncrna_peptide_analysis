#!/usr/bin/env python3
"""
Figure 5a–5c from merged NetMHC + IEDB long tables (``*_with_iedb.tsv``).

**Why this script vs wide-XLS Fig 5?** The **wide** NetMHCpan scripts (``plot_netmhc_epitopes_vs_hla_frequency.py``,
``plot_figure5b_epitope_sharing_across_alleles.py``) use **IC50-from-BA only** on ``*.xls``.
This triple script uses the **merged** ``*_with_iedb.tsv`` rows so SB can combine **IEDB**
immunogenicity / processing with **NetMHC EL / IC50** (same role as manuscript 5A–5C, richer gate).
Plot styling is aligned with those wide-XLS figures (sizes, colours) so panels are visually comparable.

**SB** (see ``scripts/netmhc_sb_core.py`` ``FIG5_IEDB_*`` and ``--sb-mode``):

  - ``full`` (default): IEDB immuno + processing + NetMHC ``EL_rank`` + IC50 (IEDB column if present, else ``BA_score``).
  - ``no_ic50``: same but **no** IC50 / BA binding gate (IEDB + EL only).
  - ``ic50_only``: **only** IC50 / BA binding at ``--ic50-max-nm`` (no IEDB immuno/proc/EL gates).

5a: allele population frequency vs **SB epitopes** (SB prediction rows per allele) by default (optional unique 9-mers).

5b / 5c: epitope-sharing histogram — **same layout** as ``plot_figure5b_epitope_sharing_across_alleles.py``,
using **color-blind-friendly** sig vs coding fills from ``figure_palettes``; x = 1 … N panel alleles (default 27), zero-filled bins.

Histogram **y** defaults to **``n_sb_rows_in_bin``** (axis label **SB epitopes per allele**; bar height is SB epitope×allele mass per bin); CSV
still includes **``n_epitopes``** (unique 9-mers per bin) and **``total_sb_row_instances``** (total SB epitopes × alleles).
Use **``--sharing-y-metric unique``** for the unique-9-mer bar height.

Allele frequencies default to ``fig5a_epitopes_vs_allele_frequency_ic50_sb.csv`` (5A only; skipped for ``--panels c``).

Use ``--panels c`` with a second ``--coding-tsv`` / ``--output-stem`` to emit **5C** for an alternate merged coding cohort without duplicating canonical **5B** (Tr-lncRNA-MPs) outputs.

**Catalog:** this pipeline is **Figure 5** (cohort allele–epitope interplay with IEDB filters);
see ``docs/figure_catalog.md``. **Figure 6** is TTN-AS1 only (`plot_figure6_ttn_as1_allele_coverage.py`).
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
import shutil

import figure_palettes as pal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy.stats import spearmanr

REPO_FIGURES = FIGURES


def _mirror(*paths: Path) -> None:
    REPO_FIGURES.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.is_file():
            shutil.copy2(p, REPO_FIGURES / p.name)


from netmhc_sb_core import (
    FIG5_IEDB_EL_RANK_MAX_DEFAULT,
    FIG5_IEDB_IC50_MAX_NM_DEFAULT,
    FIG5_IEDB_IMM_MIN_DEFAULT,
    FIG5_IEDB_PROC_MIN_DEFAULT,
    ba_score_min_for_ic50_lt,
    pick_iedb_ic50_column,
    sb_mask_fig5_defaults,
    sb_spec_from_mode,
)


def load_allele_frequencies(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    if not {"allele", "allele_frequency"} <= set(df.columns):
        raise SystemExit(f"{path} must have columns allele, allele_frequency")
    return df.set_index("allele")["allele_frequency"].astype(float)


def read_merged(path: Path, usecols: set[str]) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols)


def fig5a_table(sb: pd.DataFrame, freqs: pd.Series) -> pd.DataFrame:
    total_unique = int(sb["Peptide"].nunique())
    rows = []
    for al, g in sb.groupby("allele", sort=False):
        n = g["Peptide"].nunique()
        af = float(freqs.get(al, np.nan))
        rows.append(
            {
                "allele": al,
                "allele_frequency": af,
                "n_unique_epitopes_sb": n,
                "n_sb_row_instances_allele": int(len(g)),
                "total_unique_epitopes_sb": total_unique,
            }
        )
    return pd.DataFrame(rows).sort_values("allele")


def plot_fig5a(
    tab: pd.DataFrame,
    out_png: Path,
    title: str,
    sb_legend: str,
    *,
    y_metric: str,
) -> None:
    t = tab.dropna(subset=["allele_frequency"])
    if len(t) < 2:
        print("[fig5a] skip plot: need >=2 alleles with frequency", file=__import__("sys").stderr)
        return
    y_col = "n_sb_row_instances_allele" if y_metric == "instances" else "n_unique_epitopes_sb"
    y_label = "SB epitopes per allele" if y_metric == "instances" else "Unique 9-mers per allele"
    x = t["allele_frequency"].to_numpy(float)
    y = t[y_col].to_numpy(float)
    rho, p = spearmanr(x, y)

    _ = sb_legend
    fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=150)
    ax.scatter(x, y, s=36, alpha=0.85, color=pal.SCATTER_DEFAULT, edgecolors="none")
    for _, r in t.iterrows():
        ax.annotate(
            str(r["allele"]).replace("HLA-", ""),
            (r["allele_frequency"], r[y_col]),
            fontsize=7,
            xytext=(3, 3),
            textcoords="offset points",
            alpha=0.9,
        )
    ax.set_xlabel("Allele frequency", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    total_u = int(tab["total_unique_epitopes_sb"].iloc[0]) if "total_unique_epitopes_sb" in tab.columns and len(tab) else None
    extra = f"\nTotal unique 9-mers (global): {total_u}" if total_u is not None else ""
    ax.set_title(title + f"\nSpearman ρ = {rho:.3f}, p = {p:.2g}{extra}", fontsize=11)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig5bc_sharing(sb: pd.DataFrame) -> pd.Series:
    """Per Peptide → number of alleles with ≥1 SB row."""
    return sb.groupby("Peptide")["allele"].nunique()


def sharing_histogram_full(
    cnt: pd.Series,
    max_alleles: int,
    *,
    total_sb_row_instances: int | None = None,
) -> pd.DataFrame:
    """
    Integer bins 1 … max_alleles, zero-filled.

    - ``n_epitopes``: **unique** 9-mers with exactly that many SB alleles.
    - ``n_sb_rows_in_bin``: SB epitope×allele **mass** in the bin (= ``k * n_epitopes`` for bin ``k``); Fig 5B/5C
      use this for bar height with y-axis label **SB epitopes per allele**. The cohort total row count is **SB epitopes × alleles**.
    """
    vc = cnt.value_counts()
    full_idx = np.arange(1, max_alleles + 1, dtype=int)
    total_epitopes = int(cnt.shape[0])
    n_epis = [int(vc.get(i, 0)) for i in full_idx]
    rows_in_bin = [int(k) * int(n) for k, n in zip(full_idx, n_epis)]
    inferred = int(sum(rows_in_bin))
    tot_rows = int(total_sb_row_instances) if total_sb_row_instances is not None else inferred
    return pd.DataFrame(
        {
            "n_alleles_per_epitope": full_idx,
            "n_epitopes": n_epis,
            "n_sb_rows_in_bin": rows_in_bin,
            "total_epitopes": total_epitopes,
            "total_sb_row_instances": tot_rows,
        }
    )


def plot_fig5bc_legacy(
    hist: pd.DataFrame,
    out_png: Path,
    title: str,
    sb_legend: str,
    *,
    bar_color: str,
    y_metric: str,
) -> None:
    _ = sb_legend
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=150)
    x = hist["n_alleles_per_epitope"].to_numpy()
    y_col = "n_sb_rows_in_bin" if y_metric == "instances" else "n_epitopes"
    y_label = "SB epitopes per allele" if y_metric == "instances" else "Unique 9-mers per bin"
    h = hist[y_col].to_numpy()
    ax.bar(x, h, width=0.92, color=bar_color, edgecolor="white", linewidth=0.4, align="center")
    ax.set_xticks(x)
    ax.set_xlabel("Number of alleles", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(title, fontsize=11)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax.yaxis.grid(True, linestyle="-", alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0.5, len(x) + 0.5)
    if "total_epitopes" in hist.columns and len(hist):
        total = int(hist["total_epitopes"].iloc[0])
        inst = int(hist["total_sb_row_instances"].iloc[0]) if "total_sb_row_instances" in hist.columns else None
        if y_metric == "instances":
            summary = f"SB epitopes: {int(h.sum())}"
            if inst is not None:
                summary += f" | SB epitopes × alleles: {inst}"
            summary += f" | unique 9-mers: {total}"
        else:
            summary = f"Unique 9-mers: {total}"
            if inst is not None:
                summary += f" | SB epitopes × alleles: {inst}"
        ax.text(
            0.98,
            0.98,
            summary,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=11,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.4", "alpha": 0.92},
        )
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sig-tsv", type=Path, default=Path("data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv"))
    ap.add_argument(
        "--coding-tsv",
        type=Path,
        default=Path("data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv"),
    )
    ap.add_argument(
        "--allele-freq-csv",
        type=Path,
        default=Path("data/netmhc/figures/fig5a_epitopes_vs_allele_frequency_ic50_sb.csv"),
    )
    ap.add_argument("--out-dir", type=Path, default=Path("data/netmhc/figures"))
    ap.add_argument(
        "--imm-min",
        type=float,
        default=FIG5_IEDB_IMM_MIN_DEFAULT,
        help="IEDB immunogenicity: iedb_score > this (default from netmhc_sb_core).",
    )
    ap.add_argument(
        "--proc-min",
        type=float,
        default=FIG5_IEDB_PROC_MIN_DEFAULT,
        help="IEDB processing: iedb_processing_score > this (default from netmhc_sb_core).",
    )
    ap.add_argument(
        "--el-rank-max",
        type=float,
        default=FIG5_IEDB_EL_RANK_MAX_DEFAULT,
        help="EL %%rank cutoff (default from netmhc_sb_core).",
    )
    ap.add_argument(
        "--el-rank-lte",
        action="store_true",
        help="Use EL_rank ≤ --el-rank-max. Default: EL_rank < --el-rank-max.",
    )
    ap.add_argument(
        "--ic50-max-nm",
        type=float,
        default=FIG5_IEDB_IC50_MAX_NM_DEFAULT,
        help="IC50 strict upper bound in nM (default from netmhc_sb_core).",
    )
    ap.add_argument(
        "--sb-mode",
        choices=("full", "no_ic50", "ic50_only"),
        default="full",
        help="SB composition: full IEDB+EL+IC50; no_ic50 (IEDB+EL only); ic50_only (binding/IC50 only).",
    )
    ap.add_argument(
        "--output-stem",
        type=str,
        default="fig5abc_sb_immuno_proc_el_ic50",
        help="Prefix for output CSV/PNG filenames under --out-dir.",
    )
    ap.add_argument(
        "--n-panel-alleles",
        type=int,
        default=27,
        help="Max alleles on x-axis for 5b/5c (legacy full range with zeros).",
    )
    ap.add_argument(
        "--fig5a-y-metric",
        choices=("instances", "unique"),
        default="instances",
        help="5A scatter y-axis: SB epitopes vs unique 9-mers per allele (default instances).",
    )
    ap.add_argument(
        "--sharing-y-metric",
        choices=("instances", "unique"),
        default="instances",
        help="5B/5C histogram bar height: instances (ylabel SB epitopes per allele) vs unique 9-mers per bin (default instances).",
    )
    ap.add_argument(
        "--panels",
        choices=("abc", "c"),
        default="abc",
        help="abc: emit 5A–5C (default). c: coding-cohort epitope sharing (5C) only; uses --coding-tsv.",
    )
    ap.add_argument(
        "--no-repo-mirror",
        action="store_true",
        help="Do not copy 5A–5C PNGs into repo-root figures/.",
    )
    args = ap.parse_args()

    hdr_cod = set(pd.read_csv(args.coding_tsv, sep="\t", nrows=0).columns)
    hdr_sig = (
        set(pd.read_csv(args.sig_tsv, sep="\t", nrows=0).columns) if args.panels == "abc" else set()
    )
    iedb_ic50_col = pick_iedb_ic50_column(sorted(hdr_sig | hdr_cod))

    usecols = {
        "allele",
        "Peptide",
        "EL_rank",
        "iedb_score",
        "iedb_processing_score",
    }
    if iedb_ic50_col:
        usecols.add(iedb_ic50_col)
    else:
        usecols.add("BA_score")

    freqs = load_allele_frequencies(args.allele_freq_csv) if args.panels == "abc" else None
    ba_min = ba_score_min_for_ic50_lt(args.ic50_max_nm)

    sig = read_merged(args.sig_tsv, usecols) if args.panels == "abc" else None
    cod = read_merged(args.coding_tsv, usecols)

    el_cmp = "≤" if args.el_rank_lte else "<"
    spec_leg = sb_spec_from_mode(
        args.sb_mode,
        imm_min=args.imm_min,
        proc_min=args.proc_min,
        el_max=args.el_rank_max,
        el_lte=args.el_rank_lte,
        ic50_max_nm=args.ic50_max_nm,
    )
    if args.sb_mode == "ic50_only":
        if iedb_ic50_col:
            ic50_src = f"{iedb_ic50_col} <{args.ic50_max_nm:g} nM (IEDB)"
        else:
            ic50_src = f"IC50<{args.ic50_max_nm:g} nM via local BA_score >{ba_min:.4f}"
        sb_legend = f"SB (IC50-only): {ic50_src}"
    else:
        if not spec_leg.use_ic50:
            ic50_src = "IC50/OFF"
        elif iedb_ic50_col:
            ic50_src = f"{iedb_ic50_col} <{args.ic50_max_nm:g} nM (IEDB)"
        else:
            ic50_src = f"IC50<{args.ic50_max_nm:g} nM via local BA_score >{ba_min:.4f}"
        sb_legend = (
            f"SB: immuno>{args.imm_min:g}, processing>{args.proc_min:g}, "
            f"EL rank {el_cmp}{args.el_rank_max:g}%, {ic50_src}"
        )

    sig_sb = (
        sig.loc[
            sb_mask_fig5_defaults(
                sig,
                el_max=args.el_rank_max,
                el_lte=args.el_rank_lte,
                ba_min=ba_min,
                ic50_max_nm=args.ic50_max_nm,
                iedb_ic50_col=iedb_ic50_col,
                imm_min=args.imm_min,
                proc_min=args.proc_min,
                sb_mode=args.sb_mode,
            )
        ].copy()
        if sig is not None
        else None
    )
    cod_sb = cod.loc[
        sb_mask_fig5_defaults(
            cod,
            el_max=args.el_rank_max,
            el_lte=args.el_rank_lte,
            ba_min=ba_min,
            ic50_max_nm=args.ic50_max_nm,
            iedb_ic50_col=iedb_ic50_col,
            imm_min=args.imm_min,
            proc_min=args.proc_min,
            sb_mode=args.sb_mode,
        )
    ].copy()

    out_dir = args.out_dir
    stem = args.output_stem

    if args.panels == "abc":
        assert sig_sb is not None and freqs is not None
        tab_a = fig5a_table(sig_sb, freqs)
        tab_a.to_csv(out_dir / f"{stem}_5a_epitopes_vs_allele_frequency.csv", index=False)
        plot_fig5a(
            tab_a,
            out_dir / f"{stem}_5a_epitopes_vs_allele_frequency.png",
            "Significant lncRNA micropeptides (translated)",
            sb_legend,
            y_metric=args.fig5a_y_metric,
        )

        cnt_b = fig5bc_sharing(sig_sb)
        h5b = sharing_histogram_full(cnt_b, args.n_panel_alleles, total_sb_row_instances=len(sig_sb))
        h5b.to_csv(out_dir / f"{stem}_5b_epitope_sharing_across_alleles.csv", index=False)
        plot_fig5bc_legacy(
            h5b,
            out_dir / f"{stem}_5b_epitope_sharing.png",
            "Epitope sharing across alleles (Tr-lncRNA-MPs)",
            sb_legend,
            bar_color=pal.SIG_LNC,
            y_metric=args.sharing_y_metric,
        )

    cnt_c = fig5bc_sharing(cod_sb)
    h5c = sharing_histogram_full(cnt_c, args.n_panel_alleles, total_sb_row_instances=len(cod_sb))
    h5c.to_csv(out_dir / f"{stem}_5c_epitope_sharing_across_alleles.csv", index=False)
    plot_fig5bc_legacy(
        h5c,
        out_dir / f"{stem}_5c_epitope_sharing.png",
        "Epitope sharing across alleles (coding control)",
        sb_legend,
        bar_color=pal.SIG_LNC,
        y_metric=args.sharing_y_metric,
    )

    if not args.no_repo_mirror:
        mirror_paths = [out_dir / f"{stem}_5c_epitope_sharing.png"]
        if args.panels == "abc":
            mirror_paths = [
                out_dir / f"{stem}_5a_epitopes_vs_allele_frequency.png",
                out_dir / f"{stem}_5b_epitope_sharing.png",
                out_dir / f"{stem}_5c_epitope_sharing.png",
            ]
        _mirror(*mirror_paths)

    print(f"Wrote under {out_dir} prefix {stem}_* (panels={args.panels})")
    print(f"[stats] IC50 filter: {iedb_ic50_col or 'local BA_score (no IEDB IC50 column)'}")
    if sig_sb is not None:
        print(f"[stats] sig SB epitopes × alleles={len(sig_sb)} unique 9-mers={sig_sb['Peptide'].nunique()}")
    print(f"[stats] coding SB epitopes × alleles={len(cod_sb)} unique 9-mers={cod_sb['Peptide'].nunique()}")
    if not iedb_ic50_col:
        print(f"[stats] BA_score threshold (IC50<{args.ic50_max_nm:g} nM): > {ba_min:.6f}")


if __name__ == "__main__":
    main()
