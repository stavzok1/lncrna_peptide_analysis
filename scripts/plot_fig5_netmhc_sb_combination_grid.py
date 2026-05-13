#!/usr/bin/env python3
"""
**Catalog Figure 5 — supplement (merged NetMHC + IEDB):** SB threshold **combinations**.

This is **not** manuscript Figure 6. Manuscript **Figure 6** is TTN-AS1 allele coverage
(``plot_figure6_ttn_as1_allele_coverage.py``); see ``docs/figure_catalog.md``.

This script builds a Cartesian grid over immunogenicity / processing / EL / IC50 gates
on ``*_with_iedb.tsv`` tables, fold-change vs baseline, curated 3x3 sharing histograms,
and a heatmap slice. It complements:

  - ``plot_fig5abc_netmhc_sb_triple.py`` — single chosen SB definition (5A–5C style).
  - ``scripts/netmhc_sb_sensitivity_robustness.py`` — one-at-a-time sweeps + LOO + fold-change.

Default output folder: ``data/netmhc/figures/fig5_netmhc_sb_combinations/``.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import replace
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = _ROOT / "scripts"
for _p in (str(_ROOT), str(_SCRIPTS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import figure_palettes as pal  # noqa: E402

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from netmhc_sb_core import (
    FIG5_IEDB_EL_RANK_MAX_DEFAULT,
    FIG5_IEDB_IC50_MAX_NM_DEFAULT,
    FIG5_IEDB_IMM_MIN_DEFAULT,
    FIG5_IEDB_PROC_MIN_DEFAULT,
    SBSpec,
    ba_score_min_for_ic50_lt,
    pick_iedb_ic50_column,
    sb_mask_spec,
    spec_label,
    spec_profile_id,
)


def parse_float_list(s: str) -> list[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def fig5bc_sharing(sb: pd.DataFrame) -> pd.Series:
    return sb.groupby("Peptide")["allele"].nunique()


def sharing_histogram_full(
    cnt: pd.Series,
    max_alleles: int,
    *,
    total_sb_row_instances: int | None = None,
) -> pd.DataFrame:
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


def cohort_stats(df: pd.DataFrame, m: pd.Series) -> dict[str, int]:
    sub = df.loc[m]
    return {
        "n_rows": int(len(sub)),
        "n_unique_peptides": int(sub["Peptide"].nunique()),
    }


def plot_sharing_mini(
    hist: pd.DataFrame,
    ax,
    *,
    bar_color: str,
    ymax: float | None,
    y_column: str,
) -> None:
    x = hist["n_alleles_per_epitope"].to_numpy()
    h = hist[y_column].to_numpy()
    ax.bar(x, h, width=0.82, color=bar_color, edgecolor="black", linewidth=0.25, align="center")
    ax.set_xticks(x[:: max(1, len(x) // 13)])
    ax.set_xlabel("N alleles", fontsize=7)
    ax.set_ylabel("SB epitopes per allele" if y_column == "n_sb_rows_in_bin" else "Unique 9-mers in bin", fontsize=7)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=2))
    ax.tick_params(axis="both", labelsize=7)
    ax.grid(axis="y", alpha=0.25)
    if ymax is not None:
        ax.set_ylim(0, ymax * 1.05)
    if "total_epitopes" in hist.columns and len(hist):
        te = int(hist["total_epitopes"].iloc[0])
        tr = int(hist["total_sb_row_instances"].iloc[0]) if "total_sb_row_instances" in hist.columns else None
        lab = f"uniq={te}" + (f", rows={tr}" if tr is not None else "")
        ax.text(
            0.97,
            0.93,
            lab,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.5"},
        )


def default_panel_specs(baseline: SBSpec) -> list[tuple[str, SBSpec]]:
    b = baseline
    return [
        ("01_baseline", b),
        ("02_ic50_500nM", replace(b, ic50_max_nm=500.0)),
        ("03_ic50_1000nM", replace(b, ic50_max_nm=1000.0)),
        ("04_immuno_gt0p05", replace(b, imm_min=0.05)),
        ("05_proc_gt1p25", replace(b, proc_min=1.25)),
        ("06_ic50_500_imm0p05", replace(b, ic50_max_nm=500.0, imm_min=0.05)),
        ("07_ic500_proc125", replace(b, ic50_max_nm=500.0, proc_min=1.25)),
        ("08_moderate", replace(b, ic50_max_nm=500.0, imm_min=0.05, proc_min=1.25, el_max=5.0)),
        (
            "09_permissive",
            SBSpec(imm_min=0.0, proc_min=1.0, el_max=10.0, el_lte=b.el_lte, ic50_max_nm=1000.0),
        ),
    ]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sig-tsv", type=Path, default=Path("data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv"))
    ap.add_argument(
        "--coding-tsv",
        type=Path,
        default=Path("data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv"),
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/netmhc/figures/fig5_netmhc_sb_combinations"),
        help="Folder for combination-grid products (Figure 5 supplement, merged IEDB+NetMHC).",
    )
    ap.add_argument("--el-rank-lte", action="store_true")
    ap.add_argument("--baseline-imm", type=float, default=FIG5_IEDB_IMM_MIN_DEFAULT)
    ap.add_argument("--baseline-proc", type=float, default=FIG5_IEDB_PROC_MIN_DEFAULT)
    ap.add_argument("--baseline-el-max", type=float, default=FIG5_IEDB_EL_RANK_MAX_DEFAULT)
    ap.add_argument("--baseline-ic50-nm", type=float, default=FIG5_IEDB_IC50_MAX_NM_DEFAULT)
    ap.add_argument("--n-panel-alleles", type=int, default=27)
    ap.add_argument(
        "--sharing-y-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Mini sharing grids: y-axis SB epitopes per allele (instances) vs unique 9-mers in bin (default instances).",
    )
    ap.add_argument(
        "--imm-grid",
        type=str,
        default="0.05,0.1,0.15",
        help="Comma-separated iedb_score lower bounds (imm > value).",
    )
    ap.add_argument("--proc-grid", type=str, default="1.25,1.5,1.75")
    ap.add_argument(
        "--el-grid",
        type=str,
        default="1,2,5,10",
        help="EL %%rank upper bound ( < or <= ); includes 1 to align with default baseline EL gate.",
    )
    ap.add_argument("--ic50-grid", type=str, default="150,500,1000", help="IC50 strict upper bound (nM).")
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    hdr = set(pd.read_csv(args.sig_tsv, sep="\t", nrows=0).columns) | set(
        pd.read_csv(args.coding_tsv, sep="\t", nrows=0).columns
    )
    iedb_ic50_col = pick_iedb_ic50_column(sorted(hdr))
    usecols = {"allele", "Peptide", "EL_rank", "iedb_score", "iedb_processing_score"}
    if iedb_ic50_col:
        usecols.add(iedb_ic50_col)
    else:
        usecols.add("BA_score")

    baseline = SBSpec(
        imm_min=args.baseline_imm,
        proc_min=args.baseline_proc,
        el_max=args.baseline_el_max,
        el_lte=args.el_rank_lte,
        ic50_max_nm=args.baseline_ic50_nm,
    )

    imm_vals = parse_float_list(args.imm_grid)
    proc_vals = parse_float_list(args.proc_grid)
    el_vals = parse_float_list(args.el_grid)
    ic50_vals = parse_float_list(args.ic50_grid)
    y_col = "n_sb_rows_in_bin" if args.sharing_y_metric == "instances" else "n_epitopes"

    sig = pd.read_csv(args.sig_tsv, sep="\t", usecols=lambda c: c in usecols)
    cod = pd.read_csv(args.coding_tsv, sep="\t", usecols=lambda c: c in usecols)

    grid_rows: list[dict[str, object]] = []
    for imm, proc, el_m, ic50 in itertools.product(imm_vals, proc_vals, el_vals, ic50_vals):
        spec = SBSpec(
            imm_min=imm,
            proc_min=proc,
            el_max=el_m,
            el_lte=args.el_rank_lte,
            ic50_max_nm=ic50,
        )
        bm = ba_score_min_for_ic50_lt(spec.ic50_max_nm)
        ms = sb_mask_spec(sig, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        mc = sb_mask_spec(cod, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        ss, cs = cohort_stats(sig, ms), cohort_stats(cod, mc)
        ratio = (
            float(ss["n_unique_peptides"]) / float(cs["n_unique_peptides"])
            if cs["n_unique_peptides"] > 0
            else float("nan")
        )
        pid = spec_profile_id(spec)
        grid_rows.append(
            {
                "profile_id": pid,
                "imm_min": imm,
                "proc_min": proc,
                "el_max_pct": el_m,
                "el_lte": args.el_rank_lte,
                "ic50_max_nm": ic50,
                "filter_label": spec_label(spec),
                "sig_n_rows": ss["n_rows"],
                "sig_n_unique_peptides": ss["n_unique_peptides"],
                "coding_n_rows": cs["n_rows"],
                "coding_n_unique_peptides": cs["n_unique_peptides"],
                "ratio_sig_over_coding_unique": ratio,
                "ic50_source": iedb_ic50_col or "local_BA_score",
            }
        )

    grid_df = pd.DataFrame(grid_rows)
    grid_csv = out_dir / "counts_combination_grid.csv"
    grid_df.to_csv(grid_csv, index=False)

    mask_base = (
        (grid_df["imm_min"] == baseline.imm_min)
        & (grid_df["proc_min"] == baseline.proc_min)
        & (grid_df["el_max_pct"] == baseline.el_max)
        & (grid_df["ic50_max_nm"] == baseline.ic50_max_nm)
    )
    if not mask_base.any():
        b_bm = ba_score_min_for_ic50_lt(baseline.ic50_max_nm)
        bms = sb_mask_spec(sig, baseline, ba_min=b_bm, iedb_ic50_col=iedb_ic50_col)
        bmc = sb_mask_spec(cod, baseline, ba_min=b_bm, iedb_ic50_col=iedb_ic50_col)
        b_sig = max(cohort_stats(sig, bms)["n_unique_peptides"], 1)
        b_cod = max(cohort_stats(cod, bmc)["n_unique_peptides"], 1)
        b_ratio = (float(b_sig) / float(b_cod)) if b_cod > 0 else float("nan")
    else:
        br = grid_df.loc[mask_base].iloc[0]
        b_sig = max(int(br["sig_n_unique_peptides"]), 1)
        b_cod = max(int(br["coding_n_unique_peptides"]), 1)
        b_ratio = float(br["ratio_sig_over_coding_unique"])

    fc = grid_df.copy()
    fc["sig_pct_of_baseline_pct"] = 100.0 * fc["sig_n_unique_peptides"].astype(float) / float(b_sig)
    fc["coding_pct_of_baseline_pct"] = 100.0 * fc["coding_n_unique_peptides"].astype(float) / float(b_cod)
    fc["sig_fold_vs_baseline"] = fc["sig_n_unique_peptides"].astype(float) / float(b_sig)
    fc["coding_fold_vs_baseline"] = fc["coding_n_unique_peptides"].astype(float) / float(b_cod)
    fc["ratio_pct_of_baseline_pct"] = np.where(
        np.isfinite(b_ratio) & (b_ratio != 0.0),
        100.0 * fc["ratio_sig_over_coding_unique"].astype(float) / b_ratio,
        np.nan,
    )
    fc_path = out_dir / "fold_change_vs_baseline_grid.csv"
    fc.to_csv(fc_path, index=False)

    panels = default_panel_specs(baseline)
    ncols = 3
    nrows = int(np.ceil(len(panels) / ncols))

    sig_hmax = 0.0
    cod_hmax = 0.0
    panel_hist_cache: list[tuple[str, SBSpec, pd.DataFrame, pd.DataFrame, int, int]] = []
    for short_name, spec in panels:
        bm = ba_score_min_for_ic50_lt(spec.ic50_max_nm)
        ms = sb_mask_spec(sig, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        mc = sb_mask_spec(cod, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        sig_sb = sig.loc[ms]
        cod_sb = cod.loc[mc]
        cnt_s = fig5bc_sharing(sig_sb)
        cnt_c = fig5bc_sharing(cod_sb)
        h_s = sharing_histogram_full(cnt_s, args.n_panel_alleles, total_sb_row_instances=len(sig_sb))
        h_c = sharing_histogram_full(cnt_c, args.n_panel_alleles, total_sb_row_instances=len(cod_sb))
        pid = spec_profile_id(spec)
        h_s.to_csv(out_dir / f"sharing_{pid}_sig.csv", index=False)
        h_c.to_csv(out_dir / f"sharing_{pid}_cod.csv", index=False)
        sig_hmax = max(sig_hmax, float(h_s[y_col].max()), 1.0)
        cod_hmax = max(cod_hmax, float(h_c[y_col].max()), 1.0)
        su, cu = cohort_stats(sig, ms)["n_unique_peptides"], cohort_stats(cod, mc)["n_unique_peptides"]
        panel_hist_cache.append((short_name, spec, h_s, h_c, su, cu))

    def draw_grid(which: str, ymax: float) -> None:
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 3.2 * nrows), squeeze=False)
        color = pal.SIG_LNC if which == "sig" else pal.CODING_CONTROL
        for i, (short_name, spec, h_s, h_c, su, cu) in enumerate(panel_hist_cache):
            r, c = divmod(i, ncols)
            ax = axes[r][c]
            hist = h_s if which == "sig" else h_c
            nuniq = su if which == "sig" else cu
            plot_sharing_mini(hist, ax, bar_color=color, ymax=ymax, y_column=y_col)
            lab = spec_label(spec)
            if len(lab) > 48:
                lab = lab[:45] + "..."
            ax.set_title(f"{short_name}\nuniq={nuniq}\n{lab}", fontsize=7)
        for j in range(len(panel_hist_cache), nrows * ncols):
            r, c = divmod(j, ncols)
            axes[r][c].set_visible(False)
        cohort_lab = "sig lnc" if which == "sig" else "coding control"
        fig.suptitle(
            f"Fig 5 supplement — SB combination grid, epitope sharing ({cohort_lab}) | "
            f"IC50: {iedb_ic50_col or 'local BA'}",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(out_dir / f"fig5_netmhc_sb_combo_{which}_sharing_grid.png", dpi=220, bbox_inches="tight")
        plt.close(fig)

    draw_grid("sig", sig_hmax)
    draw_grid("cod", cod_hmax)

    names = [t[0] for t in panel_hist_cache]
    sig_u = [t[4] for t in panel_hist_cache]
    cod_u = [t[5] for t in panel_hist_cache]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    y = np.arange(len(names))
    h = 0.35
    ax.barh(y - h / 2, sig_u, height=h, label="Sig unique 9-mers", color=pal.SIG_LNC, edgecolor="black", linewidth=0.3)
    ax.barh(y + h / 2, cod_u, height=h, label="Coding unique 9-mers", color=pal.CODING_CONTROL, edgecolor="black", linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Unique peptides passing filters")
    ax.set_title("Fig 5 supplement — SB panel profiles: unique 9-mer counts")
    ax.legend()
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "fig5_netmhc_sb_combo_overview_counts.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    proc_fix = min(proc_vals)
    el_fix = min(el_vals)
    sub = grid_df[(grid_df["proc_min"] == proc_fix) & (grid_df["el_max_pct"] == el_fix)]
    pivot = None
    if len(sub) == len(imm_vals) * len(ic50_vals) and len(imm_vals) > 0 and len(ic50_vals) > 0:
        pivot = sub.pivot_table(
            index="imm_min",
            columns="ic50_max_nm",
            values="sig_n_unique_peptides",
            aggfunc="first",
        )
        if pivot.size == 0 or pivot.shape[0] < 1 or pivot.shape[1] < 1:
            pivot = None
    if pivot is not None:
        fig, ax = plt.subplots(figsize=(7.5, 5))
        im = ax.imshow(pivot.values, aspect="auto", cmap=pal.sequential_count_heatmap(), origin="lower")
        ax.set_xticks(np.arange(pivot.shape[1]))
        ax.set_xticklabels([str(x) for x in pivot.columns], rotation=45, ha="right")
        ax.set_yticks(np.arange(pivot.shape[0]))
        ax.set_yticklabels([str(x) for x in pivot.index])
        ax.set_xlabel("IC50 max (nM)")
        ax.set_ylabel("Immunogenicity min (iedb_score >)")
        ax.set_title(
            f"Fig 5 supplement — sig unique 9-mers: immuno x IC50\n(proc={proc_fix}, EL<{el_fix}%)"
        )
        fig.colorbar(im, ax=ax, label="Sig unique 9-mers")
        fig.tight_layout()
        fig.savefig(out_dir / "fig5_netmhc_sb_combo_heatmap_sig_imm_ic50.png", dpi=220, bbox_inches="tight")
        plt.close(fig)

    print(f"Wrote {grid_csv} ({len(grid_df)} combinations)")
    print(f"Wrote {fc_path}")
    print(f"Wrote fig5_netmhc_sb_combo_* under {out_dir}")
    print(f"[baseline] unique sig={b_sig}, coding={b_cod}, ratio~{b_ratio:.3f}")


if __name__ == "__main__":
    main()
