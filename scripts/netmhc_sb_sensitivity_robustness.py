#!/usr/bin/env python3
"""
Sensitivity (one-dimensional threshold sweeps) and robustness (leave-one-filter-out)
for the combined SB definition on merged NetMHC + IEDB TSVs.

Baseline defaults match ``plot_fig5abc_netmhc_sb_triple.py`` / ``netmhc_sb_core.FIG5_IEDB_*`` and
``--sb-mode`` (``full`` / ``no_ic50`` / ``ic50_only``). Plots default to **SB epitopes × alleles**
(merged-table row counts); use ``--plot-metric unique`` for unique-9-mer curves.

Writes under ``data/netmhc/figures/`` by default:

  - ``sb_threshold_sensitivity_robustness.csv`` — sweeps + LOO + baseline row.
  - ``sb_threshold_sensitivity_robustness_fold_change_vs_baseline.csv`` — %% and fold vs baseline.
  - ``sb_threshold_sensitivity_robustness.png`` — multi-panel figure.

This is **catalog Figure 5** (cohort IEDB+NetMHC SB), not Figure 6 (TTN-AS1). See ``docs/figure_catalog.md``.
"""

from __future__ import annotations

import argparse
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
from matplotlib.axes import Axes
import numpy as np
import pandas as pd

from netmhc_sb_core import (
    FIG5_IEDB_EL_RANK_MAX_DEFAULT,
    FIG5_IEDB_IC50_MAX_NM_DEFAULT,
    FIG5_IEDB_IMM_MIN_DEFAULT,
    FIG5_IEDB_PROC_MIN_DEFAULT,
    SBSpec,
    ba_score_min_for_ic50_lt,
    pick_iedb_ic50_column,
    sb_mask_spec,
    sb_spec_from_mode,
)


def cohort_stats(df: pd.DataFrame, m: pd.Series) -> dict[str, int]:
    sub = df.loc[m]
    return {
        "n_rows": int(len(sub)),
        "n_unique_peptides": int(sub["Peptide"].nunique()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sig-tsv", type=Path, default=Path("data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv"))
    ap.add_argument(
        "--coding-tsv",
        type=Path,
        default=Path("data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv"),
    )
    ap.add_argument("--out-dir", type=Path, default=Path("data/netmhc/figures"))
    ap.add_argument("--out-stem", type=str, default="sb_threshold_sensitivity_robustness")
    ap.add_argument("--el-rank-lte", action="store_true")
    ap.add_argument("--baseline-imm", type=float, default=FIG5_IEDB_IMM_MIN_DEFAULT)
    ap.add_argument("--baseline-proc", type=float, default=FIG5_IEDB_PROC_MIN_DEFAULT)
    ap.add_argument("--baseline-el-max", type=float, default=FIG5_IEDB_EL_RANK_MAX_DEFAULT)
    ap.add_argument("--baseline-ic50-nm", type=float, default=FIG5_IEDB_IC50_MAX_NM_DEFAULT)
    ap.add_argument(
        "--plot-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Y-axis / fold-change: SB epitopes × alleles vs unique 9-mers (default instances).",
    )
    ap.add_argument(
        "--sb-mode",
        choices=("full", "no_ic50", "ic50_only"),
        default="full",
        help="Baseline SB composition (see netmhc_sb_core.sb_spec_from_mode). Sweeps vary one dimension at a time.",
    )
    args = ap.parse_args()

    y_sig_col = "sig_n_rows" if args.plot_metric == "instances" else "sig_n_unique_peptides"
    y_cod_col = "coding_n_rows" if args.plot_metric == "instances" else "coding_n_unique_peptides"
    y_label = "SB epitopes × alleles" if args.plot_metric == "instances" else "Unique 9-mers"
    ratio_col = (
        "ratio_sig_over_coding_rows" if args.plot_metric == "instances" else "ratio_sig_over_coding_unique"
    )

    hdr = set(pd.read_csv(args.sig_tsv, sep="\t", nrows=0).columns) | set(
        pd.read_csv(args.coding_tsv, sep="\t", nrows=0).columns
    )
    iedb_ic50_col = pick_iedb_ic50_column(sorted(hdr))
    usecols = {"allele", "Peptide", "EL_rank", "iedb_score", "iedb_processing_score"}
    if iedb_ic50_col:
        usecols.add(iedb_ic50_col)
    else:
        usecols.add("BA_score")

    base = sb_spec_from_mode(
        args.sb_mode,
        imm_min=args.baseline_imm,
        proc_min=args.baseline_proc,
        el_max=args.baseline_el_max,
        el_lte=args.el_rank_lte,
        ic50_max_nm=args.baseline_ic50_nm,
    )

    sig = pd.read_csv(args.sig_tsv, sep="\t", usecols=lambda c: c in usecols)
    cod = pd.read_csv(args.coding_tsv, sep="\t", usecols=lambda c: c in usecols)

    rows_long: list[dict[str, object]] = []

    def add_row(
        *,
        analysis: str,
        dimension: str,
        value: float | None,
        label: str,
        df_sig: pd.DataFrame,
        df_cod: pd.DataFrame,
        spec: SBSpec,
    ) -> None:
        bm = ba_score_min_for_ic50_lt(spec.ic50_max_nm)
        ms = sb_mask_spec(df_sig, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        mc = sb_mask_spec(df_cod, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        ss, cs = cohort_stats(df_sig, ms), cohort_stats(df_cod, mc)
        ratio_u = (
            float(ss["n_unique_peptides"]) / float(cs["n_unique_peptides"])
            if cs["n_unique_peptides"] > 0
            else float("nan")
        )
        ratio_r = float(ss["n_rows"]) / float(cs["n_rows"]) if cs["n_rows"] > 0 else float("nan")
        rows_long.append(
            {
                "analysis": analysis,
                "dimension": dimension,
                "threshold_value": value,
                "scenario_label": label,
                "sig_n_rows": ss["n_rows"],
                "sig_n_unique_peptides": ss["n_unique_peptides"],
                "coding_n_rows": cs["n_rows"],
                "coding_n_unique_peptides": cs["n_unique_peptides"],
                "ratio_sig_over_coding_unique": ratio_u,
                "ratio_sig_over_coding_rows": ratio_r,
                "ic50_source": iedb_ic50_col or "local_BA_score",
            }
        )

    # Baseline
    add_row(
        analysis="baseline",
        dimension="all",
        value=None,
        label="full_SB",
        df_sig=sig,
        df_cod=cod,
        spec=base,
    )

    # Leave-one-out (relax exactly one gate)
    loo_specs: list[tuple[str, SBSpec]] = [
        ("relax_immuno", replace(base, use_imm=False)),
        ("relax_processing", replace(base, use_proc=False)),
        ("relax_EL", replace(base, use_el=False)),
        ("relax_IC50_binding", replace(base, use_ic50=False)),
    ]
    for lab, sp in loo_specs:
        add_row(
            analysis="leave_one_out",
            dimension=lab,
            value=None,
            label=lab,
            df_sig=sig,
            df_cod=cod,
            spec=sp,
        )

    # Sensitivity sweeps (one dimension at a time)
    imm_sweep = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3]
    proc_sweep = [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
    el_sweep = [1.0, 2.0, 3.0, 5.0, 10.0, 20.0]
    ic50_sweep = [50.0, 100.0, 150.0, 250.0, 500.0, 1000.0, 5000.0]

    for v in imm_sweep:
        sp = replace(base, imm_min=v)
        add_row(
            analysis="sensitivity",
            dimension="immuno_min",
            value=v,
            label=f"imm>{v}",
            df_sig=sig,
            df_cod=cod,
            spec=sp,
        )

    for v in proc_sweep:
        sp = replace(base, proc_min=v)
        add_row(
            analysis="sensitivity",
            dimension="processing_min",
            value=v,
            label=f"proc>{v}",
            df_sig=sig,
            df_cod=cod,
            spec=sp,
        )

    for v in el_sweep:
        sp = replace(base, el_max=v)
        # BA min must match ic50 gate when sweeping EL only — keep baseline IC50
        add_row(
            analysis="sensitivity",
            dimension="el_rank_max_pct",
            value=v,
            label=f"EL{'≤' if base.el_lte else '<'}{v}",
            df_sig=sig,
            df_cod=cod,
            spec=sp,
        )

    for v in ic50_sweep:
        sp = replace(base, ic50_max_nm=v)
        add_row(
            analysis="sensitivity",
            dimension="ic50_max_nm",
            value=v,
            label=f"IC50<{v:g}nM",
            df_sig=sig,
            df_cod=cod,
            spec=sp,
        )

    out = pd.DataFrame(rows_long)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.out_stem
    csv_path = args.out_dir / f"{stem}.csv"
    out.to_csv(csv_path, index=False)

    b0 = out[out["analysis"] == "baseline"].iloc[0]
    b_sig = max(int(b0[y_sig_col]), 1)
    b_cod = max(int(b0[y_cod_col]), 1)
    b_ratio = float(b0[ratio_col])
    fc = out.copy()
    fc["sig_pct_of_baseline_pct"] = 100.0 * fc[y_sig_col].astype(float) / float(b_sig)
    fc["coding_pct_of_baseline_pct"] = 100.0 * fc[y_cod_col].astype(float) / float(b_cod)
    fc["sig_fold_vs_baseline"] = fc[y_sig_col].astype(float) / float(b_sig)
    fc["coding_fold_vs_baseline"] = fc[y_cod_col].astype(float) / float(b_cod)
    fc["ratio_pct_of_baseline_pct"] = np.where(
        np.isfinite(b_ratio) & (b_ratio != 0.0),
        100.0 * fc[ratio_col].astype(float) / b_ratio,
        np.nan,
    )
    fc_path = args.out_dir / f"{stem}_fold_change_vs_baseline.csv"
    fc.to_csv(fc_path, index=False)

    loo = out[out["analysis"] == "leave_one_out"].copy()
    base_row = out[out["analysis"] == "baseline"].iloc[0]
    x_labels = ["baseline\n(full SB)"] + [s.replace("relax_", "").replace("_", "\n") for s in loo["dimension"]]
    x = np.arange(len(x_labels))
    sig_u = [int(base_row[y_sig_col])] + loo[y_sig_col].astype(int).tolist()
    cod_u = [int(base_row[y_cod_col])] + loo[y_cod_col].astype(int).tolist()
    w = 0.35

    def line_panel(dim: str, title: str, ax: Axes) -> None:
        sub = out[(out["analysis"] == "sensitivity") & (out["dimension"] == dim)].sort_values("threshold_value")
        ax.plot(sub["threshold_value"], sub[y_sig_col], "o-", color=pal.SIG_LNC, label=f"Sig ({args.plot_metric})")
        ax.plot(sub["threshold_value"], sub[y_cod_col], "s-", color=pal.CODING_CONTROL, label=f"Coding ({args.plot_metric})")
        ax.set_xlabel(title)
        ax.set_ylabel(y_label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0], hspace=0.4, wspace=0.32)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.bar(x - w / 2, sig_u, width=w, label=f"Sig lnc ({args.plot_metric})", color=pal.SIG_LNC, edgecolor="black", linewidth=0.35)
    ax0.bar(x + w / 2, cod_u, width=w, label=f"Coding control ({args.plot_metric})", color=pal.CODING_CONTROL, edgecolor="black", linewidth=0.35)
    ax0.set_xticks(x)
    ax0.set_xticklabels(x_labels, fontsize=9)
    ax0.set_ylabel(f"{y_label} passing filters")
    ax0.set_title("Robustness: leave-one-filter-out (others at baseline)")
    ax0.legend(loc="upper right", fontsize=9)
    ax0.grid(axis="y", alpha=0.3)

    ax1 = fig.add_subplot(gs[1, 0])
    line_panel("immuno_min", "Immunogenicity cutoff (iedb_score >)", ax1)
    ax2 = fig.add_subplot(gs[1, 1])
    line_panel("processing_min", "Processing cutoff (iedb_processing_score >)", ax2)
    ax3 = fig.add_subplot(gs[2, 0])
    line_panel("el_rank_max_pct", "EL %rank cutoff (strict < unless --el-rank-lte)", ax3)
    ax4 = fig.add_subplot(gs[2, 1])
    sub_ic = out[(out["analysis"] == "sensitivity") & (out["dimension"] == "ic50_max_nm")].sort_values("threshold_value")
    ax4.plot(sub_ic["threshold_value"], sub_ic[y_sig_col], "o-", color=pal.SIG_LNC, label=f"Sig ({args.plot_metric})")
    ax4.plot(sub_ic["threshold_value"], sub_ic[y_cod_col], "s-", color=pal.CODING_CONTROL, label=f"Coding ({args.plot_metric})")
    ax4.set_xscale("log")
    ax4.set_xlabel("IC50 max (nM), strict < cutoff")
    ax4.set_ylabel(y_label)
    ax4.grid(alpha=0.3)
    ax4.legend(fontsize=8)
    fig.suptitle("SB sensitivity and robustness", fontsize=12, y=0.995)
    png_path = args.out_dir / f"{stem}.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {csv_path}")
    print(f"Wrote {fc_path}")
    print(f"Wrote {png_path}")
    print(f"[info] IC50 source: {iedb_ic50_col or 'local BA_score'}")
    print("[info] baseline:", base)
    print(f"\n--- Leave-one-out ({args.plot_metric}) ---")
    print(
        out[out["analysis"] == "leave_one_out"][["scenario_label", y_sig_col, y_cod_col]].to_string(
            index=False
        )
    )
    b = out[out["analysis"] == "baseline"].iloc[0]
    bs, bc = int(b[y_sig_col]), int(b[y_cod_col])
    print(f"\n--- Baseline ({args.plot_metric}) ---")
    print(
        f"sig {y_sig_col}={bs}, coding {y_cod_col}={bc}, "
        f"ratio_unique={b['ratio_sig_over_coding_unique']:.3f}, ratio_rows={b['ratio_sig_over_coding_rows']:.3f}"
    )
    print(f"\n--- Delta {args.plot_metric} vs baseline (sig / coding) ---")
    for _, r in out[out["analysis"] == "leave_one_out"].iterrows():
        ds = int(r[y_sig_col]) - bs
        dc = int(r[y_cod_col]) - bc
        print(f"  {r['scenario_label']}: {ds:+d} / {dc:+d}")
    el_row = out[(out["analysis"] == "leave_one_out") & (out["dimension"] == "relax_EL")]
    if len(el_row) and int(el_row.iloc[0][y_sig_col]) == bs:
        print(
            "\n[note] Relaxing EL_rank did not change counts: with your other gates, "
            "no merged row passes immuno+processing+IC50 while having EL_rank at or above the baseline cut. "
            "The EL filter is redundant in this table (still interpretable as a conservative statement).",
            flush=True,
        )


if __name__ == "__main__":
    main()
