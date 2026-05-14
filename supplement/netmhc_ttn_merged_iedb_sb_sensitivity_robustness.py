#!/usr/bin/env python3
"""
**Figure 6 supplement (TTN-AS1):** one-dimensional SB threshold sweeps + leave-one-filter-out on the
**merged** NetMHC wide XLS + IEDB CSV (same join as ``--gating iedb_sb`` in
``plot_figure6_ttn_as1_allele_coverage.py``).

This mirrors the **design** of ``netmhc_sb_sensitivity_robustness.py`` (baseline → LOO relax one
gate at a time → 1D sweeps per axis), but for the **single-peptide** TTN cohort (counts only, no
sig/coding ratio).

**Contrast:** ``plot_fig6_ttn_merged_iedb_sb_combination_grid.py`` evaluates the **full Cartesian
product** of grid cutoffs; this script varies **one dimension at a time** (plus LOO).

Default outputs under ``data/netmhc/figures/fig6_ttn_merged_iedb_1d_sensitivity/`` unless
``--out-dir`` is set. ``generate_netmhc_fig5_fig6_supplement.py`` writes to
``figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_merged_iedb_1d_sensitivity_loo/``.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent.parent / "manuscript"
for _p in (str(_REPO), str(_REPO / "scripts"), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import argparse
from dataclasses import replace

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

import plot_figure6_ttn_as1_allele_coverage as ttn6  # noqa: E402


def cohort_row_stats(df: pd.DataFrame, m: pd.Series) -> dict[str, int]:
    sub = df.loc[m]
    return {
        "n_rows": int(len(sub)),
        "n_unique_peptides": int(sub["Peptide"].nunique()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--netmhc-xls",
        type=Path,
        default=Path("data/netmhc/netmhcpan_ttn_as1_108065.xls"),
    )
    ap.add_argument("--iedb-csv", type=Path, default=ttn6.TTN_IEDB_SYNTHETIC_CSV)
    ap.add_argument("--iedb-parent-input-seq-id", type=str, default=ttn6.TTN_IEDB_PARENT_SEQ_ID_DEFAULT)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/netmhc/figures/fig6_ttn_merged_iedb_1d_sensitivity"),
    )
    ap.add_argument("--out-stem", type=str, default="ttn_merged_iedb_sb_sensitivity_robustness")
    ap.add_argument("--el-rank-lte", action="store_true")
    ap.add_argument("--baseline-imm", type=float, default=FIG5_IEDB_IMM_MIN_DEFAULT)
    ap.add_argument("--baseline-proc", type=float, default=FIG5_IEDB_PROC_MIN_DEFAULT)
    ap.add_argument("--baseline-el-max", type=float, default=FIG5_IEDB_EL_RANK_MAX_DEFAULT)
    ap.add_argument("--baseline-ic50-nm", type=float, default=FIG5_IEDB_IC50_MAX_NM_DEFAULT)
    ap.add_argument(
        "--plot-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Primary y-axis for sweep panels: SB rows (instances) vs unique 9-mers.",
    )
    ap.add_argument(
        "--sb-mode",
        choices=("full", "no_ic50", "ic50_only"),
        default="full",
        help="Baseline SB composition (see netmhc_sb_core.sb_spec_from_mode).",
    )
    args = ap.parse_args()

    if not args.netmhc_xls.is_file():
        raise SystemExit(f"Missing NetMHC XLS: {args.netmhc_xls}")
    if not args.iedb_csv.is_file():
        raise SystemExit(f"Missing IEDB CSV: {args.iedb_csv}")

    y_col = "n_rows" if args.plot_metric == "instances" else "n_unique_peptides"
    y_label = "SB epitope×allele rows" if args.plot_metric == "instances" else "Unique 9-mers"

    alleles, starts, peps, ba, ba_rank, el_rank = ttn6.parse_wide_netmhc_xls_rows(args.netmhc_xls)
    long_df = ttn6.build_ttn_long_for_iedb_merge(
        starts, peps, alleles, ba, ba_rank, el_rank, str(args.iedb_parent_input_seq_id).strip()
    )
    merged = ttn6.merge_iedb_csv(long_df, args.iedb_csv)

    hdr = sorted(set(merged.columns))
    iedb_ic50_col = pick_iedb_ic50_column(hdr)

    base = sb_spec_from_mode(
        args.sb_mode,
        imm_min=args.baseline_imm,
        proc_min=args.baseline_proc,
        el_max=args.baseline_el_max,
        el_lte=args.el_rank_lte,
        ic50_max_nm=args.baseline_ic50_nm,
    )

    rows_long: list[dict[str, object]] = []

    def add_row(
        *,
        analysis: str,
        dimension: str,
        value: float | None,
        label: str,
        spec: SBSpec,
    ) -> None:
        bm = ba_score_min_for_ic50_lt(spec.ic50_max_nm)
        ms = sb_mask_spec(merged, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        st = cohort_row_stats(merged, ms)
        rows_long.append(
            {
                "analysis": analysis,
                "dimension": dimension,
                "threshold_value": value,
                "scenario_label": label,
                "n_rows": st["n_rows"],
                "n_unique_peptides": st["n_unique_peptides"],
                "ic50_source": iedb_ic50_col or "local_BA_score",
            }
        )

    add_row(analysis="baseline", dimension="all", value=None, label="full_SB", spec=base)

    loo_specs: list[tuple[str, SBSpec]] = [
        ("relax_immuno", replace(base, use_imm=False)),
        ("relax_processing", replace(base, use_proc=False)),
        ("relax_EL", replace(base, use_el=False)),
        ("relax_IC50_binding", replace(base, use_ic50=False)),
    ]
    for lab, sp in loo_specs:
        add_row(analysis="leave_one_out", dimension=lab, value=None, label=lab, spec=sp)

    imm_sweep = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3]
    proc_sweep = [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
    el_sweep = [1.0, 2.0, 3.0, 5.0, 10.0, 20.0]
    ic50_sweep = [50.0, 100.0, 150.0, 250.0, 500.0, 1000.0, 5000.0]

    for v in imm_sweep:
        add_row(
            analysis="sensitivity",
            dimension="immuno_min",
            value=v,
            label=f"imm>{v}",
            spec=replace(base, imm_min=v),
        )
    for v in proc_sweep:
        add_row(
            analysis="sensitivity",
            dimension="processing_min",
            value=v,
            label=f"proc>{v}",
            spec=replace(base, proc_min=v),
        )
    for v in el_sweep:
        add_row(
            analysis="sensitivity",
            dimension="el_rank_max_pct",
            value=v,
            label=f"EL{'≤' if base.el_lte else '<'}{v}",
            spec=replace(base, el_max=v),
        )
    for v in ic50_sweep:
        add_row(
            analysis="sensitivity",
            dimension="ic50_max_nm",
            value=v,
            label=f"IC50<{v:g}nM",
            spec=replace(base, ic50_max_nm=v),
        )

    out = pd.DataFrame(rows_long)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.out_stem
    csv_path = args.out_dir / f"{stem}.csv"
    out.to_csv(csv_path, index=False)

    b0 = out[out["analysis"] == "baseline"].iloc[0]
    b_y = max(int(b0[y_col]), 1)
    fc = out.copy()
    fc[f"{y_col}_pct_of_baseline_pct"] = 100.0 * fc[y_col].astype(float) / float(b_y)
    fc[f"{y_col}_fold_vs_baseline"] = fc[y_col].astype(float) / float(b_y)
    fc_path = args.out_dir / f"{stem}_fold_change_vs_baseline.csv"
    fc.to_csv(fc_path, index=False)

    loo = out[out["analysis"] == "leave_one_out"].copy()
    base_row = out[out["analysis"] == "baseline"].iloc[0]
    x_labels = ["baseline\n(full SB)"] + [s.replace("relax_", "").replace("_", "\n") for s in loo["dimension"]]
    x = np.arange(len(x_labels))
    y_bars = [int(base_row[y_col])] + loo[y_col].astype(int).tolist()

    def line_panel(dim: str, title: str, ax: Axes) -> None:
        sub = out[(out["analysis"] == "sensitivity") & (out["dimension"] == dim)].sort_values("threshold_value")
        ax.plot(sub["threshold_value"], sub[y_col], "o-", color=pal.SIG_LNC, lw=1.2, markersize=5)
        ax.set_xlabel(title)
        ax.set_ylabel(y_label)
        ax.grid(alpha=0.3)

    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0], hspace=0.4, wspace=0.32)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.bar(x, y_bars, color=pal.SIG_LNC, edgecolor="black", linewidth=0.35)
    ax0.set_xticks(x)
    ax0.set_xticklabels(x_labels, fontsize=9)
    ax0.set_ylabel(y_label)
    ax0.set_title("TTN merged IEDB+NetMHC — leave-one-filter-out (others at baseline)")
    ax0.grid(axis="y", alpha=0.3)

    ax1 = fig.add_subplot(gs[1, 0])
    line_panel("immuno_min", "Immunogenicity cutoff (iedb_score >)", ax1)
    ax2 = fig.add_subplot(gs[1, 1])
    line_panel("processing_min", "Processing cutoff (iedb_processing_score >)", ax2)
    ax3 = fig.add_subplot(gs[2, 0])
    line_panel("el_rank_max_pct", "EL %rank cutoff (strict < unless --el-rank-lte)", ax3)
    ax4 = fig.add_subplot(gs[2, 1])
    sub_ic = out[(out["analysis"] == "sensitivity") & (out["dimension"] == "ic50_max_nm")].sort_values(
        "threshold_value"
    )
    ax4.plot(sub_ic["threshold_value"], sub_ic[y_col], "o-", color=pal.SIG_LNC, lw=1.2, markersize=5)
    ax4.set_xscale("log")
    ax4.set_xlabel("IC50 max (nM), strict < cutoff")
    ax4.set_ylabel(y_label)
    ax4.grid(alpha=0.3)

    fig.suptitle(
        f"TTN-AS1 merged IEDB+NetMHC — SB 1D sensitivity + LOO ({args.plot_metric}; sb_mode={args.sb_mode})",
        fontsize=11,
        y=0.995,
    )
    png_path = args.out_dir / f"{stem}.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {csv_path}")
    print(f"Wrote {fc_path}")
    print(f"Wrote {png_path}")
    print(f"[info] IC50 source: {iedb_ic50_col or 'local BA_score'}")
    print(f"[info] baseline: {base}")


if __name__ == "__main__":
    main()
