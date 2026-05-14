#!/usr/bin/env python3
"""
**Figure 6 supplement (TTN-AS1):** Cartesian SB threshold grid on **merged** NetMHC wide XLS + IEDB CSV.

Mirrors the *multi-parameter* exploration in ``plot_fig5_netmhc_sb_combination_grid.py``, but for the
single-parent TTN-AS1 table built the same way as ``plot_figure6_ttn_as1_allele_coverage.py`` with
``--gating iedb_sb``: long rows keyed by ``stable_key``, merged to IEDB, then ``sb_mask_spec`` over
immuno / processing / EL / IC50 (IEDB IC50 column when available, else BA-derived gate).

**Figures (conceptually aligned with Fig 6 / Fig 5 supplement):** in addition to CSV grids, this
script writes a **five-panel sweep suite** — **(a)** four fold-vs-baseline heatmaps (imm×IC50,
proc×IC50, imm×proc, EL×IC50, each with the other two gates fixed at baseline), **(b–e)** 1D sweeps
along each gate at baseline on the other three (same design as the cohort 1D sensitivity script, but
values are taken from the Cartesian grid). Also writes one **combined** PNG with panels labeled
(a)–(e).

**Not** the same as ``plot_figure6_ttn_as1_sb_sensitivity.py`` (NetMHC-only BA_rank / IC50-from-BA sweeps).
For **1D sweeps + LOO** on this same merged TTN table, use ``netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py``.

Default output: ``data/netmhc/figures/fig6_ttn_merged_iedb_sb_combinations/`` unless ``--out-dir`` is set.
``generate_netmhc_fig5_fig6_supplement.py`` points here to
``figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_merged_iedb_cartesian_sb_grid/``.
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
import itertools

import figure_palettes as pal  # noqa: E402
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm
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
    spec_label,
    spec_profile_id,
)

import plot_figure6_ttn_as1_allele_coverage as ttn6  # noqa: E402


def parse_float_list(s: str) -> list[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def cohort_row_stats(df: pd.DataFrame, m: pd.Series) -> dict[str, int]:
    sub = df.loc[m]
    return {
        "n_rows": int(len(sub)),
        "n_unique_peptides": int(sub["Peptide"].nunique()),
    }


def _close(s: pd.Series, val: float, *, rtol: float = 1e-9) -> pd.Series:
    return np.isclose(s.astype(float), float(val), rtol=rtol, atol=0.0)


def _filter_baseline_slice(fc: pd.DataFrame, baseline: SBSpec, *, vary: str | None) -> pd.DataFrame:
    """Rows matching baseline on every axis except ``vary`` (imm|proc|el|ic50|None for full grid)."""
    m = np.ones(len(fc), dtype=bool)
    if vary != "imm":
        m &= _close(fc["imm_min"], baseline.imm_min)
    if vary != "proc":
        m &= _close(fc["proc_min"], baseline.proc_min)
    if vary != "el":
        m &= _close(fc["el_max_pct"], baseline.el_max)
    if vary != "ic50":
        m &= _close(fc["ic50_max_nm"], baseline.ic50_max_nm)
    return fc.loc[m].copy()


def _pivot_fold(
    fc: pd.DataFrame,
    *,
    row_key: str,
    col_key: str,
    fix: dict[str, float],
    fold_col: str,
) -> tuple[np.ndarray, list[float], list[float]] | None:
    sub = fc.copy()
    for fk, fv in fix.items():
        sub = sub.loc[_close(sub[fk], fv)]
    if sub.empty:
        return None
    try:
        pv = sub.pivot_table(index=row_key, columns=col_key, values=fold_col, aggfunc="first")
    except ValueError:
        return None
    if pv.size == 0 or pv.shape[0] < 1 or pv.shape[1] < 1:
        return None
    mat = np.nan_to_num(pv.to_numpy(dtype=float), nan=1.0)
    rows = [float(x) for x in pv.index]
    cols = [float(x) for x in pv.columns]
    return mat, rows, cols


def _im_fold_heatmap(
    ax,
    mat: np.ndarray,
    row_labels: list[float],
    col_labels: list[float],
    *,
    xlabel: str,
    ylabel: str,
    title: str,
    vlim: float,
) -> None:
    norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=float(vlim))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", norm=norm, origin="lower")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels([f"{c:g}" for c in col_labels], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels([f"{r:g}" for r in row_labels], fontsize=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Fold vs baseline (unique 9-mers)")


def _collect_four_matrices(fc: pd.DataFrame, baseline: SBSpec, fold_col: str) -> list[tuple[np.ndarray, list[float], list[float], str, str, str]]:
    b_imm, b_proc, b_el, b_ic = baseline.imm_min, baseline.proc_min, baseline.el_max, baseline.ic50_max_nm
    specs: list[tuple[str, str, str, str, dict[str, float], str]] = [
        ("imm_min", "ic50_max_nm", "Imm (iedb_score >)", "IC50 max (nM)", {"proc_min": b_proc, "el_max_pct": b_el}, "proc,EL@baseline"),
        ("proc_min", "ic50_max_nm", "Proc (score >)", "IC50 max (nM)", {"imm_min": b_imm, "el_max_pct": b_el}, "imm,EL@baseline"),
        ("imm_min", "proc_min", "Imm (iedb_score >)", "Proc (score >)", {"el_max_pct": b_el, "ic50_max_nm": b_ic}, "EL,IC50@baseline"),
        ("el_max_pct", "ic50_max_nm", "EL %rank cap", "IC50 max (nM)", {"imm_min": b_imm, "proc_min": b_proc}, "imm,proc@baseline"),
    ]
    out: list[tuple[np.ndarray, list[float], list[float], str, str, str]] = []
    for rk, ck, xl, yl, fixd, note in specs:
        pv = _pivot_fold(fc, row_key=rk, col_key=ck, fix=fixd, fold_col=fold_col)
        if pv is None:
            continue
        mat, rr, cc = pv
        out.append((mat, rr, cc, xl, yl, note))
    return out


def _global_vmax(mats: list[np.ndarray]) -> float:
    vals = [float(np.nanmax(m)) for m in mats if m.size]
    return max(2.0, max(vals) if vals else 2.0)


def _line_sweep(
    ax,
    fc: pd.DataFrame,
    baseline: SBSpec,
    vary: str,
    *,
    fold_u: str,
    fold_r: str,
) -> None:
    key_map = {"imm": "imm_min", "proc": "proc_min", "el": "el_max_pct", "ic50": "ic50_max_nm"}
    xk = key_map[vary]
    sub = _filter_baseline_slice(fc, baseline, vary=vary).sort_values(xk)
    if sub.empty:
        ax.set_visible(False)
        return
    x = sub[xk].to_numpy(dtype=float)
    ax.plot(x, sub[fold_u], "o-", color=pal.SIG_LNC, lw=1.2, markersize=5, label="Unique 9-mers (fold)")
    ax.plot(x, sub[fold_r], "s--", color="0.35", lw=1.0, markersize=4, label="SB rows (fold)")
    if vary == "ic50":
        ax.set_xscale("log")
    ax.axhline(1.0, color="0.6", ls=":", lw=0.8)
    xlabel_map = {
        "imm_min": "Immunogenicity min (iedb_score >)",
        "proc_min": "Processing min",
        "el_max_pct": "EL %rank cap",
        "ic50_max_nm": "IC50 max (nM)",
    }
    ax.set_xlabel(xlabel_map.get(xk, xk))
    ax.set_ylabel("Fold vs baseline")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=7)


def _draw_panel_a(out_path: Path, fc: pd.DataFrame, baseline: SBSpec, fold_col: str, *, suptitle: str) -> None:
    mats_wrapped = _collect_four_matrices(fc, baseline, fold_col)
    if not mats_wrapped:
        return
    vmax = _global_vmax([m[0] for m in mats_wrapped])
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9), squeeze=False)
    titles = [
        "(a) imm × IC50",
        "(a) proc × IC50",
        "(a) imm × proc",
        "(a) EL × IC50",
    ]
    for ax, (mat, rr, cc, xl, yl, note), tlab in zip(axes.ravel(), mats_wrapped, titles):
        _im_fold_heatmap(ax, mat, rr, cc, xlabel=xl, ylabel=yl, title=f"{tlab}\n({note})", vlim=vmax)
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _draw_line_panel(out_path: Path, fc: pd.DataFrame, baseline: SBSpec, vary: str, title: str, supt: str) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    _line_sweep(ax, fc, baseline, vary, fold_u="n_unique_9mers_fold_vs_baseline", fold_r="n_sb_rows_fold_vs_baseline")
    ax.set_title(title, fontsize=10)
    fig.suptitle(supt, fontsize=9, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _draw_combined_suite(
    out_path: Path,
    fc: pd.DataFrame,
    baseline: SBSpec,
    fold_col: str,
    *,
    suptitle: str,
    ic50_col_note: str,
) -> None:
    mats_wrapped = _collect_four_matrices(fc, baseline, fold_col)
    if len(mats_wrapped) < 4:
        print(
            "[plot_fig6_ttn_merged_iedb_sb_combination_grid] Skipping combined suite PNG: "
            "need four 2D slices (check grids include baseline-aligned combinations).",
            flush=True,
        )
        return
    fig = plt.figure(figsize=(11, 14))
    gs = gridspec.GridSpec(4, 2, figure=fig, height_ratios=[1.0, 1.0, 0.85, 0.85], hspace=0.45, wspace=0.35)

    vmax = _global_vmax([m[0] for m in mats_wrapped]) if mats_wrapped else 2.0
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    panel_titles = [
        "(a) imm × IC50 @ baseline proc, EL",
        "(a) proc × IC50 @ baseline imm, EL",
        "(a) imm × proc @ baseline EL, IC50",
        "(a) EL × IC50 @ baseline imm, proc",
    ]
    for k, (mat, rr, cc, xl, yl, note) in enumerate(mats_wrapped):
        r, c = positions[k]
        ax = fig.add_subplot(gs[r, c])
        _im_fold_heatmap(ax, mat, rr, cc, xlabel=xl, ylabel=yl, title=f"{panel_titles[k]}\n({note})", vlim=vmax)

    ax_b = fig.add_subplot(gs[2, 0])
    _line_sweep(ax_b, fc, baseline, "imm", fold_u="n_unique_9mers_fold_vs_baseline", fold_r="n_sb_rows_fold_vs_baseline")
    ax_b.set_title("(b) Sweep immunogenicity min (others @ baseline)", fontsize=10)

    ax_c = fig.add_subplot(gs[2, 1])
    _line_sweep(ax_c, fc, baseline, "proc", fold_u="n_unique_9mers_fold_vs_baseline", fold_r="n_sb_rows_fold_vs_baseline")
    ax_c.set_title("(c) Sweep processing min (others @ baseline)", fontsize=10)

    ax_d = fig.add_subplot(gs[3, 0])
    _line_sweep(ax_d, fc, baseline, "el", fold_u="n_unique_9mers_fold_vs_baseline", fold_r="n_sb_rows_fold_vs_baseline")
    ax_d.set_title("(d) Sweep EL %rank cap (others @ baseline)", fontsize=10)

    ax_e = fig.add_subplot(gs[3, 1])
    _line_sweep(ax_e, fc, baseline, "ic50", fold_u="n_unique_9mers_fold_vs_baseline", fold_r="n_sb_rows_fold_vs_baseline")
    ax_e.set_title("(e) Sweep IC50 max nM — log x (others @ baseline)", fontsize=10)

    fig.suptitle(f"{suptitle}\nIC50 column: {ic50_col_note}", fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.07, right=0.96, top=0.93, bottom=0.05, hspace=0.42, wspace=0.38)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--netmhc-xls",
        type=Path,
        default=Path("data/netmhc/netmhcpan_ttn_as1_108065.xls"),
    )
    ap.add_argument(
        "--iedb-csv",
        type=Path,
        default=ttn6.TTN_IEDB_SYNTHETIC_CSV,
    )
    ap.add_argument(
        "--iedb-parent-input-seq-id",
        type=str,
        default=ttn6.TTN_IEDB_PARENT_SEQ_ID_DEFAULT,
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/netmhc/figures/fig6_ttn_merged_iedb_sb_combinations"),
        help="Output folder for CSV + PNG.",
    )
    ap.add_argument("--el-rank-lte", action="store_true")
    ap.add_argument("--baseline-imm", type=float, default=FIG5_IEDB_IMM_MIN_DEFAULT)
    ap.add_argument("--baseline-proc", type=float, default=FIG5_IEDB_PROC_MIN_DEFAULT)
    ap.add_argument("--baseline-el-max", type=float, default=FIG5_IEDB_EL_RANK_MAX_DEFAULT)
    ap.add_argument("--baseline-ic50-nm", type=float, default=FIG5_IEDB_IC50_MAX_NM_DEFAULT)
    ap.add_argument("--imm-grid", type=str, default="0.05,0.1,0.15")
    ap.add_argument("--proc-grid", type=str, default="1.25,1.5,1.75")
    ap.add_argument("--el-grid", type=str, default="1,2,5,10")
    ap.add_argument("--ic50-grid", type=str, default="150,500,1000")
    ap.add_argument(
        "--no-suite-pngs",
        action="store_true",
        help="Skip multi-panel sweep PNGs (CSV outputs are still written).",
    )
    args = ap.parse_args()

    if not args.netmhc_xls.is_file():
        raise SystemExit(f"Missing NetMHC XLS: {args.netmhc_xls}")
    if not args.iedb_csv.is_file():
        raise SystemExit(f"Missing IEDB CSV: {args.iedb_csv}")

    alleles, starts, peps, ba, ba_rank, el_rank = ttn6.parse_wide_netmhc_xls_rows(args.netmhc_xls)
    long_df = ttn6.build_ttn_long_for_iedb_merge(
        starts, peps, alleles, ba, ba_rank, el_rank, str(args.iedb_parent_input_seq_id).strip()
    )
    merged = ttn6.merge_iedb_csv(long_df, args.iedb_csv)

    hdr = sorted(set(merged.columns))
    iedb_ic50_col = pick_iedb_ic50_column(hdr)

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
        ms = sb_mask_spec(merged, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
        st = cohort_row_stats(merged, ms)
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
                "n_sb_rows": st["n_rows"],
                "n_unique_9mers": st["n_unique_peptides"],
                "ic50_source": iedb_ic50_col or "local_BA_score",
            }
        )

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    grid_df = pd.DataFrame(grid_rows)
    grid_csv = out_dir / "fig6_ttn_merged_iedb_counts_cartesian_grid.csv"
    grid_df.to_csv(grid_csv, index=False)

    mask_base = (
        (grid_df["imm_min"] == baseline.imm_min)
        & (grid_df["proc_min"] == baseline.proc_min)
        & (grid_df["el_max_pct"] == baseline.el_max)
        & (grid_df["ic50_max_nm"] == baseline.ic50_max_nm)
    )
    if mask_base.any():
        br = grid_df.loc[mask_base].iloc[0]
        b_rows = max(int(br["n_sb_rows"]), 1)
        b_uni = max(int(br["n_unique_9mers"]), 1)
    else:
        bm0 = ba_score_min_for_ic50_lt(baseline.ic50_max_nm)
        mb = sb_mask_spec(merged, baseline, ba_min=bm0, iedb_ic50_col=iedb_ic50_col)
        bst = cohort_row_stats(merged, mb)
        b_rows = max(int(bst["n_rows"]), 1)
        b_uni = max(int(bst["n_unique_9mers"]), 1)

    fc = grid_df.copy()
    fc["n_sb_rows_pct_of_baseline"] = 100.0 * fc["n_sb_rows"].astype(float) / float(b_rows)
    fc["n_unique_9mers_pct_of_baseline"] = 100.0 * fc["n_unique_9mers"].astype(float) / float(b_uni)
    fc["n_sb_rows_fold_vs_baseline"] = fc["n_sb_rows"].astype(float) / float(b_rows)
    fc["n_unique_9mers_fold_vs_baseline"] = fc["n_unique_9mers"].astype(float) / float(b_uni)
    fc_path = out_dir / "fig6_ttn_merged_iedb_fold_change_vs_baseline_cartesian_grid.csv"
    fc.to_csv(fc_path, index=False)

    ic50_note = iedb_ic50_col or "local BA"
    supt = "Fig 6 supplement — TTN merged IEDB+NetMHC: Cartesian SB grid (fold vs baseline)"

    if not args.no_suite_pngs:
        stem = "fig6_ttn_merged_iedb_cartesian_sweep"
        _draw_panel_a(
            out_dir / f"{stem}_a_heatmaps_fold_unique.png",
            fc,
            baseline,
            "n_unique_9mers_fold_vs_baseline",
            suptitle=f"{supt} — (a) heatmaps",
        )
        _draw_line_panel(
            out_dir / f"{stem}_b_imm.png",
            fc,
            baseline,
            "imm",
            "(b) Immunogenicity min",
            supt,
        )
        _draw_line_panel(
            out_dir / f"{stem}_c_proc.png",
            fc,
            baseline,
            "proc",
            "(c) Processing min",
            supt,
        )
        _draw_line_panel(
            out_dir / f"{stem}_d_el.png",
            fc,
            baseline,
            "el",
            "(d) EL %rank cap",
            supt,
        )
        _draw_line_panel(
            out_dir / f"{stem}_e_ic50.png",
            fc,
            baseline,
            "ic50",
            "(e) IC50 max (nM)",
            supt,
        )
        _draw_combined_suite(
            out_dir / f"{stem}_abcde_combined.png",
            fc,
            baseline,
            "n_unique_9mers_fold_vs_baseline",
            suptitle=supt,
            ic50_col_note=ic50_note,
        )

    print(f"Wrote {grid_csv} ({len(grid_df)} combinations)")
    print(f"Wrote {fc_path}")
    print(f"[baseline] n_sb_rows={b_rows}, n_unique_9mers={b_uni}")
    if not args.no_suite_pngs:
        print(f"Wrote sweep suite PNGs under {out_dir} (see fig6_ttn_merged_iedb_cartesian_sweep_*.png)")


if __name__ == "__main__":
    main()
