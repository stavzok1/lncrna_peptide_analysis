"""
**Catalog Figure 6 — sensitivity (TTN-AS1):** how counts depend on the strong-binder (SB) rule.

Complements ``plot_figure6_ttn_as1_allele_coverage.py`` (main multi-panel figure). This
script **does not** redraw panels A–E; it tabulates summary metrics across SB parameter
grids and writes a small multi-panel diagnostic figure.

Default inputs match Figure 6: ``data/netmhc/netmhcpan_ttn_as1_108065.xls`` and the same
79-aa parent sequence (unless ``--parent-fasta`` is set). These sweeps are **NetMHC-only**
(BA_rank / IC50 from BA_score / optional EL); they do **not** apply IEDB immunogenicity or
processing thresholds (those are Figure 5 cohort / merged-TSV only).

Outputs (default): ``data/netmhc/figures/fig6_ttn_as1_sensitivity/``

  - ``ttn_as1_sb_sensitivity.csv`` — one row per configuration.
  - ``ttn_as1_sb_sensitivity_fold_change_vs_baseline.csv`` — adds %% and fold vs default SB.
  - ``fig6_ttn_as1_sb_sensitivity_overview.png`` — line plots of key metrics vs cutoffs.
    IC50 panels default to **no EL gate** (single curves, linear IC50 axis), matching the
    historical figure; use ``--overview-dual-ic50-el-lines`` to overlay the +EL variant.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_REPO / "manuscript")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES

ROOT = REPO_ROOT


import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_figure6_ttn_as1_allele_coverage as ttn

import figure_palettes as pal

NET = ROOT / "data" / "netmhc"
DEFAULT_OUT = NET / "figures" / "fig6_ttn_as1_sensitivity"


def parse_float_list(s: str) -> list[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def load_parent_sequence(parent_fasta: Path | None) -> str:
    if parent_fasta is None:
        return (
            "ISATDRICENTSMSRLGIILRHHLASPASHFKMIANDSTSSITDWLIPLYFHAVPGGQCDNWSARRTRNFEWILGYSRL"
        )
    from Bio import SeqIO

    rec = next(SeqIO.parse(parent_fasta, "fasta"))
    return str(rec.seq).upper().strip()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--netmhc-xls",
        type=Path,
        default=NET / "netmhcpan_ttn_as1_108065.xls",
    )
    ap.add_argument("--parent-fasta", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--ba-rank-grid",
        type=str,
        default="0.25,0.5,1.0,2.0",
        help="Comma-separated BA_rank %% cutoffs (SB if BA_rank <= value). Used when --include-ba-rank.",
    )
    ap.add_argument(
        "--ic50-grid",
        type=str,
        default="50,150,500,1000",
        help="Comma-separated IC50 (nM) upper bounds when using --include-ic50-criterion.",
    )
    ap.add_argument(
        "--no-ba-rank-sweep",
        action="store_true",
        help="Skip BA_rank sweeps (keep IC50 sweeps unless --no-ic50-sweep).",
    )
    ap.add_argument(
        "--no-ic50-sweep",
        action="store_true",
        help="Skip IC50-from-BA_score sweeps (keep BA_rank unless --no-ba-rank-sweep).",
    )
    ap.add_argument(
        "--require-el-grid",
        type=str,
        default="0,1",
        help="Comma-separated 0/1: whether to require EL_rank <= cutoff (1=yes). "
        "Drives rows in the CSV; overview IC50 panels default to require_el_rank=False only "
        "(see --overview-dual-ic50-el-lines).",
    )
    ap.add_argument(
        "--overview-dual-ic50-el-lines",
        action="store_true",
        help="On the overview PNG only: plot two IC50 curves (no EL vs +EL_rank gate). "
        "Default is a single curve (no EL), matching the historical figure layout.",
    )
    args = ap.parse_args()

    full = load_parent_sequence(args.parent_fasta)
    alleles, starts, peps, ba, ba_rank, el_rank = ttn.parse_wide_netmhc_xls_rows(args.netmhc_xls)

    req_el_flags = [bool(int(x)) for x in parse_float_list(args.require_el_grid) if int(float(x)) in (0, 1)]
    if not req_el_flags:
        req_el_flags = [False, True]

    rows: list[dict[str, object]] = []
    for require_el in req_el_flags:
        if not args.no_ba_rank_sweep:
            for pct in parse_float_list(args.ba_rank_grid):
                sb = ttn.strong_binder_mask(
                    ba,
                    ba_rank,
                    "ba_rank",
                    150.0,
                    pct,
                    el_rank=el_rank,
                    require_el_rank=require_el,
                    el_rank_pct=None,
                )
                epitope_cov, allele_cov = ttn.per_position_metrics(full, starts, peps, sb)
                sb_u = sorted(set(ttn.collect_sb_epitopes(peps, sb, None)))
                rows.append(
                    {
                        "sb_criterion": "ba_rank",
                        "ba_rank_pct_cutoff": pct,
                        "ic50_nm_cutoff": None,
                        "require_el_rank": require_el,
                        "el_rank_uses_same_cutoff_as_ba_rank": True,
                        "n_unique_sb_9mers": len(sb_u),
                        "n_positions_allele_cov_gt0": int((allele_cov > 0).sum()),
                        "frac_positions_with_any_sb_allele": float((allele_cov > 0).mean()),
                        "max_allele_diversity_per_position": int(allele_cov.max()),
                        "mean_allele_diversity_per_position": float(allele_cov.mean()),
                        "max_distinct_epitopes_per_position": int(epitope_cov.max()),
                    }
                )
        if not args.no_ic50_sweep:
            for ic50 in parse_float_list(args.ic50_grid):
                sb = ttn.strong_binder_mask(
                    ba,
                    ba_rank,
                    "ic50",
                    ic50,
                    0.5,
                    el_rank=el_rank,
                    require_el_rank=require_el,
                    el_rank_pct=None,
                )
                epitope_cov, allele_cov = ttn.per_position_metrics(full, starts, peps, sb)
                sb_u = sorted(set(ttn.collect_sb_epitopes(peps, sb, None)))
                rows.append(
                    {
                        "sb_criterion": "ic50_from_ba_score",
                        "ba_rank_pct_cutoff": None,
                        "ic50_nm_cutoff": ic50,
                        "require_el_rank": require_el,
                        "el_rank_uses_same_cutoff_as_ba_rank": False,
                        "n_unique_sb_9mers": len(sb_u),
                        "n_positions_allele_cov_gt0": int((allele_cov > 0).sum()),
                        "frac_positions_with_any_sb_allele": float((allele_cov > 0).mean()),
                        "max_allele_diversity_per_position": int(allele_cov.max()),
                        "mean_allele_diversity_per_position": float(allele_cov.mean()),
                        "max_distinct_epitopes_per_position": int(epitope_cov.max()),
                    }
                )

    out = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "ttn_as1_sb_sensitivity.csv"
    out.to_csv(csv_path, index=False)

    # Baseline: BA_rank 0.5, no EL requirement (matches plot_figure6 defaults)
    base_mask = (out["sb_criterion"] == "ba_rank") & (out["ba_rank_pct_cutoff"] == 0.5) & (~out["require_el_rank"])
    if base_mask.any():
        b = out.loc[base_mask].iloc[0]
    else:
        b = out.iloc[0]
    b_n = max(int(b["n_unique_sb_9mers"]), 1)
    b_cov = max(float(b["frac_positions_with_any_sb_allele"]), 1e-9)
    fc = out.copy()
    fc["n_unique_pct_of_baseline_pct"] = 100.0 * fc["n_unique_sb_9mers"].astype(float) / float(b_n)
    fc["n_unique_fold_vs_baseline"] = fc["n_unique_sb_9mers"].astype(float) / float(b_n)
    fc["frac_covered_pct_of_baseline_pct"] = 100.0 * fc["frac_positions_with_any_sb_allele"].astype(float) / b_cov
    fc_path = args.out_dir / "ttn_as1_sb_sensitivity_fold_change_vs_baseline.csv"
    fc.to_csv(fc_path, index=False)

    # Overview figure: BA_rank sweeps (no EL vs +EL) and IC50 sweeps.  The CSV may contain
    # duplicate IC50 x-values when --require-el-grid includes 1 (EL gate on).  By default the
    # overview plots IC50 for require_el_rank=False only (single lines, linear x), matching the
    # historical PNG; pass --overview-dual-ic50-el-lines for two IC50 curves + legend.
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    br = out[(out["sb_criterion"] == "ba_rank") & (~out["require_el_rank"])].sort_values("ba_rank_pct_cutoff")
    if len(br):
        axes[0, 0].plot(br["ba_rank_pct_cutoff"], br["n_unique_sb_9mers"], "o-", color=pal.OI_SKY_BLUE)
        axes[0, 0].set_xlabel("BA_rank cutoff (% rank, lower=stricter)")
        axes[0, 0].set_ylabel("Unique SB 9-mers")
        axes[0, 0].set_title("TTN-AS1: SB count vs BA_rank (no EL filter)")
        axes[0, 0].grid(alpha=0.3)
    br_el = out[(out["sb_criterion"] == "ba_rank") & (out["require_el_rank"])].sort_values("ba_rank_pct_cutoff")
    if len(br_el):
        axes[0, 1].plot(br_el["ba_rank_pct_cutoff"], br_el["n_unique_sb_9mers"], "s-", color=pal.OI_VERMILLION)
        axes[0, 1].set_xlabel("BA_rank cutoff")
        axes[0, 1].set_ylabel("Unique SB 9-mers")
        axes[0, 1].set_title("TTN-AS1: SB count vs BA_rank (+ EL_rank <= same cutoff)")
        axes[0, 1].grid(alpha=0.3)
    ic_base = out[out["sb_criterion"] == "ic50_from_ba_score"]
    if len(ic_base):
        if args.overview_dual_ic50_el_lines:
            for req_el, sty, lab in (
                (False, ("o", "-"), "IC50 SB, no EL gate"),
                (True, ("^", "--"), "IC50 SB + EL_rank ≤ 0.5 %rank"),
            ):
                ic = ic_base.loc[ic_base["require_el_rank"] == req_el].sort_values("ic50_nm_cutoff")
                if len(ic) == 0:
                    continue
                mk, ls = sty
                axes[1, 0].plot(
                    ic["ic50_nm_cutoff"],
                    ic["n_unique_sb_9mers"],
                    marker=mk,
                    linestyle=ls,
                    color=pal.OI_BLUISH_GREEN,
                    label=lab,
                )
                axes[1, 1].plot(
                    ic["ic50_nm_cutoff"],
                    ic["frac_positions_with_any_sb_allele"],
                    marker=mk,
                    linestyle=ls,
                    color=pal.OI_REDDISH_PURPLE,
                    label=lab,
                )
            axes[1, 0].legend(fontsize=7, loc="best")
            axes[1, 1].legend(fontsize=7, loc="best")
        else:
            ic = ic_base.loc[~ic_base["require_el_rank"]].sort_values("ic50_nm_cutoff")
            if len(ic) == 0:
                ic = ic_base.sort_values("ic50_nm_cutoff")
            axes[1, 0].plot(ic["ic50_nm_cutoff"], ic["n_unique_sb_9mers"], "o-", color=pal.OI_BLUISH_GREEN)
            axes[1, 1].plot(ic["ic50_nm_cutoff"], ic["frac_positions_with_any_sb_allele"], "s-", color=pal.OI_REDDISH_PURPLE)
        axes[1, 0].set_xlabel("IC50 cutoff (nM)")
        axes[1, 0].set_ylabel("Unique SB 9-mers")
        axes[1, 0].set_title("TTN-AS1: SB count vs IC50 (from BA_score)")
        axes[1, 0].grid(alpha=0.3)
        axes[1, 1].set_xlabel("IC50 cutoff (nM)")
        axes[1, 1].set_ylabel("Fraction positions with >=1 SB allele")
        axes[1, 1].set_ylim(0, 1.05)
        axes[1, 1].set_title("TTN-AS1: positional coverage vs IC50")
        axes[1, 1].grid(alpha=0.3)
    else:
        axes[1, 0].set_visible(False)
        axes[1, 1].set_visible(False)

    fig.suptitle("TTN-AS1 (smPEP 108065): SB definition sensitivity", fontsize=12)
    fig.tight_layout()
    png = args.out_dir / "fig6_ttn_as1_sb_sensitivity_overview.png"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {csv_path} ({len(out)} rows)")
    print(f"Wrote {fc_path}")
    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
