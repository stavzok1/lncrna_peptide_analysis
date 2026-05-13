"""
Figure 5B / 5C: distribution of epitope promiscuity across HLA alleles (NetMHCpan).

For each **unique 9-mer** (epitope), count how many alleles it is a **strong binder**
(IC50 < threshold nM, derived from ``BA_score`` in wide XLS — same convention as
``plot_netmhc_epitopes_vs_hla_frequency.py``). Multiple XLS rows with the same 9-mer
(e.g. duplicated sliding-window entries) are merged with **logical OR** per allele.

- **Default / Fig 5B:** ``data/netmhc/netmhcpan_sig_lnc.xls`` (lncRNA significant 9-mers).
- **``--coding-control`` / Fig 5C:** ``data/netmhc/netmhcpan_coding_proportional_whole.xls`` (Fig 5 coding
  sample: **whole** proteome proteins matched to the significant-MP length mix; see
  ``prepare_netmhc_tr_vs_coding_epitopes.py --coding-control-mode proportional_whole``).

Bar chart: x = number of alleles per epitope (1 … max). Default **y** = **SB epitopes per allele** (bar height is
``n_alleles_per_epitope`` × count of epitopes in that bin); use ``--histogram-metric unique``
for counts of unique 9-mers only.

PNG and CSV are also copied to repo-root ``figures/``.
"""
from __future__ import annotations

import argparse
import shutil
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from plot_netmhc_epitopes_vs_hla_frequency import ba_score_to_ic50_nm, parse_wide_netmhc_xls

ROOT = Path(__file__).resolve().parent
import figure_palettes as pal

NET = ROOT / "data" / "netmhc"
FIGS = NET / "figures"
REPO_FIGURES = ROOT / "figures"


def _mirror_to_repo_figures(*paths: Path) -> None:
    REPO_FIGURES.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.is_file():
            shutil.copy2(p, REPO_FIGURES / p.name)


def _coding_dataset_label(xls: Path) -> str:
    name = xls.name.lower()
    if "proportional_whole" in name:
        return "Coding control (whole proteins, proportional sample)"
    if "coding_control" in name:
        return "Coding control (random proteome fragments)"
    return "Coding control"
def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 5B/5C: epitope sharing across alleles (SB by IC50).")
    ap.add_argument(
        "--coding-control",
        action="store_true",
        help="Use proteome control NetMHC XLS (Fig 5C): netmhcpan_coding_proportional_whole.xls.",
    )
    ap.add_argument("--netmhc-xls", type=Path, default=None, help="Wide NetMHCpan XLS (default: sig or control).")
    ap.add_argument("--ic50-nm", type=float, default=150.0)
    ap.add_argument(
        "--histogram-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Bar heights: y-axis SB epitopes per allele (instances) vs unique 9-mers per bin (default instances).",
    )
    ap.add_argument("--out-png", type=Path, default=None)
    ap.add_argument("--out-csv", type=Path, default=None)
    ap.add_argument(
        "--dataset-label",
        default=None,
        help="Legend text on plot (defaults: Tr-lncRNA-MPs vs proteome-sampled control).",
    )
    ap.add_argument(
        "--no-repo-mirror",
        action="store_true",
        help="Do not copy outputs into repo-root figures/.",
    )
    args = ap.parse_args()

    if args.netmhc_xls is None:
        args.netmhc_xls = (
            NET / "netmhcpan_coding_proportional_whole.xls"
            if args.coding_control
            else NET / "netmhcpan_sig_lnc.xls"
        )
    if args.out_png is None:
        args.out_png = (
            FIGS / "fig5c_epitope_sharing_across_alleles.png"
            if args.coding_control
            else FIGS / "fig5b_epitope_sharing_across_alleles.png"
        )
    if args.out_csv is None:
        args.out_csv = (
            FIGS / "fig5c_epitope_sharing_across_alleles.csv"
            if args.coding_control
            else FIGS / "fig5b_epitope_sharing_across_alleles.csv"
        )
    if args.dataset_label is None:
        args.dataset_label = (
            _coding_dataset_label(args.netmhc_xls) if args.coding_control else "Tr-lncRNA-MPs"
        )
    alleles, ba_mat, peps = parse_wide_netmhc_xls(args.netmhc_xls)
    ic50 = np.vectorize(ba_score_to_ic50_nm)(ba_mat)
    strong = ic50 < args.ic50_nm

    # OR across rows that share the same 9-mer sequence
    agg: dict[str, np.ndarray] = {}
    for pep, row in zip(peps, strong):
        if pep not in agg:
            agg[pep] = row.copy()
        else:
            agg[pep] |= row

    n_alleles_per_epitope: list[int] = []
    for mask in agg.values():
        n = int(mask.sum())
        if n >= 1:
            n_alleles_per_epitope.append(n)

    if not n_alleles_per_epitope:
        raise SystemExit("No epitopes passed IC50 threshold for any allele.")

    ctr = Counter(n_alleles_per_epitope)
    kmax = max(ctr)
    xs = np.arange(1, kmax + 1, dtype=int)
    heights_unique = np.array([ctr.get(int(k), 0) for k in xs], dtype=int)
    rows_in_bin = np.array(xs * heights_unique, dtype=int)
    heights = rows_in_bin if args.histogram_metric == "instances" else heights_unique
    total_unique = int(heights_unique.sum())
    total_sb_rows = int(rows_in_bin.sum())

    FIGS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "n_alleles_per_epitope": xs,
            "n_epitopes": heights_unique,
            "n_sb_rows_in_bin": rows_in_bin,
            "total_epitopes": total_unique,
            "total_sb_row_instances": total_sb_rows,
        }
    ).to_csv(args.out_csv, index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=150)
    bar_c = pal.SIG_LNC
    ax.bar(xs, heights, color=bar_c, edgecolor="white", linewidth=0.4, width=0.92)
    ax.set_xlabel("Number of alleles per epitope")
    ax.set_ylabel(
        "SB epitopes per allele" if args.histogram_metric == "instances" else "Unique 9-mers per bin"
    )
    ax.set_title(
        "Distribution of epitope sharing across alleles (Tr-lncRNA-MPs)"
        if not args.coding_control
        else "Distribution of epitope sharing across alleles"
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax.set_xticks(xs)
    ax.set_xticklabels([str(int(x)) for x in xs])
    ax.yaxis.grid(True, linestyle="-", alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(0.98, 0.96, args.dataset_label, transform=ax.transAxes, ha="right", va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)
    if not args.no_repo_mirror:
        _mirror_to_repo_figures(args.out_png, args.out_csv)
    print(
        f"Unique 9-mers (>=1 SB allele): {len(n_alleles_per_epitope)} | "
        f"alleles in panel: {len(alleles)} | max sharing: {kmax}\n"
        f"Wrote {args.out_png}\n{args.out_csv}"
    )


if __name__ == "__main__":
    main()
