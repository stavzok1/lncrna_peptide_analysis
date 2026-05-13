"""
Figure 5D / 5E: per-allele **SB epitopes** (default) or unique strong-binding 9-mers (NetMHCpan wide XLS).

**Strong binder:** predicted IC50 < ``--ic50-nm`` (default 150 nM) from ``BA_score``, same as
``plot_netmhc_epitopes_vs_hla_frequency.py``.

- **Fig 5D (default):** ``netmhcpan_sig_lnc.xls`` — significant lncRNA-derived 9-mers.
- **Fig 5E:** ``--coding-control`` → ``netmhcpan_coding_proportional_whole.xls`` — whole-protein coding
  control sample (proportional to significant MP lengths; same cohort as Fig 5C).

Bars are sorted **descending** by the chosen count metric. The title states whether the panel is
**Tr-lncRNA MPs** (5D) or the **proteome control sample** (5E). PNG and CSV are also copied to
repo-root ``figures/``.
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from plot_netmhc_epitopes_vs_hla_frequency import (
    ba_score_to_ic50_nm,
    display_allele,
    parse_wide_netmhc_xls,
)

import figure_palettes as pal

NET = ROOT / "data" / "netmhc"
FIGS = NET / "figures"
REPO_FIGURES = ROOT / "figures"


def _mirror_to_repo_figures(*paths: Path) -> None:
    REPO_FIGURES.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.is_file():
            shutil.copy2(p, REPO_FIGURES / p.name)


def _coding_sample_title_suffix(xls: Path) -> str:
    name = xls.name.lower()
    if "proportional_whole" in name:
        return "Coding control (whole proteins, proportional sample; 9-mers)"
    if "coding_control" in name:
        return "Coding control (random proteome fragments; 9-mers)"
    return "Coding control (9-mers)"
def epitope_counts_per_allele(
    alleles: list[str], strong: np.ndarray, peps: list[str]
) -> list[tuple[str, int, int]]:
    """Return (display_allele, n_unique_epitopes, n_sb_row_instances) per allele index, unsorted."""
    n_alleles = len(alleles)
    out: list[tuple[str, int, int]] = []
    for i in range(n_alleles):
        uniq = {peps[j] for j in range(len(peps)) if strong[j, i]}
        n_inst = int(strong[:, i].sum())
        out.append((display_allele(alleles[i]), len(uniq), n_inst))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 5D/5E: strong epitopes per allele (sorted bars).")
    ap.add_argument(
        "--coding-control",
        action="store_true",
        help="Use coding / proteome control XLS (Fig 5E).",
    )
    ap.add_argument("--netmhc-xls", type=Path, default=None)
    ap.add_argument("--ic50-nm", type=float, default=150.0)
    ap.add_argument(
        "--count-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Bar height: SB prediction rows vs unique 9-mers per allele (default instances).",
    )
    ap.add_argument("--out-png", type=Path, default=None)
    ap.add_argument("--out-csv", type=Path, default=None)
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
            FIGS / "fig5e_epitopes_per_allele.png"
            if args.coding_control
            else FIGS / "fig5d_epitopes_per_allele.png"
        )
    if args.out_csv is None:
        args.out_csv = (
            FIGS / "fig5e_epitopes_per_allele.csv"
            if args.coding_control
            else FIGS / "fig5d_epitopes_per_allele.csv"
        )
    panel = "E" if args.coding_control else "D"
    alleles, ba_mat, peps = parse_wide_netmhc_xls(args.netmhc_xls)
    ic50 = np.vectorize(ba_score_to_ic50_nm)(ba_mat)
    strong = ic50 < args.ic50_nm

    pairs = epitope_counts_per_allele(alleles, strong, peps)
    sort_i = 2 if args.count_metric == "instances" else 1
    pairs.sort(key=lambda x: -x[sort_i])
    labels = [p[0] for p in pairs]
    n_uniq = np.array([p[1] for p in pairs], dtype=int)
    n_inst = np.array([p[2] for p in pairs], dtype=int)
    counts = n_inst if args.count_metric == "instances" else n_uniq

    FIGS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "allele": labels,
            "n_sb_row_instances": n_inst,
            "n_unique_strong_epitopes": n_uniq,
        }
    ).to_csv(args.out_csv, index=False)

    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(10.0, 0.35 * n), 5.2), dpi=150)
    x = np.arange(n)
    ax.bar(x, counts, color=pal.BAR_SINGLE_SERIES, edgecolor="white", linewidth=0.35, width=0.88)
    ax.set_ylabel(
        "SB epitopes per allele" if args.count_metric == "instances" else "Unique 9-mers per allele"
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax.set_xlabel("Allele")
    sample_desc = (
        "Tr-lncRNA-MPs"
        if not args.coding_control
        else _coding_sample_title_suffix(args.netmhc_xls)
    )
    ax.set_title(
        f"{'SB epitopes' if args.count_metric == 'instances' else 'Unique 9-mers'} per allele — {sample_desc}"
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=7)
    ax.yaxis.grid(True, linestyle="-", alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.5, n - 0.5)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.28)
    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)
    if not args.no_repo_mirror:
        _mirror_to_repo_figures(args.out_png, args.out_csv)
    print(f"Wrote {args.out_png}\n{args.out_csv} ({n} alleles, panel {panel})")


if __name__ == "__main__":
    main()
