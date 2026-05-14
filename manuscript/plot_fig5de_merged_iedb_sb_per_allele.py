#!/usr/bin/env python3
"""
Figure **5D / 5E**-style per-allele bar charts using **merged** ``*_with_iedb.tsv`` cohort tables
(same SB definitions as ``plot_fig5abc_netmhc_sb_triple.py``).

Use this when you want **5D/5E layout** (sorted allele bars) under different ``--sb-mode`` /
immuno / processing / EL / IC50 settings without re-running the wide NetMHC ``*.xls`` pipeline.

Examples::

    # Full IEDB+EL+IC50 stack (defaults from netmhc_sb_core)
    python scripts/plot_fig5de_merged_iedb_sb_per_allele.py --out-dir data/netmhc/figures

    # Second coding cohort (distinct 5E/5D filenames; 5D is sig-lnc, duplicated across runs):
    python scripts/plot_fig5de_merged_iedb_sb_per_allele.py \\
        --coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv \\
        --output-stem fig5de_merged_iedb_sb_random_fragments

    # IC50 binding only (no IEDB immuno/proc/EL gates)
    python scripts/plot_fig5de_merged_iedb_sb_per_allele.py --sb-mode ic50_only --output-stem fig5de_merged_ic50_only

    # Write all three SB modes into subfolders
    python scripts/plot_fig5de_merged_iedb_sb_per_allele.py --write-all-sb-modes
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES
from figure_export import add_publication_args, save_figure_bundle

ROOT = REPO_ROOT


import argparse
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

import figure_palettes as pal  # noqa: E402

from netmhc_sb_core import (  # noqa: E402
    FIG5_IEDB_EL_RANK_MAX_DEFAULT,
    FIG5_IEDB_IC50_MAX_NM_DEFAULT,
    FIG5_IEDB_IMM_MIN_DEFAULT,
    FIG5_IEDB_PROC_MIN_DEFAULT,
    ba_score_min_for_ic50_lt,
    pick_iedb_ic50_column,
    sb_mask_fig5_defaults,
)

REPO_FIGURES = FIGURES


def _mirror(*paths: Path) -> None:
    REPO_FIGURES.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.is_file():
            shutil.copy2(p, REPO_FIGURES / p.name)


def _norm_display_allele(name: str) -> str:
    """HLA-A01:01 -> HLA-A*01:01 for tick labels (matches wide Fig5 scripts)."""
    import re

    u = str(name).strip().upper().replace("*", "")
    m = re.fullmatch(r"HLA-([AB])(\d{2}:\d{2})", u)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return str(name)


def per_allele_table(sb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for al, g in sb.groupby("allele", sort=False):
        rows.append(
            {
                "allele": _norm_display_allele(al),
                "n_unique_epitopes_sb": int(g["Peptide"].nunique()),
                "n_sb_row_instances_allele": int(len(g)),
            }
        )
    return pd.DataFrame(rows)


def plot_bars(
    tab: pd.DataFrame,
    out_png: Path,
    *,
    title: str,
    subtitle: str,
    count_metric: str,
    publication_dir: Path | None = None,
    publication_tiff_kind: str = "color",
    figures_root: Path = REPO_FIGURES,
) -> None:
    ycol = "n_sb_row_instances_allele" if count_metric == "instances" else "n_unique_epitopes_sb"
    ylab = "SB epitopes per allele" if count_metric == "instances" else "Unique 9-mers per allele"
    t = tab.sort_values(ycol, ascending=False).reset_index(drop=True)
    labels = t["allele"].tolist()
    counts = t[ycol].to_numpy(dtype=int)
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(10.0, 0.35 * n), 5.2), dpi=150)
    x = np.arange(n)
    ax.bar(x, counts, color=pal.BAR_SINGLE_SERIES, edgecolor="white", linewidth=0.35, width=0.88)
    ax.set_ylabel(ylab)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax.set_xlabel("Allele")
    ax.set_title(f"{title}\n{subtitle}" if subtitle else title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=7)
    ax.yaxis.grid(True, linestyle="-", alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.28)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    save_figure_bundle(
        fig,
        out_png,
        png_dpi=150,
        publication_dir=publication_dir,
        publication_tiff_kind=publication_tiff_kind,
        figures_root=figures_root,
        bbox_inches="tight",
    )
    plt.close(fig)


def run_one(
    *,
    args: argparse.Namespace,
    sb_mode: str,
    out_dir: Path,
    stem: str,
) -> None:
    hdr_sig = set(pd.read_csv(args.sig_tsv, sep="\t", nrows=0).columns)
    hdr_cod = set(pd.read_csv(args.coding_tsv, sep="\t", nrows=0).columns)
    iedb_ic50_col = pick_iedb_ic50_column(sorted(hdr_sig | hdr_cod))
    usecols = {"allele", "Peptide", "EL_rank", "iedb_score", "iedb_processing_score"}
    if iedb_ic50_col:
        usecols.add(iedb_ic50_col)
    else:
        usecols.add("BA_score")

    ba_min = ba_score_min_for_ic50_lt(args.ic50_max_nm)
    sig = pd.read_csv(args.sig_tsv, sep="\t", usecols=lambda c: c in usecols)
    cod = pd.read_csv(args.coding_tsv, sep="\t", usecols=lambda c: c in usecols)

    sig_sb = sig.loc[
        sb_mask_fig5_defaults(
            sig,
            el_max=args.el_rank_max,
            el_lte=args.el_rank_lte,
            ba_min=ba_min,
            ic50_max_nm=args.ic50_max_nm,
            iedb_ic50_col=iedb_ic50_col,
            imm_min=args.imm_min,
            proc_min=args.proc_min,
            sb_mode=sb_mode,
        )
    ]
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
            sb_mode=sb_mode,
        )
    ]

    tab_s = per_allele_table(sig_sb)
    tab_c = per_allele_table(cod_sb)
    out_dir.mkdir(parents=True, exist_ok=True)
    tab_s.to_csv(out_dir / f"{stem}_5d_sig_per_allele.csv", index=False)
    tab_c.to_csv(out_dir / f"{stem}_5e_coding_per_allele.csv", index=False)

    plot_bars(
        tab_s,
        out_dir / f"{stem}_5d_sig_per_allele.png",
        title="Per-allele SB load (Tr-lncRNA-MPs, merged IEDB+NetMHC)",
        subtitle="",
        count_metric=args.count_metric,
        publication_dir=args.publication_dir,
        publication_tiff_kind=args.publication_tiff_kind,
        figures_root=REPO_FIGURES,
    )
    plot_bars(
        tab_c,
        out_dir / f"{stem}_5e_coding_per_allele.png",
        title="Per-allele SB load (coding control, merged IEDB+NetMHC)",
        subtitle="",
        count_metric=args.count_metric,
        publication_dir=args.publication_dir,
        publication_tiff_kind=args.publication_tiff_kind,
        figures_root=REPO_FIGURES,
    )
    if not getattr(args, "no_repo_mirror", False):
        _mirror(
            out_dir / f"{stem}_5d_sig_per_allele.png",
            out_dir / f"{stem}_5e_coding_per_allele.png",
        )
    print(f"[{sb_mode}] Wrote {stem}_5d_* and {stem}_5e_* under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sig-tsv", type=Path, default=ROOT / "data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv")
    ap.add_argument(
        "--coding-tsv",
        type=Path,
        default=ROOT / "data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv",
    )
    ap.add_argument("--out-dir", type=Path, default=ROOT / "data/netmhc/figures")
    ap.add_argument(
        "--output-stem",
        type=str,
        default="fig5de_merged_iedb_sb_proportional_whole",
        help="Filename prefix; use e.g. fig5de_merged_iedb_sb_random_fragments with the fragment merged TSV.",
    )
    ap.add_argument(
        "--sb-mode",
        choices=("full", "no_ic50", "ic50_only"),
        default="full",
        help="SB mask mode (same as plot_fig5abc_netmhc_sb_triple.py).",
    )
    ap.add_argument(
        "--write-all-sb-modes",
        action="store_true",
        help="Also write no_ic50 and ic50_only variants into subfolders sb_full, sb_no_ic50, sb_ic50_only.",
    )
    ap.add_argument("--imm-min", type=float, default=FIG5_IEDB_IMM_MIN_DEFAULT)
    ap.add_argument("--proc-min", type=float, default=FIG5_IEDB_PROC_MIN_DEFAULT)
    ap.add_argument("--el-rank-max", type=float, default=FIG5_IEDB_EL_RANK_MAX_DEFAULT)
    ap.add_argument("--el-rank-lte", action="store_true")
    ap.add_argument("--ic50-max-nm", type=float, default=FIG5_IEDB_IC50_MAX_NM_DEFAULT)
    ap.add_argument(
        "--count-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Bar height: SB epitopes vs unique 9-mers per allele.",
    )
    ap.add_argument(
        "--no-repo-mirror",
        action="store_true",
        help="Do not copy 5D/5E PNGs into repo-root figures/.",
    )
    add_publication_args(ap)
    args = ap.parse_args()

    if args.write_all_sb_modes:
        for mode, sub in (("full", "sb_full"), ("no_ic50", "sb_no_ic50"), ("ic50_only", "sb_ic50_only")):
            od = args.out_dir / sub
            stem = f"{args.output_stem}_{mode}"
            run_one(args=args, sb_mode=mode, out_dir=od, stem=stem)
    else:
        run_one(args=args, sb_mode=args.sb_mode, out_dir=args.out_dir, stem=args.output_stem)


if __name__ == "__main__":
    main()
