#!/usr/bin/env python3
"""
Rebuild a **clean** NetMHC manuscript figure tree under ``figures/manuscript_netmhc/``.

Layout (each figure number has its own subtree; SB variants and instance/unique splits below):

- ``fig5_wide_ic50_lt_150nm/{instances,unique}/`` — wide XLS 5A–5E (IC50-from-BA gate only).
- ``fig5_merged/sb_{full,no_ic50,ic50_only}/{instances,unique}/`` — merged ``*_with_iedb.tsv`` 5A–5C.
- ``fig5_merged_de/sb_{...}/{instances,unique}/`` — merged 5D–5E per-allele bars (**two** stems per folder: ``*_proportional_whole`` + ``*_random_fragments`` coding cohorts)..
- ``fig6/netmhc_default/`` — NetMHC-only gating; split panels + instances + unique companion.
- ``fig6/iedb_sb_{full,no_ic50,ic50_only}/`` — ``--gating iedb_sb`` with matching ``--sb-mode``.
- ``sensitivity/cohort_sb_{...}_{instances,unique}/`` — ``netmhc_sb_sensitivity_robustness.py``.
- ``sensitivity/fig6_netmhc_sb_sweeps/`` — ``plot_figure6_ttn_as1_sb_sensitivity.py``.

Optional purges (see ``--help``) remove legacy flat exports under ``figures/`` (top-level ``fig5*`` / ``fig6*``),
files under ``figures/supplementary/netmhc/``, and/or PNG+CSV artifacts under ``data/netmhc/figures/`` while **keeping**
the allele-frequency CSVs used as ``--freq-file`` for wide 5A.

Run from repository root::

    python supplement/regenerate_manuscript_netmhc_figures.py --clean --purge-repo-figures-netmhc \\
        --purge-data-netmhc-figures
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import (
    REPO_ROOT,
    DATA,
    FIGURES,
    FIGURES_SUPPLEMENTARY_NETMHC,
    NETMHC_DATA,
    NETMHC_FIGURES,
    NETMHC_HLA27_ALLELE_FREQ_CSV,
    MANUSCRIPT_DIR,
    SUPPLEMENT_DIR,
)

ROOT = REPO_ROOT
MANUSCRIPT = MANUSCRIPT_DIR
SUPPLEMENT = SUPPLEMENT_DIR

import argparse
import shutil
import subprocess

SB_MODES = ("full", "no_ic50", "ic50_only")
METRICS = ("instances", "unique")

FREQ_CSV_KEEP = frozenset(
    {
        "fig5a_epitopes_vs_allele_frequency_ic50_sb.csv",
        "epitopes_vs_allele_frequency_ic50_sb.csv",
        "hla_european27_allele_frequencies.csv",
    }
)

TTN_IEDB_CSV = NETMHC_DATA / "ttn_as1_iedb_companion_synthetic.csv"
TTN_IEDB_PARENT_ID = "108065|TTN-AS1|synthetic"


def resolve_freq_file() -> Path:
    if NETMHC_HLA27_ALLELE_FREQ_CSV.is_file():
        return NETMHC_HLA27_ALLELE_FREQ_CSV
    for name in (
        "fig5a_epitopes_vs_allele_frequency_ic50_sb.csv",
        "epitopes_vs_allele_frequency_ic50_sb.csv",
    ):
        p = NETMHC_FIGURES / name
        if p.is_file():
            return p
    raise FileNotFoundError(
        "No allele frequency table found. Need one of:\n"
        f"  {NETMHC_HLA27_ALLELE_FREQ_CSV}\n"
        f"  {NETMHC_FIGURES / 'fig5a_epitopes_vs_allele_frequency_ic50_sb.csv'}\n"
        f"  {NETMHC_FIGURES / 'epitopes_vs_allele_frequency_ic50_sb.csv'}"
    )


def run(cmd: list[str], *, dry_run: bool) -> None:
    print("+", " ".join(cmd), flush=True)
    if dry_run:
        return
    r = subprocess.run(cmd, cwd=str(ROOT))
    if r.returncode != 0:
        raise SystemExit(f"Command failed ({r.returncode}): {' '.join(cmd)}")


def purge_repo_netmhc_flat_figures() -> None:
    fig = FIGURES
    if not fig.is_dir():
        return
    for pat in ("fig5*", "fig6*"):
        for p in fig.glob(pat):
            if p.is_file():
                p.unlink()
                print(f"Removed {p.relative_to(ROOT)}")
    supp_nm = FIGURES_SUPPLEMENTARY_NETMHC
    if supp_nm.is_dir():
        for p in sorted(supp_nm.rglob("*"), reverse=True):
            if p.is_file():
                p.unlink()
                print(f"Removed {p.relative_to(ROOT)}")
            elif p.is_dir():
                try:
                    p.rmdir()
                except OSError:
                    pass


def purge_data_netmhc_figure_artifacts() -> None:
    base = NETMHC_FIGURES
    if not base.is_dir():
        return
    for path in sorted(base.rglob("*")):
        if not path.is_file():
            continue
        suf = path.suffix.lower()
        if suf not in (".png", ".csv"):
            continue
        if suf == ".csv" and path.name in FREQ_CSV_KEEP:
            continue
        path.unlink()
        print(f"Removed {path.relative_to(ROOT)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--out-root",
        type=Path,
        default=ROOT / "figures" / "manuscript_netmhc",
        help="All regenerated products go under this directory.",
    )
    ap.add_argument(
        "--clean",
        action="store_true",
        help="Delete --out-root recursively before regenerating.",
    )
    ap.add_argument(
        "--purge-repo-figures-netmhc",
        action="store_true",
        help="Remove loose fig5*/fig6* files directly under repo figures/ (not under subfolders).",
    )
    ap.add_argument(
        "--purge-data-netmhc-figures",
        action="store_true",
        help="Delete *.png and *.csv under data/netmhc/figures/ recursively, except freq CSV allowlist.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only.",
    )
    ap.add_argument(
        "--skip-wide",
        action="store_true",
        help="Skip wide-XLS Figure 5A–5E (requires --freq-file inputs).",
    )
    ap.add_argument(
        "--skip-merged",
        action="store_true",
        help="Skip merged-TSV Figure 5ABC / 5DE.",
    )
    ap.add_argument(
        "--skip-fig6",
        action="store_true",
        help="Skip TTN Figure 6 panels.",
    )
    ap.add_argument(
        "--skip-sensitivity",
        action="store_true",
        help="Skip cohort + TTN sensitivity scripts.",
    )
    ap.add_argument(
        "--with-combination-grid",
        action="store_true",
        help="Also run plot_fig5_netmhc_sb_combination_grid.py (slow supplement).",
    )
    args = ap.parse_args()
    dry = bool(args.dry_run)
    out_root: Path = args.out_root
    py = sys.executable

    if args.purge_repo_figures_netmhc and not dry:
        purge_repo_netmhc_flat_figures()
    if args.purge_data_netmhc_figures and not dry:
        purge_data_netmhc_figure_artifacts()

    if args.clean and not dry:
        if out_root.exists():
            shutil.rmtree(out_root)
    if not dry:
        out_root.mkdir(parents=True, exist_ok=True)

    freq = resolve_freq_file()
    print(f"Using allele freq table: {freq.relative_to(ROOT)}", flush=True)

    # --- Wide Fig 5 (IC50-from-BA only) ---
    if not args.skip_wide:
        base = out_root / "fig5_wide_ic50_lt_150nm"
        for metric in METRICS:
            d = base / metric
            if not dry:
                d.mkdir(parents=True, exist_ok=True)
            stem = f"fig5a_wide_{metric}"
            run(
                [
                    py,
                    str(MANUSCRIPT / "plot_netmhc_epitopes_vs_hla_frequency.py"),
                    "--freq-file",
                    str(freq),
                    "--y-metric",
                    metric,
                    "--out-png",
                    str(d / f"{stem}.png"),
                    "--out-csv",
                    str(d / f"{stem}.csv"),
                    "--no-repo-mirror",
                ],
                dry_run=dry,
            )
            run(
                [
                    py,
                    str(SUPPLEMENT / "plot_figure5b_epitope_sharing_across_alleles.py"),
                    "--histogram-metric",
                    metric,
                    "--out-png",
                    str(d / f"fig5b_wide_{metric}.png"),
                    "--out-csv",
                    str(d / f"fig5b_wide_{metric}.csv"),
                    "--no-repo-mirror",
                ],
                dry_run=dry,
            )
            run(
                [
                    py,
                    str(SUPPLEMENT / "plot_figure5b_epitope_sharing_across_alleles.py"),
                    "--coding-control",
                    "--histogram-metric",
                    metric,
                    "--out-png",
                    str(d / f"fig5c_wide_{metric}.png"),
                    "--out-csv",
                    str(d / f"fig5c_wide_{metric}.csv"),
                    "--no-repo-mirror",
                ],
                dry_run=dry,
            )
            run(
                [
                    py,
                    str(SUPPLEMENT / "plot_figure5de_epitopes_per_allele.py"),
                    "--count-metric",
                    metric,
                    "--out-png",
                    str(d / f"fig5d_wide_{metric}.png"),
                    "--out-csv",
                    str(d / f"fig5d_wide_{metric}.csv"),
                    "--no-repo-mirror",
                ],
                dry_run=dry,
            )
            run(
                [
                    py,
                    str(SUPPLEMENT / "plot_figure5de_epitopes_per_allele.py"),
                    "--coding-control",
                    "--count-metric",
                    metric,
                    "--out-png",
                    str(d / f"fig5e_wide_{metric}.png"),
                    "--out-csv",
                    str(d / f"fig5e_wide_{metric}.csv"),
                    "--no-repo-mirror",
                ],
                dry_run=dry,
            )

    # --- Merged Fig 5 ABC / DE (IEDB + NetMHC SB modes) ---
    if not args.skip_merged:
        for sb in SB_MODES:
            for metric in METRICS:
                d_abc = out_root / "fig5_merged" / f"sb_{sb}" / metric
                d_de = out_root / "fig5_merged_de" / f"sb_{sb}" / metric
                if not dry:
                    d_abc.mkdir(parents=True, exist_ok=True)
                    d_de.mkdir(parents=True, exist_ok=True)
                stem_abc = f"fig5abc_sb_{sb}_{metric}"
                stem_de = f"fig5de_sb_{sb}_{metric}"
                run(
                    [
                        py,
                        str(MANUSCRIPT / "plot_fig5abc_netmhc_sb_triple.py"),
                        "--allele-freq-csv",
                        str(freq),
                        "--sb-mode",
                        sb,
                        "--out-dir",
                        str(d_abc),
                        "--output-stem",
                        stem_abc,
                        "--fig5a-y-metric",
                        metric,
                        "--sharing-y-metric",
                        metric,
                        "--no-repo-mirror",
                    ],
                    dry_run=dry,
                )
                run(
                    [
                        py,
                        str(MANUSCRIPT / "plot_fig5abc_netmhc_sb_triple.py"),
                        "--coding-tsv",
                        str(NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv"),
                        "--sb-mode",
                        sb,
                        "--out-dir",
                        str(d_abc),
                        "--output-stem",
                        f"{stem_abc}_coding_control",
                        "--sharing-y-metric",
                        metric,
                        "--panels",
                        "c",
                        "--no-repo-mirror",
                    ],
                    dry_run=dry,
                )
                run(
                    [
                        py,
                        str(MANUSCRIPT / "plot_fig5de_merged_iedb_sb_per_allele.py"),
                        "--sb-mode",
                        sb,
                        "--out-dir",
                        str(d_de),
                        "--coding-tsv",
                        str(NETMHC_DATA / "netmhcpan_coding_proportional_whole_with_iedb.tsv"),
                        "--output-stem",
                        f"{stem_de}_proportional_whole",
                        "--count-metric",
                        metric,
                        "--no-repo-mirror",
                    ],
                    dry_run=dry,
                )
                run(
                    [
                        py,
                        str(MANUSCRIPT / "plot_fig5de_merged_iedb_sb_per_allele.py"),
                        "--sb-mode",
                        sb,
                        "--out-dir",
                        str(d_de),
                        "--coding-tsv",
                        str(NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv"),
                        "--output-stem",
                        f"{stem_de}_random_fragments",
                        "--count-metric",
                        metric,
                        "--no-repo-mirror",
                    ],
                    dry_run=dry,
                )

    # --- Fig 6 TTN ---
    if not args.skip_fig6:
        d6n = out_root / "fig6" / "netmhc_default"
        if not dry:
            d6n.mkdir(parents=True, exist_ok=True)
        run(
            [
                py,
                str(MANUSCRIPT / "plot_figure6_ttn_as1_allele_coverage.py"),
                "--gating",
                "netmhc",
                "--split-panels",
                "--also-write-unique",
                "-o",
                str(d6n / "fig6_ttn_split.png"),
            ],
            dry_run=dry,
        )

        if not TTN_IEDB_CSV.is_file():
            print(f"Skip IEDB-gated Fig6: missing {TTN_IEDB_CSV}", flush=True)
        else:
            for sb in SB_MODES:
                d6i = out_root / "fig6" / f"iedb_sb_{sb}"
                if not dry:
                    d6i.mkdir(parents=True, exist_ok=True)
                run(
                    [
                        py,
                        str(MANUSCRIPT / "plot_figure6_ttn_as1_allele_coverage.py"),
                        "--gating",
                        "iedb_sb",
                        "--sb-mode",
                        sb,
                        "--iedb-csv",
                        str(TTN_IEDB_CSV),
                        "--iedb-parent-input-seq-id",
                        TTN_IEDB_PARENT_ID,
                        "--split-panels",
                        "--also-write-unique",
                        "-o",
                        str(d6i / "fig6_ttn_split.png"),
                    ],
                    dry_run=dry,
                )

    # --- Sensitivity ---
    if not args.skip_sensitivity:
        for sb in SB_MODES:
            for metric in METRICS:
                sd = out_root / "sensitivity" / "cohort" / f"sb_{sb}_{metric}"
                if not dry:
                    sd.mkdir(parents=True, exist_ok=True)
                run(
                    [
                        py,
                        str(SUPPLEMENT / "netmhc_sb_sensitivity_robustness.py"),
                        "--sb-mode",
                        sb,
                        "--plot-metric",
                        metric,
                        "--out-dir",
                        str(sd),
                        "--out-stem",
                        f"sb_sensitivity_{sb}_{metric}",
                    ],
                    dry_run=dry,
                )
        s6 = out_root / "sensitivity" / "fig6_netmhc_sb_sweeps"
        if not dry:
            s6.mkdir(parents=True, exist_ok=True)
        run(
            [
                py,
                str(SUPPLEMENT / "plot_figure6_ttn_as1_sb_sensitivity.py"),
                "--out-dir",
                str(s6),
            ],
            dry_run=dry,
        )

    if args.with_combination_grid:
        cg = out_root / "supplement" / "fig5_sb_combination_grid"
        if not dry:
            cg.mkdir(parents=True, exist_ok=True)
        run(
            [
                py,
                str(SUPPLEMENT / "plot_fig5_netmhc_sb_combination_grid.py"),
                "--out-dir",
                str(cg),
            ],
            dry_run=dry,
        )

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
