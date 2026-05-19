"""
**Figure 5–6 NetMHC supplement** — non-redundant coverage in one tree.

Writes under ``figures/supplementary/netmhc_fig5_fig6_supplement/``:

- ``fig5_merged_cohort_1d_sensitivity_loo/`` — **Fig 5** merged sig+coding cohorts: **one gate at a time**
  + leave-one-filter-out (``supplement/netmhc_sb_sensitivity_robustness.py``).
- ``fig5_merged_cohort_cartesian_sb_grid/`` — **Fig 5** merged cohorts: **Cartesian** immuno × proc ×
  EL × IC50 + sharing panels (``supplement/plot_fig5_netmhc_sb_combination_grid.py``).
- ``fig6_ttn_wide_netmhc_sb_sweeps/`` — **Fig 6** TTN wide XLS **only** (``supplement/plot_figure6_ttn_as1_sb_sensitivity.py``).
- ``fig6_ttn_merged_iedb_1d_sensitivity_loo/`` — **Fig 6** TTN XLS + IEDB merge: **1D sweeps + LOO** on the
  merged table (``supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py``), parallel in role to
  Fig 5’s ``netmhc_sb_sensitivity_robustness.py``.
- ``fig6_ttn_merged_iedb_cartesian_sb_grid/`` — **Fig 6** same merge: **Cartesian** threshold grid
  (``supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py``).

Run ``generate_netmhc_figure_bundle.py`` first so merged ``*_with_iedb.tsv`` inputs exist for the Fig 5 steps.

``generate_netmhc_supplement.py`` (legacy wrapper) may invoke this script; prefer calling it directly.

Usage::

    python generate_netmhc_fig5_fig6_supplement.py
    python generate_netmhc_fig5_fig6_supplement.py --strict
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import FIGURES_SUPPLEMENTARY_NETMHC_FIG5_FIG6, REPO_ROOT, SUPPLEMENT_DIR

SUP = SUPPLEMENT_DIR


def run_supp(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(SUP / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="Exit non-zero if any step fails.")
    args = ap.parse_args()

    root = FIGURES_SUPPLEMENTARY_NETMHC_FIG5_FIG6
    d_1d = root / "fig5_merged_cohort_1d_sensitivity_loo"
    d_f5grid = root / "fig5_merged_cohort_cartesian_sb_grid"
    d_ttn_nm = root / "fig6_ttn_wide_netmhc_sb_sweeps"
    d_ttn_1d = root / "fig6_ttn_merged_iedb_1d_sensitivity_loo"
    d_ttn_ie = root / "fig6_ttn_merged_iedb_cartesian_sb_grid"
    for p in (d_1d, d_f5grid, d_ttn_nm, d_ttn_1d, d_ttn_ie):
        p.mkdir(parents=True, exist_ok=True)

    failures: list[tuple[str, int]] = []

    def step(script: str, extra: list[str]) -> None:
        if run_supp(script, extra) != 0:
            failures.append((script, 1))

    step("netmhc_sb_sensitivity_robustness.py", ["--out-dir", str(d_1d)])
    step("plot_fig5_netmhc_sb_combination_grid.py", ["--out-dir", str(d_f5grid)])
    step("plot_figure6_ttn_as1_sb_sensitivity.py", ["--out-dir", str(d_ttn_nm)])
    step("netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py", ["--out-dir", str(d_ttn_1d)])
    step("plot_fig6_ttn_merged_iedb_sb_combination_grid.py", ["--out-dir", str(d_ttn_ie)])

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with failures:", failures)
    else:
        print(f"Fig 5–6 NetMHC supplement figures done under {root}")


if __name__ == "__main__":
    main()
