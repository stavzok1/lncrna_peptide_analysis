"""
NetMHC **supplement** figures: legacy wide-XLS cohort 5B–5E (IC50-from-BA), and — unless skipped —
the **Fig 5–6 supplement bundle** (five subfolders: Fig 5 1D+LOO, Fig 5 Cartesian, Fig 6 NetMHC sweeps, Fig 6 merged IEDB 1D+LOO, Fig 6 merged IEDB Cartesian) under ``figures/supplementary/netmhc_fig5_fig6_supplement/``
(via ``generate_netmhc_fig5_fig6_supplement.py``).

Run ``generate_netmhc_figure_bundle.py`` first for canonical merged panels.

Usage::

    python generate_netmhc_supplement.py
    python generate_netmhc_supplement.py --include-wide-xls-fig5 --strict
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import NETMHC_DATA, NETMHC_FIGURES, NETMHC_HLA27_ALLELE_FREQ_CSV, REPO_ROOT, SUPPLEMENT_DIR

SUP = SUPPLEMENT_DIR
MS = REPO_ROOT / "manuscript"


def run(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(SUP / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_manuscript(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(MS / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_root(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(REPO_ROOT / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="Exit non-zero if any step fails.")
    ap.add_argument(
        "--include-wide-xls-fig5",
        action="store_true",
        help="Run legacy wide-XLS Fig 5A–5E (IC50-from-BA; requires cohort ``*.xls`` where referenced).",
    )
    ap.add_argument(
        "--skip-sensitivity",
        action="store_true",
        help="Skip generate_netmhc_fig5_fig6_supplement.py (cohort sensitivity, Fig 5 combination grid, Fig 6 TTN SB sweeps).",
    )
    args = ap.parse_args()
    failures: list[tuple[str, int]] = []

    def step(script: str, extra: list[str]) -> None:
        if run(script, extra) != 0:
            failures.append((script, 1))

    def step_ms(script: str, extra: list[str]) -> None:
        if run_manuscript(script, extra) != 0:
            failures.append((script, 1))

    if args.include_wide_xls_fig5:
        extra_5a: list[str] = []
        if NETMHC_HLA27_ALLELE_FREQ_CSV.is_file():
            extra_5a.extend(["--freq-file", str(NETMHC_HLA27_ALLELE_FREQ_CSV)])
        else:
            print(
                "Warning: missing bundled HLA frequency CSV; wide 5A relies on script defaults / existing sidecars.",
                flush=True,
            )
        step_ms("plot_netmhc_epitopes_vs_hla_frequency.py", extra_5a)
        step("plot_figure5b_epitope_sharing_across_alleles.py", [])
        step("plot_figure5b_epitope_sharing_across_alleles.py", ["--coding-control"])
        step(
            "plot_figure5b_epitope_sharing_across_alleles.py",
            [
                "--coding-control",
                "--netmhc-xls",
                str(NETMHC_DATA / "netmhcpan_coding_control.xls"),
                "--out-png",
                str(NETMHC_FIGURES / "fig5c_epitope_sharing_fragments_control.png"),
                "--out-csv",
                str(NETMHC_FIGURES / "fig5c_epitope_sharing_fragments_control.csv"),
            ],
        )
        step("plot_figure5de_epitopes_per_allele.py", [])
        step("plot_figure5de_epitopes_per_allele.py", ["--coding-control"])
        step(
            "plot_figure5de_epitopes_per_allele.py",
            [
                "--coding-control",
                "--netmhc-xls",
                str(NETMHC_DATA / "netmhcpan_coding_control.xls"),
                "--out-png",
                str(NETMHC_FIGURES / "fig5e_epitopes_per_allele_fragments_control.png"),
                "--out-csv",
                str(NETMHC_FIGURES / "fig5e_epitopes_per_allele_fragments_control.csv"),
            ],
        )

    if not args.skip_sensitivity:
        f56_extra: list[str] = []
        if args.strict:
            f56_extra.append("--strict")
        if run_root("generate_netmhc_fig5_fig6_supplement.py", f56_extra) != 0:
            failures.append(("generate_netmhc_fig5_fig6_supplement.py", 1))

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with failures:", failures)
    else:
        print("NetMHC supplement bundle done.")


if __name__ == "__main__":
    main()
