"""
Canonical NetMHC manuscript figures (outputs mirrored to repo ``figures/`` where applicable).

Runs **merged 5A–5C / 5D–5E** (IEDB + NetMHC SB defaults on ``*_with_iedb.tsv``) and **Figure 6**
TTN-AS1 coverage (split panels + unique companion). **Wide-XLS Fig 5A** (IC50-from-BA only) is
**not** part of this bundle; use ``generate_netmhc_supplement.py --include-wide-xls-fig5``.

Supplement jobs (wide 5A–5E legacy, cohort sensitivity, SB combination grid, Fig 6 SB sweeps) live
under ``supplement/`` and are run via ``generate_netmhc_supplement.py``.

Does **not** replace ``generate_catalog_figures.py`` (Figures 2–4A).

Usage::

    python generate_netmhc_figure_bundle.py
    python generate_netmhc_figure_bundle.py --strict
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import FIGURES, MANUSCRIPT_DIR, NETMHC_DATA, REPO_ROOT

MS = MANUSCRIPT_DIR


def run(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(MS / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero if any step fails.",
    )
    ap.add_argument(
        "--skip-iedb-pipeline",
        action="store_true",
        help="Skip scripts that require ``*_with_iedb.tsv`` merges.",
    )
    args = ap.parse_args()

    failures: list[tuple[str, int]] = []

    def step(script: str, extra: list[str]) -> None:
        code = run(script, extra)
        if code != 0:
            failures.append((script, code))

    if not args.skip_iedb_pipeline:
        step("plot_fig5abc_netmhc_sb_triple.py", [])
        step(
            "plot_fig5abc_netmhc_sb_triple.py",
            [
                "--panels",
                "c",
                "--coding-tsv",
                str(NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv"),
                "--output-stem",
                "fig5abc_sb_immuno_proc_el_ic50_coding_control",
            ],
        )
        step(
            "plot_fig5de_merged_iedb_sb_per_allele.py",
            [
                "--coding-tsv",
                str(NETMHC_DATA / "netmhcpan_coding_proportional_whole_with_iedb.tsv"),
                "--output-stem",
                "fig5de_merged_iedb_sb_proportional_whole",
            ],
        )
        step(
            "plot_fig5de_merged_iedb_sb_per_allele.py",
            [
                "--coding-tsv",
                str(NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv"),
                "--output-stem",
                "fig5de_merged_iedb_sb_random_fragments",
            ],
        )

    ttn_out = FIGURES / "fig6_ttn_as1_split.png"
    step(
        "plot_figure6_ttn_as1_allele_coverage.py",
        [
            "--split-panels",
            "--also-write-unique",
            "-o",
            str(ttn_out),
        ],
    )

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with failures:", failures)
    else:
        print("NetMHC canonical figure bundle done. See docs/figure_catalog.md and docs/netmhc_figure_commands.md.")


if __name__ == "__main__":
    main()
