"""
Canonical NetMHC manuscript figures (outputs mirrored to repo ``figures/`` where applicable).

Runs **merged 5A–5C / 5D–5E** (IEDB + NetMHC SB defaults on ``*_with_iedb.tsv``) and **Figure 6**
TTN-AS1 coverage (split panels, **instances** only by default). Pass ``--also-write-unique`` to also
emit ``*_unique_*`` companion PNGs under ``figures/``. **Wide-XLS Fig 5A** (IC50-from-BA only) is
**not** part of this bundle; use ``generate_netmhc_supplement.py --include-wide-xls-fig5``.

Random-fragment merged-coding **5C** and **5D–5E** mirrors are written under
``figures/supplementary/netmhc/coding_fragments_random_sample/`` unless you pass
``--canonical-main-text-only`` (main-text bundle: proportional-whole coding cohort only, no
random-fragment mirrors).

Supplement jobs (wide 5A–5E legacy, cohort sensitivity, SB combination grid, Fig 6 SB sweeps) live
under ``supplement/`` and are run via ``generate_netmhc_fig5_fig6_supplement.py`` (or the deprecated ``generate_netmhc_supplement.py`` wrapper).

Does **not** replace ``generate_catalog_figures.py`` (Figures 2–4A).

Usage::

    python generate_netmhc_figure_bundle.py
    python generate_netmhc_figure_bundle.py --strict
    python generate_netmhc_figure_bundle.py --canonical-main-text-only
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import (
    FIGURES,
    FIGURES_SUPPLEMENTARY_NETMHC_CODING_FRAGMENTS,
    MANUSCRIPT_DIR,
    NETMHC_DATA,
    REPO_ROOT,
)

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
    ap.add_argument(
        "--canonical-main-text-only",
        action="store_true",
        help="Omit random coding-fragment NetMHC mirrors (supplementary 5C/5DE). Fig 6 remains "
        "instances-only unless you also pass --also-write-unique.",
    )
    ap.add_argument(
        "--supplement-mirrors-only",
        action="store_true",
        help="Only random-fragment NetMHC mirrors under figures/supplementary/netmhc/ "
        "(no main-text 5A–5E or Fig 6). Used by generate_supplementary_figures.py.",
    )
    ap.add_argument(
        "--also-write-unique",
        action="store_true",
        help="With Fig 6 split panels, also write *_unique_* companion PNGs next to *_instances_* "
        "(same directory as -o). Default is instances only.",
    )
    args = ap.parse_args()

    if args.supplement_mirrors_only and args.canonical_main_text_only:
        raise SystemExit("Use at most one of --canonical-main-text-only and --supplement-mirrors-only.")

    failures: list[tuple[str, int]] = []

    def step(script: str, extra: list[str]) -> None:
        code = run(script, extra)
        if code != 0:
            failures.append((script, code))

    if not args.skip_iedb_pipeline and not args.supplement_mirrors_only:
        step("plot_fig5abc_netmhc_sb_triple.py", [])
        if not args.canonical_main_text_only:
            step(
                "plot_fig5abc_netmhc_sb_triple.py",
                [
                    "--panels",
                    "c",
                    "--coding-tsv",
                    str(NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv"),
                    "--output-stem",
                    "figS5c_random_fragments",
                    "--repo-mirror-dir",
                    str(FIGURES_SUPPLEMENTARY_NETMHC_CODING_FRAGMENTS),
                ],
            )
        if not args.supplement_mirrors_only:
            step(
                "plot_fig5de_merged_iedb_sb_per_allele.py",
                [
                    "--coding-tsv",
                    str(NETMHC_DATA / "netmhcpan_coding_proportional_whole_with_iedb.tsv"),
                    "--output-stem",
                    "fig5de_merged_whole",
                ],
            )
        if not args.canonical_main_text_only:
            step(
                "plot_fig5de_merged_iedb_sb_per_allele.py",
                [
                    "--coding-tsv",
                    str(NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv"),
                    "--output-stem",
                    "figS5de_random_fragments",
                    "--repo-mirror-dir",
                    str(FIGURES_SUPPLEMENTARY_NETMHC_CODING_FRAGMENTS),
                ],
            )

    if not args.supplement_mirrors_only:
        ttn_out = FIGURES / "fig6_ttn_as1_split.png"
        ttn_args: list[str] = ["--split-panels"]
        if args.also_write_unique:
            ttn_args.append("--also-write-unique")
        ttn_args.extend(["-o", str(ttn_out)])
        step("plot_figure6_ttn_as1_allele_coverage.py", ttn_args)

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with failures:", failures)
    else:
        print("NetMHC canonical figure bundle done. See docs/figure_catalog.md, docs/netmhc_figure_commands.md, and figures/supplementary/README.md.")


if __name__ == "__main__":
    main()
