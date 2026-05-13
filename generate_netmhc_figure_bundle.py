"""
Orchestrate **catalog Figure 5** (NetMHCpan cohort: sig lnc vs coding) and **Figure 6**
(TTN-AS1 smPEP 108065) scripts from the repository root.

Does **not** replace ``generate_catalog_figures.py`` (Figures 2–4A); run that separately for
SmProt / TCGA composition figures.

Usage::

    python generate_netmhc_figure_bundle.py
    python generate_netmhc_figure_bundle.py --strict

Prerequisites: merged TSVs for IEDB+NetMHC pipeline (see ``data/netmhc/README_netmhc.md``),
allele-frequency CSV for merged 5A (and wide 5A). Optional: wide XLS cohort **5B–5E**
(IC50-from-BA) only when ``--include-wide-xls-fig5`` is set. **Default manuscript 5D–5E**
are merged ``*_with_iedb.tsv`` + ``sb_mode=full`` (two stems: proportional-whole vs random-fragment coding).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(ROOT / script), *args]
    print("+", " ".join(cmd))
    return subprocess.call(cmd, cwd=str(ROOT))


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
        "--include-wide-xls-fig5",
        action="store_true",
        help=(
            "Also run wide-XLS cohort scripts: ``plot_figure5b_epitope_sharing_across_alleles.py`` "
            "(5B/5C; IC50-from-BA), ``plot_figure5de_epitopes_per_allele.py`` (5D/5E main + "
            "``--coding-control`` + random-fragment XLS 5E). Default bundle uses merged "
            "``*_with_iedb.tsv`` strict SB for cohort **5B–5E** (5D/5E never from wide XLS unless this flag is set)."
        ),
    )
    args = ap.parse_args()

    failures: list[tuple[str, int]] = []

    def step(script: str, extra: list[str]) -> None:
        code = run(script, extra)
        if code != 0:
            failures.append((script, code))

    # --- Figure 5 (wide XLS cohort): 5A always; 5B–5E wide only with legacy flag ---
    step("plot_netmhc_epitopes_vs_hla_frequency.py", [])
    if args.include_wide_xls_fig5:
        step("plot_figure5b_epitope_sharing_across_alleles.py", [])
        step("plot_figure5b_epitope_sharing_across_alleles.py", ["--coding-control"])
        # Legacy random-fragment coding control (wide XLS); supplement only.
        step(
            "plot_figure5b_epitope_sharing_across_alleles.py",
            [
                "--coding-control",
                "--netmhc-xls",
                str(ROOT / "data" / "netmhc" / "netmhcpan_coding_control.xls"),
                "--out-png",
                str(ROOT / "data" / "netmhc" / "figures" / "fig5c_epitope_sharing_fragments_control.png"),
                "--out-csv",
                str(ROOT / "data" / "netmhc" / "figures" / "fig5c_epitope_sharing_fragments_control.csv"),
            ],
        )
        step("plot_figure5de_epitopes_per_allele.py", [])
        step("plot_figure5de_epitopes_per_allele.py", ["--coding-control"])
        step(
            "plot_figure5de_epitopes_per_allele.py",
            [
                "--coding-control",
                "--netmhc-xls",
                str(ROOT / "data" / "netmhc" / "netmhcpan_coding_control.xls"),
                "--out-png",
                str(ROOT / "data" / "netmhc" / "figures" / "fig5e_epitopes_per_allele_fragments_control.png"),
                "--out-csv",
                str(ROOT / "data" / "netmhc" / "figures" / "fig5e_epitopes_per_allele_fragments_control.csv"),
            ],
        )

    if not args.skip_iedb_pipeline:
        # Canonical merged 5A–5C (sig lnc + proportional-whole coding); SB = IEDB + NetMHC defaults, sb_mode=full.
        step("scripts/plot_fig5abc_netmhc_sb_triple.py", [])
        # Second merged coding cohort: 5C only (length-matched control TSV), distinct stem.
        step(
            "scripts/plot_fig5abc_netmhc_sb_triple.py",
            [
                "--panels",
                "c",
                "--coding-tsv",
                str(ROOT / "data" / "netmhc" / "netmhcpan_coding_control_with_iedb.tsv"),
                "--output-stem",
                "fig5abc_sb_immuno_proc_el_ic50_coding_control",
            ],
        )
        # Merged 5D–5E (IEDB + NetMHC SB default): sig lnc + two coding cohorts (distinct stems / filenames).
        step(
            "scripts/plot_fig5de_merged_iedb_sb_per_allele.py",
            [
                "--coding-tsv",
                str(ROOT / "data" / "netmhc" / "netmhcpan_coding_proportional_whole_with_iedb.tsv"),
                "--output-stem",
                "fig5de_merged_iedb_sb_proportional_whole",
            ],
        )
        step(
            "scripts/plot_fig5de_merged_iedb_sb_per_allele.py",
            [
                "--coding-tsv",
                str(ROOT / "data" / "netmhc" / "netmhcpan_coding_control_with_iedb.tsv"),
                "--output-stem",
                "fig5de_merged_iedb_sb_random_fragments",
            ],
        )
        step("scripts/netmhc_sb_sensitivity_robustness.py", [])
        step("scripts/plot_fig5_netmhc_sb_combination_grid.py", [])

    # --- Figure 6 (TTN-AS1 single peptide) ---
    ttn_out = ROOT / "figures" / "fig6_ttn_as1_split.png"
    step(
        "plot_figure6_ttn_as1_allele_coverage.py",
        [
            "--split-panels",
            "--also-write-unique",
            "-o",
            str(ttn_out),
        ],
    )
    step("plot_figure6_ttn_as1_sb_sensitivity.py", [])

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with failures:", failures)
    else:
        print("NetMHC figure bundle done. See docs/figure_catalog.md and docs/netmhc_figure_commands.md.")


if __name__ == "__main__":
    main()
