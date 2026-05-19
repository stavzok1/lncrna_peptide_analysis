"""
Build **all supplementary** figure outputs under ``figures/supplementary/``.

Does **not** write main-text panels at ``figures/`` root — use
``generate_canonical_manuscript_figures.py`` for those.

**Steps (in order):**

1. OpenTSNE Figure 1 supplement (``embedding/``)
2. Sample PCA (``pca/``)
3. Fig 2–3 for both peptide modes (``tcga_matrix/``, ``all_smprot_filtered/``) + all-filtered Fig 3C
4. NetMHC random-fragment cohort mirrors (``netmhc/coding_fragments_random_sample/``)
5. NetMHC Fig 5–6 sensitivity tree (``netmhc_fig5_fig6_supplement/``)
6. Optional Fig 6 **unique** split panels (``figure6_ttn_as1/``)

**Prerequisites:** merged ``*_with_iedb.tsv`` under ``data/netmhc/`` for NetMHC steps.

Usage::

    python generate_supplementary_figures.py
    python generate_supplementary_figures.py --strict
    python generate_supplementary_figures.py --include-fig6-unique
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import (
    DATA,
    FIGURES_SUPPLEMENTARY_EMBEDDING,
    FIGURES_SUPPLEMENTARY_FIG6_TTN,
    FIGURES_SUPPLEMENTARY_PCA,
    MANUSCRIPT_DIR,
    REPO_ROOT,
)

MS = MANUSCRIPT_DIR
SUP = REPO_ROOT / "supplement"
ALL_FILTERED_FAA = DATA / "smprot_all_filtered_peptides.faa"
PEPTIDE_MODES = ("tcga_matrix", "all_smprot_filtered")


def run_ms(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(MS / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_sup(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(SUP / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_root(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(REPO_ROOT / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_fig2_3_catalog(*, strict: bool) -> list[tuple[str, int]]:
    """Fig 2–3 supplement trees for tcga_matrix and all_smprot_filtered (+ all-filtered 3C)."""
    failures: list[tuple[str, int]] = []

    for mode in PEPTIDE_MODES:
        if mode == "all_smprot_filtered" and not ALL_FILTERED_FAA.exists():
            print(
                f"Skip {mode}: missing {ALL_FILTERED_FAA}\n"
                "  Build with: python pipeline/export_tcga_filtered_peptides_fasta.py "
                "--peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa"
            )
            if strict:
                print(f"Missing required FASTA: {ALL_FILTERED_FAA}", file=sys.stderr)
                sys.exit(1)
            continue

        for script, extra in (
            ("plot_tr_de_peptide_fractions_by_transition.py", ["--peptide-gene-set", mode]),
            ("plot_aa_frequency_tcga_vs_proteome.py", ["--peptide-set", mode]),
            ("plot_dipeptide_volcano_lnc_vs_proteome.py", ["--peptide-set", mode]),
        ):
            code = run_ms(script, extra)
            if code != 0:
                failures.append((f"{mode}:{script}", code))

    if ALL_FILTERED_FAA.exists():
        code_3c = run_ms(
            "plot_figure3cd_dipeptide_log2fc_heatmaps.py",
            ["--only-all-smprot-filtered-3c"],
        )
        if code_3c != 0:
            failures.append(("all_filtered_3c", code_3c))
    elif strict:
        print(f"Missing required FASTA for all-filtered Fig 3C: {ALL_FILTERED_FAA}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"Skip all-filtered Fig 3C: missing {ALL_FILTERED_FAA}")

    return failures


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="Exit non-zero on any failed subprocess.")
    ap.add_argument(
        "--skip-netmhc-fig5-fig6-supplement",
        action="store_true",
        help="Skip generate_netmhc_fig5_fig6_supplement.py (sensitivity / Cartesian grids).",
    )
    ap.add_argument(
        "--skip-netmhc-random-fragments",
        action="store_true",
        help="Skip random-fragment NetMHC mirrors under netmhc/coding_fragments_random_sample/.",
    )
    ap.add_argument(
        "--include-fig6-unique",
        action="store_true",
        help="Write Fig 6 unique split panels under figures/supplementary/figure6_ttn_as1/.",
    )
    args = ap.parse_args()

    strict = ["--strict"] if args.strict else []
    failures: list[tuple[str, int]] = []

    def record(tag: str, code: int) -> None:
        if code != 0:
            failures.append((tag, code))

    FIGURES_SUPPLEMENTARY_EMBEDDING.mkdir(parents=True, exist_ok=True)
    record(
        "opentsne",
        run_ms(
            "plot_figure1_tsne_stage_lncrna.py",
            [
                "--embedding",
                "opentsne4",
                "--out-dir",
                str(FIGURES_SUPPLEMENTARY_EMBEDDING),
                "--filename-prefix",
                "figS1_opentsne4_tsne_stage_lncrna_samples",
            ],
        ),
    )

    record(
        "pca",
        run_sup(
            "plot_supplement_pca_stage_samples.py",
            ["--out-dir", str(FIGURES_SUPPLEMENTARY_PCA)],
        ),
    )

    failures.extend(run_fig2_3_catalog(strict=args.strict))

    if not args.skip_netmhc_random_fragments:
        record(
            "netmhc_fragments",
            run_root("generate_netmhc_figure_bundle.py", [*strict, "--supplement-mirrors-only"]),
        )

    if not args.skip_netmhc_fig5_fig6_supplement:
        record(
            "netmhc_fig5_fig6_supp",
            run_root("generate_netmhc_fig5_fig6_supplement.py", strict),
        )

    if args.include_fig6_unique:
        FIGURES_SUPPLEMENTARY_FIG6_TTN.mkdir(parents=True, exist_ok=True)
        ttn_sup = FIGURES_SUPPLEMENTARY_FIG6_TTN / "fig6_ttn_as1_split.png"
        record(
            "fig6_unique",
            run_ms(
                "plot_figure6_ttn_as1_allele_coverage.py",
                [
                    "--split-panels",
                    "--coverage-output",
                    "unique",
                    "-o",
                    str(ttn_sup),
                ],
            ),
        )

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        return 1
    if failures:
        print("Completed with non-fatal failures:", failures)
    else:
        print(
            "Supplementary figures done. Outputs under figures/supplementary/. "
            "See figures/supplementary/README.md."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
