"""
Build **all supplementary** figure outputs under ``figures/supplementary/``.

Does **not** write main-text panels at ``figures/`` root — use
``generate_canonical_manuscript_figures.py`` for those.

**Steps (in order):**

1. OpenTSNE Fig 1B (``embedding/``)
2. Sample PCA (``pca/``)
3. Fig 2–3 catalog for both peptide modes (``tcga_matrix/``, ``all_smprot_filtered/``) + all-filtered Fig 3C when FASTA exists
4. NetMHC random-fragment cohort mirrors (``netmhc/coding_fragments_random_sample/``)
5. NetMHC Fig 5–6 sensitivity tree (``netmhc_fig5_fig6_supplement/``)
6. Optional Fig 6 **unique** split panels (``figure6_ttn_as1/``)

**Prerequisites:** same inputs as main-text NetMHC (merged ``*_with_iedb.tsv`` under ``data/netmhc/``).
Run ``generate_canonical_manuscript_figures.py`` first if main-text Fig 5–6 are not yet built.

Usage::

    python generate_supplementary_figures.py
    python generate_supplementary_figures.py --strict
    python generate_supplementary_figures.py --include-fig6-unique
    python generate_supplementary_figures.py --skip-netmhc-fig5-fig6-supplement
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import (
    FIGURES_SUPPLEMENTARY_EMBEDDING,
    FIGURES_SUPPLEMENTARY_FIG6_TTN,
    FIGURES_SUPPLEMENTARY_PCA,
    MANUSCRIPT_DIR,
    REPO_ROOT,
)

MS = MANUSCRIPT_DIR
SUP = REPO_ROOT / "supplement"


def run_ms(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(MS / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_sup(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(SUP / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_root(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(REPO_ROOT / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


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
            "plot_figure1b_tsne_stage_lncrna.py",
            [
                "--embedding",
                "opentsne4",
                "--out-dir",
                str(FIGURES_SUPPLEMENTARY_EMBEDDING),
                "--filename-prefix",
                "figS1b_opentsne4_tsne_stage_lncrna_samples",
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

    record("catalog_supp", run_root("generate_catalog_figures.py", [*strict, "--supplement-only"]))

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
