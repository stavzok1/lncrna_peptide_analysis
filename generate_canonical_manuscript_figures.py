"""
Build **main-text** manuscript figures only (outputs under ``figures/``).

Does **not** write supplementary panels under ``figures/supplementary/`` — use
``generate_supplementary_figures.py`` for those.

**Panels:**

- **Fig 1B** — sklearn 2D t-SNE, combined A|B (`fig1b_*_dims12_panels_AB.png`)
- **Fig 2** (tcga-matrix) + **Fig 2B** at ``figures/fig2b_stage_E_L_combined.png``
- **Fig 3A–3B–3C** (TCGA-matrix) + **Fig 3D** (Tr MPs)
- **Fig 4A**
- **Fig 5** merged **5A–5C** + **5D–5E** (proportional-whole coding cohort only)
- **Fig 6** TTN-AS1 split **instances** panels

NetMHC + Fig 6 are delegated to ``generate_netmhc_figure_bundle.py --canonical-main-text-only``.

Usage::

    python generate_canonical_manuscript_figures.py
    python generate_canonical_manuscript_figures.py --strict
    python generate_canonical_manuscript_figures.py --skip-netmhc --skip-iedb-pipeline
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import MANUSCRIPT_DIR, REPO_ROOT

MS = MANUSCRIPT_DIR


def run_ms(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(MS / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def run_root(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(REPO_ROOT / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="Exit non-zero on any failed subprocess.")
    ap.add_argument("--skip-netmhc", action="store_true", help="Skip generate_netmhc_figure_bundle.py entirely.")
    ap.add_argument(
        "--skip-iedb-pipeline",
        action="store_true",
        help="Forwarded to generate_netmhc_figure_bundle.py (skip merged *_with_iedb.tsv steps and Fig 6).",
    )
    args = ap.parse_args()

    failures: list[tuple[str, int]] = []

    def record(tag: str, script: str, code: int) -> None:
        if code != 0:
            failures.append((tag, script, code))

    code = run_ms("plot_figure1b_tsne_stage_lncrna.py", ["--embedding", "sklearn2_pca34"])
    record("fig1b", "plot_figure1b_tsne_stage_lncrna.py", code)

    for script, extra in (
        ("plot_tr_de_peptide_fractions_by_transition.py", ["--peptide-gene-set", "tcga_matrix"]),
        ("plot_aa_frequency_tcga_vs_proteome.py", ["--peptide-set", "tcga_matrix"]),
        ("plot_dipeptide_volcano_lnc_vs_proteome.py", ["--peptide-set", "tcga_matrix"]),
    ):
        c = run_ms(script, extra)
        record("catalog", script, c)

    c3 = run_ms("plot_figure3cd_dipeptide_log2fc_heatmaps.py", ["--only-tcga-matrix-3c"])
    record("fig3cd", "plot_figure3cd_dipeptide_log2fc_heatmaps.py", c3)

    c4 = run_ms("plot_figure4a_tis_vs_ribo_tr_mps.py", [])
    record("fig4a", "plot_figure4a_tis_vs_ribo_tr_mps.py", c4)

    if not args.skip_netmhc:
        bundle_extra = ["--canonical-main-text-only"]
        if args.skip_iedb_pipeline:
            bundle_extra.append("--skip-iedb-pipeline")
        if args.strict:
            bundle_extra.append("--strict")
        c_bundle = run_root("generate_netmhc_figure_bundle.py", bundle_extra)
        record("netmhc", "generate_netmhc_figure_bundle.py", c_bundle)

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with non-fatal failures:", failures)
    else:
        print("Canonical manuscript figures done. Outputs under figures/. See docs/figure_catalog.md.")


if __name__ == "__main__":
    main()
