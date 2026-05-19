"""
Build **main-text** manuscript figures in one pass (single driver).

**Canonical** (repo ``figures/``):

- **Fig 1B** — default ``sklearn2_pca34`` (**2D sklearn t-SNE only**, two PNGs under ``figures/``); use
  ``--fig1b-embedding opentsne4`` for four t-SNE panels (dims 1–2 and 3–4).
- **Fig 2** (tcga-matrix mode only) + **Fig 2B** copy at ``figures/fig2b_stage_E_L_combined.png``
- **Fig 3A–3B–3C** (TCGA-matrix lncRNA MPs only) + **Fig 3D** (Tr MPs; TCGA-matrix peptide set)
- **Fig 4A** — ``significant_lnc_peptides.tsv`` (~501 MPs; default in ``plot_figure4a_tis_vs_ribo_tr_mps.py``)
- **Fig 5** merged — **5A–5C** + **5D–5E** for **proportional-whole** coding cohort only (no random-fragment mirrors)
- **Fig 6** TTN-AS1 — ``--split-panels`` **instances** only (under ``figures/fig6_ttn_as1_split_instances_*.png``)

**Supplementary** (``figures/supplementary/``):

- **Fig 1B** — **OpenTSNE** four panels (dims 1–2 and 3–4) under ``supplementary/embedding/``, unless
  OpenTSNE is already the canonical embedding or you pass ``--skip-opentsne-supplement``.
- **Fig 6** — split panels for **unique** coverage under ``supplementary/figure6_ttn_as1/`` only when you pass
  ``--write-fig6-unique-supplement`` (default: **do not** write Fig 6 unique anywhere).

NetMHC + Fig 6 core are delegated to ``generate_netmhc_figure_bundle.py --canonical-main-text-only``.

Usage::

    python generate_canonical_manuscript_figures.py
    python generate_canonical_manuscript_figures.py --strict
    python generate_canonical_manuscript_figures.py --fig1b-embedding opentsne4
    python generate_canonical_manuscript_figures.py --skip-netmhc --skip-iedb-pipeline
    python generate_canonical_manuscript_figures.py --write-fig6-unique-supplement
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import (
    FIGURES_SUPPLEMENTARY_EMBEDDING,
    FIGURES_SUPPLEMENTARY_FIG6_TTN,
    MANUSCRIPT_DIR,
    REPO_ROOT,
)

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
    ap.add_argument(
        "--fig1b-embedding",
        choices=("sklearn2_pca34", "opentsne4"),
        default="sklearn2_pca34",
        help="Embedding for Fig 1B under figures/: sklearn = 2D t-SNE only (two PNGs); opentsne4 = four PNGs.",
    )
    ap.add_argument(
        "--skip-opentsne-supplement",
        action="store_true",
        help="Do not write OpenTSNE 1–2 / 3–4 panels under figures/supplementary/embedding/.",
    )
    ap.add_argument(
        "--write-fig6-unique-supplement",
        action="store_true",
        help="Write Fig 6 split unique panels under figures/supplementary/figure6_ttn_as1/ (default: off).",
    )
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

    code = run_ms(
        "plot_figure1b_tsne_stage_lncrna.py",
        ["--embedding", args.fig1b_embedding],
    )
    record("fig1b", "plot_figure1b_tsne_stage_lncrna.py", code)

    want_opentsne_supp = (not args.skip_opentsne_supplement) and args.fig1b_embedding != "opentsne4"
    if want_opentsne_supp:
        FIGURES_SUPPLEMENTARY_EMBEDDING.mkdir(parents=True, exist_ok=True)
        code_ot = run_ms(
            "plot_figure1b_tsne_stage_lncrna.py",
            [
                "--embedding",
                "opentsne4",
                "--out-dir",
                str(FIGURES_SUPPLEMENTARY_EMBEDDING),
                "--filename-prefix",
                "figS1b_opentsne4_tsne_stage_lncrna_samples",
            ],
        )
        record("fig1b_supp", "plot_figure1b_tsne_stage_lncrna.py (opentsne4)", code_ot)

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

    if args.write_fig6_unique_supplement:
        FIGURES_SUPPLEMENTARY_FIG6_TTN.mkdir(parents=True, exist_ok=True)
        ttn_sup = FIGURES_SUPPLEMENTARY_FIG6_TTN / "fig6_ttn_as1_split.png"
        c6u = run_ms(
            "plot_figure6_ttn_as1_allele_coverage.py",
            [
                "--split-panels",
                "--coverage-output",
                "unique",
                "-o",
                str(ttn_sup),
            ],
        )
        record("fig6_supp", "plot_figure6_ttn_as1_allele_coverage.py (unique)", c6u)

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with non-fatal failures:", failures)
    else:
        print(
            "Canonical manuscript figures done. Main outputs under figures/; alternates under "
            "figures/supplementary/. See figures/supplementary/README.md."
        )


if __name__ == "__main__":
    main()
