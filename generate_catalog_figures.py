"""
Build catalog **Figure 2** and **Figure 3** outputs (and optionally main-text panels).

**Default (full catalog):** Fig 1B + Fig 2–3 for both peptide modes under
``figures/supplementary/<mode>/``, plus Fig 3C–3D and 4A at ``figures/`` root.

**``--supplement-only``:** only supplementary trees (both modes) + all-filtered Fig 3C when the
FASTA exists — no Fig 1B, 3D, or 4A at repo root. Used by ``generate_supplementary_figures.py``.

For a clean split, prefer:

- ``generate_canonical_manuscript_figures.py`` — main text at ``figures/``
- ``generate_supplementary_figures.py`` — everything under ``figures/supplementary/``

``all_smprot_filtered`` requires ``data/smprot_all_filtered_peptides.faa`` (see
``pipeline/export_tcga_filtered_peptides_fasta.py``).
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import DATA, MANUSCRIPT_DIR, REPO_ROOT

MS = MANUSCRIPT_DIR
ALL_FILTERED_FAA = DATA / "smprot_all_filtered_peptides.faa"


def run_script(script: str, extra: list[str]) -> int:
    cmd = [sys.executable, str(MS / script), *extra]
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit with non-zero status if any subprocess fails (including all_smprot_filtered).",
    )
    ap.add_argument(
        "--only",
        choices=("tcga_matrix", "all_smprot_filtered", "both"),
        default="both",
        help="Which peptide mode(s) to run.",
    )
    ap.add_argument(
        "--supplement-only",
        action="store_true",
        help="Only Fig 2–3 supplement trees (+ all-filtered 3C); skip Fig 1B, 3D, and 4A at figures/ root.",
    )
    args = ap.parse_args()

    modes: list[str]
    if args.only == "both":
        modes = ["tcga_matrix", "all_smprot_filtered"]
    else:
        modes = [args.only]

    failures: list[tuple[str, str, int]] = []

    if not args.supplement_only:
        code_1b = run_script("plot_figure1b_tsne_stage_lncrna.py", [])
        if code_1b != 0:
            failures.append(("fig1b", "plot_figure1b_tsne_stage_lncrna.py", code_1b))

    for mode in modes:
        if mode == "all_smprot_filtered" and not ALL_FILTERED_FAA.exists():
            print(
                f"Skip {mode}: missing {ALL_FILTERED_FAA}\n"
                "  Build with: python pipeline/export_tcga_filtered_peptides_fasta.py "
                "--peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa"
            )
            if args.strict:
                print(f"Missing required FASTA: {ALL_FILTERED_FAA}", file=sys.stderr)
                sys.exit(1)
            continue

        for script, extra in (
            ("plot_tr_de_peptide_fractions_by_transition.py", ["--peptide-gene-set", mode]),
            ("plot_aa_frequency_tcga_vs_proteome.py", ["--peptide-set", mode]),
            ("plot_dipeptide_volcano_lnc_vs_proteome.py", ["--peptide-set", mode]),
        ):
            code = run_script(script, extra)
            if code != 0:
                failures.append((mode, script, code))

    if args.supplement_only:
        if ALL_FILTERED_FAA.exists():
            code_3c_af = run_script(
                "plot_figure3cd_dipeptide_log2fc_heatmaps.py",
                ["--only-all-smprot-filtered-3c"],
            )
            if code_3c_af != 0:
                failures.append(("shared", "plot_figure3cd_dipeptide_log2fc_heatmaps.py (all-filtered 3C)", code_3c_af))
        elif args.strict:
            print(f"Missing required FASTA for all-filtered 3C: {ALL_FILTERED_FAA}", file=sys.stderr)
            sys.exit(1)
        else:
            print(f"Skip all-filtered Fig 3C: missing {ALL_FILTERED_FAA}")
    else:
        code_3cd = run_script("plot_figure3cd_dipeptide_log2fc_heatmaps.py", [])
        if code_3cd != 0:
            failures.append(("shared", "plot_figure3cd_dipeptide_log2fc_heatmaps.py", code_3cd))

        code_4a = run_script("plot_figure4a_tis_vs_ribo_tr_mps.py", [])
        if code_4a != 0:
            failures.append(("shared", "plot_figure4a_tis_vs_ribo_tr_mps.py", code_4a))

    if failures and args.strict:
        print("Failures:", failures)
        sys.exit(1)
    if failures:
        print("Completed with non-fatal failures (re-run with --strict to exit non-zero):", failures)
    else:
        print("Done. Outputs under", REPO_ROOT / "figures")


if __name__ == "__main__":
    main()
