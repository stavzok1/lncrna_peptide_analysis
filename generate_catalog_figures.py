"""
Build **Figure 1B** (t-SNE), then catalog **Figure 2** and **Figure 3** outputs for both peptide modes:

  - ``figures/tcga_matrix/`` — TCGA-matrix SmProt gene list / FASTA
  - ``figures/all_smprot_filtered/`` — full ``smprot_filtered.tsv`` genes / all-filtered FASTA

Runs **first**:

  0. ``manuscript/plot_figure1b_tsne_stage_lncrna.py`` — Fig. 1B (requires ``data/primary_exp_stage_lncRNA.csv``; **openTSNE**).

Then, **in order for each mode**:

  1. ``manuscript/plot_tr_de_peptide_fractions_by_transition.py --peptide-gene-set <mode>``
  2. ``manuscript/plot_aa_frequency_tcga_vs_proteome.py --peptide-set <mode>``
  3. ``manuscript/plot_dipeptide_volcano_lnc_vs_proteome.py --peptide-set <mode>``

Then **once** (shared under ``figures/``):

  4. ``manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py`` — Fig. 3C and Fig. 3D.

  5. ``manuscript/plot_figure4a_tis_vs_ribo_tr_mps.py`` — Fig. 4A.

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
    args = ap.parse_args()

    modes: list[str]
    if args.only == "both":
        modes = ["tcga_matrix", "all_smprot_filtered"]
    else:
        modes = [args.only]

    failures: list[tuple[str, str, int]] = []

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
