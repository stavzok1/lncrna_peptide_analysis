"""
Build catalog Figure 2 and Figure 3 outputs for both peptide modes:

  - ``figures/tcga_matrix/`` — TCGA-matrix SmProt gene list / FASTA
  - ``figures/all_smprot_filtered/`` — full ``smprot_filtered.tsv`` genes / all-filtered FASTA

Runs, in order for each mode:

  1. ``plot_tr_de_peptide_fractions_by_transition.py --peptide-gene-set <mode>``
  2. ``plot_aa_frequency_tcga_vs_proteome.py --peptide-set <mode>``
  3. ``plot_dipeptide_volcano_lnc_vs_proteome.py --peptide-set <mode>``

Then **once** (shared under ``figures/``):

  4. ``plot_figure3cd_dipeptide_log2fc_heatmaps.py`` — Fig. 3C (separate TCGA-matrix and
     all-filtered log2FC heatmaps vs proteome) and Fig. 3D (Tr-lncRNA TCGA MPs vs proteome).

  5. ``plot_figure4a_tis_vs_ribo_tr_mps.py`` — Fig. 4A (TIS vs Ribo-seq p-values for Tr MPs).

``all_smprot_filtered`` requires ``data/smprot_all_filtered_peptides.faa`` (see
``export_tcga_filtered_peptides_fasta.py``). Failures for that mode are reported but
do not stop the TCGA-matrix mode unless ``--strict`` is passed.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ALL_FILTERED_FAA = ROOT / "data" / "smprot_all_filtered_peptides.faa"


def run_script(script: str, extra: list[str]) -> int:
    cmd = [sys.executable, str(ROOT / script), *extra]
    print("+", " ".join(cmd))
    return subprocess.call(cmd, cwd=str(ROOT))


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate figures/tcga_matrix and figures/all_smprot_filtered.")
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
    for mode in modes:
        if mode == "all_smprot_filtered" and not ALL_FILTERED_FAA.exists():
            print(
                f"Skip {mode}: missing {ALL_FILTERED_FAA}\n"
                "  Build with: python export_tcga_filtered_peptides_fasta.py "
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
        print("Done. Outputs under", ROOT / "figures")


if __name__ == "__main__":
    main()
