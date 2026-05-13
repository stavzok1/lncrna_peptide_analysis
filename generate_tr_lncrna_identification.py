"""
Tr-lncRNA identification from TCGA lncRNA expression matrices.

Runs, in order:

1. **Z-score screen** — ``pipeline/tr_lncrna_de_analysis.py``  
   Per cancer × (AJCC stage or M_stage) transition: genes with expression filter, log2FC
   from group means, z across genes, **|z| ≥ 3**. Writes ``tr_lncrna_output/tr_lncrnas_*_detail.csv``,
   ``tr_genes_union.txt`` (z union), summary JSON, diagnostic PNGs.

2. **Limma DE** — ``tr_limma_de.R`` (R: limma + eBayes, BH-FDR < 0.05).  
   Always writes ``tr_lncrna_output/limma/`` tables. **Only if** ``tr_lncrna_output/tr_genes_union.txt``
   exists (from step 1) does R also write ``limma_z_intersection_genes.txt`` and, when non-empty,
   ``data/canonical_significant_lncRNAs.txt``. If that file is missing, run step 1 first, then re-run R.

**Prerequisites:** ``data/primary_exp_stage_lncRNA.csv``, ``data/primary_exp_metastasis_lncRNA.csv``.  
**R:** install once with ``Rscript install_r_dependencies.R``; ``Rscript`` on ``PATH`` (or ``--rscript``).

**Downstream (not run here):** SmProt peptide table — ``python pipeline/build_significant_lncs_smprot.py``
(needs ``data/SmProt2.txt``, ``data/lncrna_genes_small.csv``, etc.).

Usage::

    python generate_tr_lncrna_identification.py
    python generate_tr_lncrna_identification.py --strict
    python generate_tr_lncrna_identification.py --skip-r   # z-scores only
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from orchestrate_subprocess import call_echo
from repo_paths import REPO_ROOT

ROOT = REPO_ROOT
PY = sys.executable
Z_SCRIPT = ROOT / "pipeline" / "tr_lncrna_de_analysis.py"
R_SCRIPT = ROOT / "tr_limma_de.R"


def run(cmd: list[str], *, cwd: Path) -> int:
    return call_echo(cmd, cwd=cwd)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="Exit non-zero if any step fails.")
    ap.add_argument("--skip-r", action="store_true", help="Skip tr_limma_de.R (z-score outputs only).")
    ap.add_argument(
        "--rscript",
        type=str,
        default="Rscript",
        help="Rscript executable name or path (default: Rscript).",
    )
    args = ap.parse_args()
    failures: list[tuple[str, int]] = []

    stage_csv = ROOT / "data" / "primary_exp_stage_lncRNA.csv"
    meta_csv = ROOT / "data" / "primary_exp_metastasis_lncRNA.csv"
    for p in (stage_csv, meta_csv):
        if not p.is_file():
            msg = f"Missing required matrix: {p}"
            print(msg, file=sys.stderr)
            if args.strict:
                sys.exit(1)
            print("Abort (non-strict): cannot run z-score step.", file=sys.stderr)
            return

    code_z = run([PY, str(Z_SCRIPT)], cwd=ROOT)
    if code_z != 0:
        failures.append(("tr_lncrna_de_analysis.py", code_z))

    if not args.skip_r:
        rscript = shutil.which(args.rscript) or args.rscript
        code_r = run([rscript, str(R_SCRIPT), str(ROOT)], cwd=ROOT)
        if code_r != 0:
            failures.append(("tr_limma_de.R", code_r))
    elif args.strict:
        print("[warn] --skip-r: limma outputs and canonical_significant_lncRNAs.txt not updated.", flush=True)

    if failures and args.strict:
        print("Failures:", failures, file=sys.stderr)
        sys.exit(1)
    if failures:
        print("Completed with failures:", failures)
    else:
        print("Tr-lncRNA identification done.")
        print("  Z detail:", ROOT / "tr_lncrna_output" / "tr_lncrnas_stage_detail.csv")
        print("  Limma:  ", ROOT / "tr_lncrna_output" / "limma")
        print("  Canonical genes:", ROOT / "data" / "canonical_significant_lncRNAs.txt", "(if limma∩z non-empty)")


if __name__ == "__main__":
    main()
