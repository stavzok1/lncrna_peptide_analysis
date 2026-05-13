#!/usr/bin/env python3
"""
Figure 5 (merged NetMHC + IEDB-style gates): **significant lncRNA** vs **proportional-whole coding**.

1. Builds ``data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv`` from the wide
   ``netmhcpan_coding_proportional_whole.xls`` using ``merge_netmhcpan_xls_with_iedb.py
   --synthetic-iedb-pass`` (no real IEDB API run for ``COD|coding_pwhole|…`` headers;
   synthetic imm/proc + BA→IC50 so ``plot_fig5abc_netmhc_sb_triple.py`` full SB matches
   manuscript defaults).

2. Runs ``plot_fig5abc_netmhc_sb_triple.py`` and ``plot_fig5de_merged_iedb_sb_per_allele.py``
   writing under ``figures/manuscript_netmhc/fig5_sig_vs_proportional_whole/`` by default.

Prerequisites
-------------
- ``data/netmhc/netmhcpan_coding_proportional_whole.xls`` (local NetMHCpan on proportional FASTA)
- ``data/netmhc/ninemers_coding_proportional_whole.fasta`` (same run / prep as the XLS)
- ``data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv`` (significant cohort, real IEDB merge)

Example::

  python supplement/run_fig5_sig_vs_proportional_coding.py
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import FIGURES, MANUSCRIPT_DIR, NETMHC_DATA, REPO_ROOT, SCRIPTS_DIR

ROOT = REPO_ROOT


import argparse
import subprocess


def run(cmd: list[str], *, cwd: Path) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=str(cwd))


def main() -> None:
    root = REPO_ROOT
    net = NETMHC_DATA
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--coding-xls",
        type=Path,
        default=net / "netmhcpan_coding_proportional_whole.xls",
        help="Wide NetMHCpan XLS for proportional-whole coding 9-mers.",
    )
    ap.add_argument(
        "--coding-ninemers-fasta",
        type=Path,
        default=net / "ninemers_coding_proportional_whole.fasta",
        help="FASTA that matches the XLS row order (COD|coding_pwhole|… headers).",
    )
    ap.add_argument(
        "--coding-out-tsv",
        type=Path,
        default=net / "netmhcpan_coding_proportional_whole_with_iedb.tsv",
        help="Output merged long TSV for the coding cohort.",
    )
    ap.add_argument(
        "--sig-tsv",
        type=Path,
        default=net / "netmhcpan_sig_lnc_with_iedb.tsv",
        help="Merged significant lncRNA TSV (real IEDB).",
    )
    ap.add_argument(
        "--allele-freq-csv",
        type=Path,
        default=net / "figures" / "fig5a_epitopes_vs_allele_frequency_ic50_sb.csv",
        help="Allele frequency table for Fig 5A (same as main manuscript pipeline).",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=FIGURES / "manuscript_netmhc" / "fig5_sig_vs_proportional_whole",
        help="Folder for Fig 5 PNG/CSV outputs.",
    )
    ap.add_argument(
        "--output-stem",
        type=str,
        default="fig5abc_sig_vs_proportional_whole",
        help="Filename prefix for 5A–5C (and base for 5D–5E stems).",
    )
    ap.add_argument(
        "--skip-merge",
        action="store_true",
        help="Do not rebuild the coding merged TSV (use existing --coding-out-tsv).",
    )
    args = ap.parse_args()

    if not args.coding_xls.is_file():
        raise SystemExit(f"Missing {args.coding_xls}")
    if not args.coding_ninemers_fasta.is_file():
        raise SystemExit(f"Missing {args.coding_ninemers_fasta}")
    if not args.sig_tsv.is_file():
        raise SystemExit(f"Missing {args.sig_tsv}")
    if not args.allele_freq_csv.is_file():
        raise SystemExit(f"Missing {args.allele_freq_csv} (generate via wide-XLS Fig 5A pipeline if needed)")

    py = sys.executable
    merge_py = SCRIPTS_DIR / "merge_netmhcpan_xls_with_iedb.py"
    figabc = MANUSCRIPT_DIR / "plot_fig5abc_netmhc_sb_triple.py"
    figde = MANUSCRIPT_DIR / "plot_fig5de_merged_iedb_sb_per_allele.py"

    if not args.skip_merge:
        run(
            [
                py,
                str(merge_py),
                "--netmhc-xls",
                str(args.coding_xls),
                "--synthetic-iedb-pass",
                "--join",
                "stable",
                "--ninemers-fasta",
                str(args.coding_ninemers_fasta),
                "--out-tsv",
                str(args.coding_out_tsv),
            ],
            cwd=root,
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    readme = args.out_dir / "README.txt"
    readme.write_text(
        "Figure 5 outputs: significant lncRNA (netmhcpan_sig_lnc_with_iedb.tsv) vs "
        "proportional-whole coding (netmhcpan_coding_proportional_whole_with_iedb.tsv).\n"
        "Coding cohort uses synthetic IEDB pass columns from merge_netmhcpan_xls_with_iedb.py "
        "--synthetic-iedb-pass (not real IEDB API scores).\n"
        f"Generated by: python supplement/run_fig5_sig_vs_proportional_coding.py\n",
        encoding="utf-8",
    )

    run(
        [
            py,
            str(figabc),
            "--sig-tsv",
            str(args.sig_tsv),
            "--coding-tsv",
            str(args.coding_out_tsv),
            "--allele-freq-csv",
            str(args.allele_freq_csv),
            "--out-dir",
            str(args.out_dir),
            "--output-stem",
            args.output_stem,
        ],
        cwd=root,
    )
    run(
        [
            py,
            str(figde),
            "--sig-tsv",
            str(args.sig_tsv),
            "--coding-tsv",
            str(args.coding_out_tsv),
            "--out-dir",
            str(args.out_dir),
            "--output-stem",
            "fig5de_sig_vs_proportional_whole",
            "--no-repo-mirror",
        ],
        cwd=root,
    )
    print(f"Done. Outputs under {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
