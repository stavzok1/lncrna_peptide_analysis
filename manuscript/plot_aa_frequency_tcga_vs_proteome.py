"""
Pooled 1-mer (amino-acid) frequency bars: lncRNA micropeptides (FASTA) vs human coding
reference proteome (``data/known_proteins.fasta``).

**Manuscript Figure 3A (default):** ``data/smprot_tcga_filtered_peptides.faa`` (TCGA-matrix
SmProt export).

**Alternate:** ``--peptide-set all_smprot_filtered`` uses ``data/smprot_all_filtered_peptides.faa``
(build with ``export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv
--out-aa data/smprot_all_filtered_peptides.faa``).

For each standard amino acid, tests whether its proportion differs between the two
pooled sequences using Fisher's exact test on a 2×2 table (residue vs not-residue ×
set). P-values are Benjamini-Hochberg adjusted across the 20 tests (default for stars);
use ``--significance raw`` for uncorrected Fisher p < alpha on each residue.

Writes figure + CSV + short report under ``figures/<peptide_set>/`` by default. For the
default TCGA-matrix run (no ``--peptide-fa`` / no custom ``--out-dir``), a copy is also
written to ``figures/fig3a.png``.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES

ROOT = REPO_ROOT


import argparse
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats

import figure_palettes as pal

DATA = ROOT / "data"
TCGA_MATRIX_PEPTIDE_FA = DATA / "smprot_tcga_filtered_peptides.faa"
ALL_SMPROT_FILTERED_FA = DATA / "smprot_all_filtered_peptides.faa"
DEFAULT_PROTEOME_FA = DATA / "known_proteins.fasta"
FIGURES = ROOT / "figures"

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA20)


def count_aa_in_sequence(seq: str, counts: np.ndarray) -> int:
    """Uppercase, U→C; count standard 20. Returns number of counted residues."""
    s = str(seq).upper().replace("U", "C")
    n = 0
    for ch in s:
        if ch in AA_SET:
            counts[AA20.index(ch)] += 1
            n += 1
    return n


def count_fasta(path: Path) -> tuple[np.ndarray, int]:
    """Return (length-20 counts, total aa counted)."""
    c = np.zeros(20, dtype=np.int64)
    total = 0
    for rec in SeqIO.parse(path, "fasta"):
        total += count_aa_in_sequence(str(rec.seq), c)
    return c, int(total)


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Return BH-adjusted p-values (same length as p)."""
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, dtype=np.float64)
    ranked[order] = np.arange(1, n + 1)
    q = p * n / ranked
    # cumulative min from largest rank downward
    adj = np.empty(n, dtype=np.float64)
    adj_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    return adj


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "AA frequency bar plot (Fig. 3A): lncRNA peptide FASTA vs known_proteins.fasta + FDR. "
            "Default is TCGA-matrix filtered FASTA."
        )
    )
    ap.add_argument(
        "--peptide-set",
        choices=("tcga_matrix", "all_smprot_filtered"),
        default="tcga_matrix",
        help=(
            "tcga_matrix: data/smprot_tcga_filtered_peptides.faa (default). "
            "all_smprot_filtered: data/smprot_all_filtered_peptides.faa (full SmProt-filtered export)."
        ),
    )
    ap.add_argument(
        "--peptide-fa",
        type=Path,
        default=None,
        help="Override FASTA path (ignores --peptide-set when set).",
    )
    ap.add_argument(
        "--proteome-fa",
        type=Path,
        default=DEFAULT_PROTEOME_FA,
        help="Reference proteome FASTA (default: data/known_proteins.fasta).",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: figures/<peptide_set>/ under repo root when --peptide-fa is unset).",
    )
    ap.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance threshold for stars on x labels (see --significance).",
    )
    ap.add_argument(
        "--significance",
        choices=("fdr", "raw"),
        default="fdr",
        help="fdr: star when BH-FDR < alpha (default). raw: star when Fisher p < alpha.",
    )
    args = ap.parse_args()

    if args.peptide_fa is not None:
        peptide_fa = args.peptide_fa
        out_stem = args.peptide_fa.stem
    elif args.peptide_set == "tcga_matrix":
        peptide_fa = TCGA_MATRIX_PEPTIDE_FA
        # Historical filenames used "tcga_filtered" for this FASTA.
        out_stem = "tcga_filtered"
    else:
        peptide_fa = ALL_SMPROT_FILTERED_FA
        out_stem = "all_smprot_filtered"

    if not peptide_fa.exists():
        hint = ""
        if peptide_fa == ALL_SMPROT_FILTERED_FA:
            hint = (
                "\nCreate it with:\n  python export_tcga_filtered_peptides_fasta.py "
                "--peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa"
            )
        raise SystemExit(f"Missing peptide FASTA: {peptide_fa}{hint}")
    if not args.proteome_fa.exists():
        raise SystemExit(f"Missing proteome FASTA: {args.proteome_fa}")

    if args.out_dir is not None:
        out_dir = args.out_dir
    elif args.peptide_fa is None:
        out_dir = FIGURES / args.peptide_set
    else:
        out_dir = FIGURES / args.peptide_fa.stem

    out_dir.mkdir(parents=True, exist_ok=True)

    c_pep, n_pep = count_fasta(peptide_fa)
    c_ref, n_ref = count_fasta(args.proteome_fa)

    if n_pep <= 0 or n_ref <= 0:
        raise SystemExit("No standard amino acids counted in one or both inputs.")

    f_pep = c_pep.astype(np.float64) / n_pep
    f_ref = c_ref.astype(np.float64) / n_ref

    p_raw = np.zeros(20, dtype=np.float64)
    odds = np.full(20, np.nan, dtype=np.float64)
    for i in range(20):
        a = int(c_pep[i])
        b = int(n_pep - c_pep[i])
        c_ = int(c_ref[i])
        d = int(n_ref - c_ref[i])
        # Rows: peptide / proteome; cols: this AA / other AA
        table = [[a, b], [c_, d]]
        oddsr, p_two = stats.fisher_exact(table, alternative="two-sided")
        p_raw[i] = p_two
        odds[i] = oddsr

    p_fdr = benjamini_hochberg(p_raw)
    sig_fdr = p_fdr < args.alpha
    sig_raw = p_raw < args.alpha
    sig = sig_fdr if args.significance == "fdr" else sig_raw

    df = pd.DataFrame(
        {
            "aa": AA20,
            "count_lnc_peptide_fasta": c_pep,
            "count_known_proteins": c_ref,
            "freq_lnc_peptide_fasta": f_pep,
            "freq_known_proteins": f_ref,
            "delta_lnc_minus_known": f_pep - f_ref,
            "fisher_p_two_sided": p_raw,
            "fisher_p_fdr_bh": p_fdr,
            "significant_raw_p_lt_alpha": sig_raw,
            "significant_fdr_bh_lt_alpha": sig_fdr,
            "odds_ratio_peptide_vs_proteome": odds,
            "peptide_fasta": str(peptide_fa),
            "peptide_set": args.peptide_set if args.peptide_fa is None else "custom_fasta",
        }
    )
    csv_path = out_dir / f"aa_frequency_{out_stem}_vs_known_proteins_stats.csv"
    df.to_csv(csv_path, index=False)

    if out_stem == "tcga_filtered":
        lnc_bar_label = "TCGA-matrix lncRNA peptides"
    elif out_stem == "all_smprot_filtered":
        lnc_bar_label = "All SmProt-filtered peptides"
    else:
        lnc_bar_label = f"lncRNA peptides ({out_stem})"

    x = np.arange(20)
    w = 0.36
    fig, ax = plt.subplots(figsize=(11, 4.8), dpi=150)
    ax.bar(
        x - w / 2,
        f_pep,
        width=w,
        label=lnc_bar_label,
        color=pal.AA_FREQ_TCGA,
        edgecolor="white",
        linewidth=0.35,
    )
    ax.bar(
        x + w / 2,
        f_ref,
        width=w,
        label="Known proteins",
        color=pal.AA_FREQ_PROTEOME,
        edgecolor="white",
        linewidth=0.35,
    )
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Amino acid")
    ax.set_xticks(x)
    labels = [f"{a}*" if sig[i] else a for i, a in enumerate(AA20)]
    ax.set_xticklabels(labels)
    ax.legend(frameon=False, loc="upper right")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, linestyle=":", alpha=0.6)
    ax.set_ylim(0, max(f_pep.max(), f_ref.max()) * 1.12)
    ax.set_title(
        "Amino acid frequency comparison\n"
        f"{lnc_bar_label}: {n_pep:,} AAs | Known proteins: {n_ref:,} AAs"
    )
    fig.tight_layout()
    png_path = out_dir / f"aa_frequency_{out_stem}_vs_known_proteins.png"
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)

    if args.peptide_fa is None and args.peptide_set == "tcga_matrix" and args.out_dir is None:
        fig3a_root = FIGURES / "fig3a.png"
        shutil.copy2(png_path, fig3a_root)
        print(f"Copied to {fig3a_root} (Fig 3A)")

    report = out_dir / f"aa_frequency_{out_stem}_vs_known_proteins_report.txt"
    lines = [
        f"peptide_fasta: {peptide_fa}",
        f"proteome_fasta: {args.proteome_fa}",
        f"n_aa_lnc_peptide_fasta: {n_pep:,}",
        f"n_aa_known_proteins: {n_ref:,}",
        f"per-residue test: Fisher exact (2×2), two-sided",
        "multiple testing: Benjamini-Hochberg FDR over 20 amino acids (reported in CSV)",
        f"significance for x-axis '*': {args.significance} p < {args.alpha}",
        f"n_significant_for_plot ({args.significance}): {int(sig.sum())}",
        f"csv: {csv_path}",
        f"figure: {png_path}",
    ]
    if args.peptide_fa is None and args.peptide_set == "tcga_matrix" and args.out_dir is None:
        lines.append(f"fig3a_root_copy: {FIGURES / 'fig3a.png'}")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
