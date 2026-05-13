"""
First-amino-acid (N-terminal) distribution for TCGA-filtered lncRNA micropeptides.

Sets from ``data/smprot_tcga_filtered_peptides.faa``:
  - **Filtered (all):** every record in the FASTA (TCGA-panel filtered MPs).
  - **Significant (canonical):** records whose ``GeneSymbol`` is in
    ``data/canonical_significant_lncRNAs.txt`` (limma ∩ z gene list).

When ``data/significant_lnc_peptides.faa`` exists, adds **SmProt significant exportable**
(first AA across that FASTA; same peptides as ``significant_lnc_peptides.tsv`` after sync).

Writes a bar chart + CSV + short text report under ``tr_lncrna_output/figures/``.
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



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats

import figure_palettes as pal

DATA = ROOT / "data"
MP_FAA = DATA / "smprot_tcga_filtered_peptides.faa"
SIG_FAA = DATA / "significant_lnc_peptides.faa"
CANONICAL = DATA / "canonical_significant_lncRNAs.txt"
LIMMA_Z = ROOT / "tr_lncrna_output" / "limma" / "limma_z_intersection_genes.txt"
OUT_DIR = ROOT / "tr_lncrna_output" / "figures"
OUT_FIG = OUT_DIR / "starting_aa_filtered_vs_significant.png"
OUT_CSV = OUT_DIR / "starting_aa_distribution_counts.csv"
OUT_TXT = OUT_DIR / "starting_aa_distribution_report.txt"

AA20 = list("ACDEFGHIKLMNPQRSTVWY")


def load_tr_genes() -> set[str]:
    path = CANONICAL if CANONICAL.exists() else LIMMA_Z
    if not path.exists():
        raise FileNotFoundError(f"Need {CANONICAL} or {LIMMA_Z}")
    return {ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()}


def first_aa_counts(faa_path: Path, gene_filter: set[str] | None) -> tuple[np.ndarray, int, int]:
    """Return (counts length 20 for AA20 order, n_used, n_nonstandard_or_empty)."""
    c = np.zeros(20, dtype=np.int64)
    n_used = 0
    n_skip = 0
    for rec in SeqIO.parse(faa_path, "fasta"):
        parts = rec.description.split("|")
        gene = parts[1].strip() if len(parts) > 1 else ""
        if gene_filter is not None and gene not in gene_filter:
            continue
        seq = str(rec.seq).upper().strip()
        if not seq:
            n_skip += 1
            continue
        a0 = seq[0]
        if a0 in AA20:
            c[AA20.index(a0)] += 1
            n_used += 1
        else:
            n_skip += 1
    return c, int(n_used), int(n_skip)


def first_aa_sig_export_fasta(faa_path: Path) -> tuple[np.ndarray, int, int]:
    """N-terminal AA counts for every record in significant_lnc_peptides.faa."""
    c = np.zeros(20, dtype=np.int64)
    n_used = 0
    n_skip = 0
    for rec in SeqIO.parse(faa_path, "fasta"):
        seq = str(rec.seq).upper().strip()
        if not seq:
            n_skip += 1
            continue
        a0 = seq[0]
        if a0 in AA20:
            c[AA20.index(a0)] += 1
            n_used += 1
        else:
            n_skip += 1
    return c, int(n_used), int(n_skip)


def main() -> None:
    if not MP_FAA.exists():
        raise FileNotFoundError(f"Missing {MP_FAA}")

    tr = load_tr_genes()
    c_all, n_all, skip_all = first_aa_counts(MP_FAA, None)
    c_sig, n_sig, skip_sig = first_aa_counts(MP_FAA, tr)

    p_all = c_all / max(n_all, 1)
    p_sig = c_sig / max(n_sig, 1)

    c_smprot = n_smprot = skip_smprot = None
    p_smprot = None
    if SIG_FAA.exists():
        c_smprot, n_smprot, skip_smprot = first_aa_sig_export_fasta(SIG_FAA)
        p_smprot = c_smprot / max(n_smprot, 1)

    # Homogeneity across sets present
    cols = [c_all.astype(float), c_sig.astype(float)]
    if c_smprot is not None:
        cols.append(c_smprot.astype(float))
    table = np.column_stack(cols)
    chi2, p_homog, dof, expected = stats.chi2_contingency(table)
    # vs uniform (each set separately)
    chi2_uni_all, p_uni_all = stats.chisquare(c_all)
    chi2_uni_sig, p_uni_sig = stats.chisquare(c_sig)
    chi2_uni_smp, p_uni_smp = (None, None)
    if c_smprot is not None:
        chi2_uni_smp, p_uni_smp = stats.chisquare(c_smprot)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame({"aa": AA20, "n_filtered_all": c_all, "pct_filtered_all": 100.0 * p_all})
    df["n_significant_canonical"] = c_sig
    df["pct_significant_canonical"] = 100.0 * p_sig
    df["delta_pct_sig_minus_all"] = 100.0 * (p_sig - p_all)
    if c_smprot is not None and p_smprot is not None:
        df["n_significant_smprot_export"] = c_smprot
        df["pct_significant_smprot_export"] = 100.0 * p_smprot
        df["delta_pct_smprot_minus_all"] = 100.0 * (p_smprot - p_all)
    df.to_csv(OUT_CSV, index=False)

    x = np.arange(20)
    fig, ax = plt.subplots(figsize=(11, 4.8))
    if c_smprot is not None and p_smprot is not None:
        w = 0.24
        ax.bar(x - w, 100 * p_all, width=w, label=f"Filtered (all FASTA, n={n_all})", color=pal.OI_SKY_BLUE, edgecolor=pal.OI_BLUE)
        ax.bar(x, 100 * p_sig, width=w, label=f"Canonical genes (n={n_sig})", color=pal.OI_ORANGE, edgecolor=pal.OI_VERMILLION)
        ax.bar(
            x + w,
            100 * p_smprot,
            width=w,
            label=f"SmProt significant exportable (n={n_smprot})",
            color=pal.OI_BLUISH_GREEN,
            edgecolor="#006666",
        )
        ymax = max(float(p_all.max()), float(p_sig.max()), float(p_smprot.max()))
    else:
        w = 0.38
        ax.bar(x - w / 2, 100 * p_all, width=w, label=f"Filtered (all FASTA, n={n_all})", color=pal.OI_SKY_BLUE, edgecolor=pal.OI_BLUE)
        ax.bar(x + w / 2, 100 * p_sig, width=w, label=f"Significant canonical genes (n={n_sig})", color=pal.OI_ORANGE, edgecolor=pal.OI_VERMILLION)
        ymax = max(float(p_all.max()), float(p_sig.max()))
    ax.set_xticks(x)
    ax.set_xticklabels(AA20)
    ax.set_xlabel("First amino acid (N-terminal)")
    ax.set_ylabel("Percent of peptides")
    ttl = "Starting amino acid distribution: filtered vs canonical"
    if c_smprot is not None:
        ttl += " vs SmProt significant (exportable)"
    ax.set_title(ttl)
    ax.legend(loc="upper right")
    ax.set_ylim(0, max(100 * ymax * 1.15, 5))
    ax.yaxis.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    plt.close(fig)

    lines = [
        "Starting AA distribution (TCGA-filtered peptide FASTA + optional significant export FASTA)",
        f"FASTA (filtered): {MP_FAA}",
        f"Peptides with standard AA at position 1 - all filtered: n={n_all} (skipped non-standard/empty: {skip_all})",
        f"Peptides from GeneSymbol in canonical list ({len(tr)} genes): n={n_sig} (skipped: {skip_sig})",
    ]
    if c_smprot is not None:
        lines.append(
            f"SmProt significant exportable ({SIG_FAA.name}): n={n_smprot} (skipped: {skip_smprot})"
        )
    lines += [
        "",
        f"Chi-square homogeneity (20x{table.shape[1]} contingency):",
        f"  chi2 = {chi2:.3f}, df = {dof}, p = {p_homog:.4g}",
        "  Interpretation: small p suggests first-AA frequencies differ across columns.",
        "",
        "Goodness-of-fit to uniform (each set separately):",
        f"  All filtered: chi2 = {chi2_uni_all:.3f}, p = {p_uni_all:.4g}",
        f"  Canonical subset: chi2 = {chi2_uni_sig:.3f}, p = {p_uni_sig:.4g}",
    ]
    if chi2_uni_smp is not None:
        lines.append(f"  SmProt significant exportable: chi2 = {chi2_uni_smp:.3f}, p = {p_uni_smp:.4g}")
    lines += [
        "",
        "Top 5 first AAs (filtered all):",
    ]
    top_all = df.nlargest(5, "n_filtered_all")[["aa", "n_filtered_all", "pct_filtered_all"]]
    for _, r in top_all.iterrows():
        lines.append(f"  {r['aa']}: n={int(r['n_filtered_all'])}, {r['pct_filtered_all']:.2f}%")
    lines.append("")
    lines.append("Top 5 first AAs (canonical subset):")
    top_sig = df.nlargest(5, "n_significant_canonical")[["aa", "n_significant_canonical", "pct_significant_canonical"]]
    for _, r in top_sig.iterrows():
        lines.append(f"  {r['aa']}: n={int(r['n_significant_canonical'])}, {r['pct_significant_canonical']:.2f}%")
    if c_smprot is not None and "n_significant_smprot_export" in df.columns:
        lines.append("")
        lines.append("Top 5 first AAs (SmProt significant exportable):")
        top_smp = df.nlargest(5, "n_significant_smprot_export")[
            ["aa", "n_significant_smprot_export", "pct_significant_smprot_export"]
        ]
        for _, r in top_smp.iterrows():
            lines.append(
                f"  {r['aa']}: n={int(r['n_significant_smprot_export'])}, {r['pct_significant_smprot_export']:.2f}%"
            )

    report = "\n".join(lines)
    OUT_TXT.write_text(report, encoding="utf-8")

    print(report)
    print(f"\nWrote: {OUT_FIG}\nWrote: {OUT_CSV}")


if __name__ == "__main__":
    main()
