"""
Figure 3C–3D: 20×20 dipeptide log2 fold-change heatmaps vs ``data/known_proteins.fasta``.

Uses the same counting and smoothing as ``plot_dipeptide_mp_figure2.py`` (overlapping
2-mers, standard 20 aa, U→C, ``freq_matrix`` with +0.5 over 400 cells, ``log2(mp/proteome)``).

- **Figure 3C:** **two separate PNG files** when both FASTAs exist — (1) TCGA-matrix
  filtered MPs vs proteome; (2) all SmProt-filtered MPs vs proteome. If the all-filtered
  FASTA is missing, only the TCGA-matrix 3C file is written.

- **Figure 3D:** **Tr-lncRNA** MPs — same TCGA-matrix FASTA restricted to ``GeneSymbol``
  in ``data/canonical_significant_lncRNAs.txt`` (fallback: ``limma_z_intersection_genes.txt``)
  vs proteome.

Default output directory: ``figures/`` at the repository root.
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

_SUP = ROOT / "supplement"
if str(_SUP) not in sys.path:
    sys.path.insert(0, str(_SUP))


import argparse

import matplotlib.pyplot as plt
import numpy as np

import plot_dipeptide_mp_figure2 as dp

import figure_palettes as pal

DATA = ROOT / "data"
TCGA_FAA = DATA / "smprot_tcga_filtered_peptides.faa"
ALL_FILTERED_FAA = DATA / "smprot_all_filtered_peptides.faa"
PROTEOME_FA = DATA / "known_proteins.fasta"
FIGURES = ROOT / "figures"


def plot_log2_panel(
    ax,
    log2_mat: np.ndarray,
    title: str,
    subtitle: str,
    *,
    vmax: float,
) -> None:
    im = ax.imshow(log2_mat, cmap=pal.diverging_log2fc_cmap(), aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(20))
    ax.set_yticks(np.arange(20))
    ax.set_xticklabels(list(dp.AA20), fontsize=8)
    ax.set_yticklabels(list(dp.AA20), fontsize=8)
    ax.set_xlabel("Second amino acid")
    ax.set_ylabel("First amino acid")
    ax.set_title(f"{title}\n{subtitle}", fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Log2 ratio (peptide / proteome)")


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 3C–3D dipeptide log2FC heatmaps vs proteome.")
    ap.add_argument("--tcga-fa", type=Path, default=TCGA_FAA)
    ap.add_argument("--all-filtered-fa", type=Path, default=ALL_FILTERED_FAA)
    ap.add_argument("--proteome-fa", type=Path, default=PROTEOME_FA)
    ap.add_argument("--out-dir", type=Path, default=FIGURES, help="Directory for PNG outputs.")
    args = ap.parse_args()

    if not args.proteome_fa.exists():
        raise SystemExit(f"Missing proteome FASTA: {args.proteome_fa}")
    if not args.tcga_fa.exists():
        raise SystemExit(f"Missing TCGA-matrix peptide FASTA: {args.tcga_fa}")

    print("Counting proteome dipeptides (may take a minute)...")
    c_ref, aa_ref, di_ref = dp.count_proteome(args.proteome_fa)
    f_ref = dp.freq_matrix(c_ref)

    tr_genes = dp.load_tr_genes()
    rec_tcga = dp.parse_mp_fasta(args.tcga_fa)
    c_tcga_all, aa_tcga, di_tcga = dp.count_mp_subset(rec_tcga, gene_filter=None)
    c_tcga_tr, aa_tr, di_tr = dp.count_mp_subset(rec_tcga, gene_filter=tr_genes)
    log_tcga = dp.log2_ratio_mat(dp.freq_matrix(c_tcga_all), f_ref)
    log_tr = dp.log2_ratio_mat(dp.freq_matrix(c_tcga_tr), f_ref)

    all_filtered_exists = args.all_filtered_fa.exists()
    log_all_filtered: np.ndarray | None = None
    aa_all = di_all = 0
    if all_filtered_exists:
        rec_all = dp.parse_mp_fasta(args.all_filtered_fa)
        c_allf, aa_all, di_all = dp.count_mp_subset(rec_all, gene_filter=None)
        log_all_filtered = dp.log2_ratio_mat(dp.freq_matrix(c_allf), f_ref)

    def vmax_for(mat: np.ndarray) -> float:
        return max(float(np.nanpercentile(np.abs(mat), 99)), 0.5)

    vmax_tcga_3c = vmax_for(log_tcga)
    vmax_all_3c = vmax_for(log_all_filtered) if log_all_filtered is not None else vmax_tcga_3c
    vmax_3d = max(float(np.nanpercentile(np.abs(log_tr), 99)), 0.5)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- Figure 3C (separate files: TCGA-matrix vs all-filtered) ---
    fig_tcga, ax_tcga = plt.subplots(figsize=(6.5, 5.2), dpi=150)
    plot_log2_panel(
        ax_tcga,
        log_tcga,
        "TCGA-matrix filtered lncRNA MPs vs proteome",
        f"{aa_tcga:,} aa in MPs; {aa_ref:,} aa proteome",
        vmax=vmax_tcga_3c,
    )
    fig_tcga.tight_layout()
    out_tcga = args.out_dir / "fig3c_dipeptide_log2fc_tcga_matrix_vs_proteome.png"
    fig_tcga.savefig(out_tcga, bbox_inches="tight")
    plt.close(fig_tcga)
    print(f"Wrote {out_tcga}")

    if log_all_filtered is not None:
        fig_all, ax_all = plt.subplots(figsize=(6.5, 5.2), dpi=150)
        plot_log2_panel(
            ax_all,
            log_all_filtered,
            "All SmProt-filtered lncRNA MPs vs proteome",
            f"{aa_all:,} aa in MPs; {aa_ref:,} aa proteome",
            vmax=vmax_all_3c,
        )
        fig_all.tight_layout()
        out_all = args.out_dir / "fig3c_dipeptide_log2fc_all_smprot_filtered_vs_proteome.png"
        fig_all.savefig(out_all, bbox_inches="tight")
        plt.close(fig_all)
        print(f"Wrote {out_all}")
    else:
        print(f"Skip 3C all-filtered: missing {args.all_filtered_fa}")

    # --- Figure 3D ---
    fig_d, ax_d = plt.subplots(figsize=(6.5, 5.2), dpi=150)
    plot_log2_panel(
        ax_d,
        log_tr,
        "(D) Tr-lncRNA MPs (canonical genes, TCGA-matrix FASTA) vs proteome",
        f"{aa_tr:,} aa in MPs; {aa_ref:,} aa proteome; {len(tr_genes)} canonical genes",
        vmax=vmax_3d,
    )
    fig_d.tight_layout()
    out_d = args.out_dir / "fig3d_dipeptide_log2fc_tr_lncrna_tcga_vs_proteome.png"
    fig_d.savefig(out_d, bbox_inches="tight")
    plt.close(fig_d)
    print(f"Wrote {out_d}")


if __name__ == "__main__":
    main()
