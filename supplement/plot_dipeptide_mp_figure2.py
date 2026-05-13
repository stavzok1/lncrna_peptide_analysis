"""
Fig. 2 style: 20x20 dipeptide (overlapping) composition vs human reference proteome.

Inputs:
  - data/known_proteins.fasta  (reference proteome)
  - data/smprot_tcga_filtered_peptides.faa  (filtered MPs; TCGA lncRNA panel;
    headers >smPEP_ID|GeneSymbol|...)
  - data/significant_lnc_peptides.faa  (optional): SmProt significant peptides that
    have an exportable translated sequence (same rows as ``significant_lnc_peptides.tsv``).

Panels:
  (A) log2 ratio: all lncRNA-MPs vs proteome
  (B) log2 ratio: Tr-lncRNA-MPs (canonical gene list) vs proteome
  (C–D) Cohen's h (two-proportion effect): h = 2*arcsin(sqrt(p_mp)) - 2*arcsin(sqrt(p_ref))
        with smoothed p = (c+0.5)/(N+200) per cell (Cohen 1988). Not a Wald z-score.
  When the significant-peptide FASTA exists, two extra rows compare that set to the proteome
  (log2 ratio and Cohen's h).
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
from Bio import SeqIO

import figure_palettes as pal

DATA = ROOT / "data"
PROTEOME_FA = DATA / "known_proteins.fasta"
MP_FAA = DATA / "smprot_tcga_filtered_peptides.faa"
SIG_FAA = DATA / "significant_lnc_peptides.faa"
CANONICAL = DATA / "canonical_significant_lncRNAs.txt"
LIMMA_Z = ROOT / "tr_lncrna_output" / "limma" / "limma_z_intersection_genes.txt"
OUT_DIR = ROOT / "tr_lncrna_output" / "figures"
OUT_FIG = OUT_DIR / "fig2_dipeptide_mp_composition.png"

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA20)}


def load_tr_genes() -> set[str]:
    path = CANONICAL if CANONICAL.exists() else LIMMA_Z
    if not path.exists():
        raise FileNotFoundError(f"Need {CANONICAL} or {LIMMA_Z}")
    return {ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()}


def count_dipeptides_from_string(seq: str, mat: np.ndarray) -> tuple[int, int]:
    """Add overlapping dipeptide counts from one protein string. Returns (aa_used, dipep_used)."""
    s = seq.upper().replace("U", "C")
    aa_n = di = 0
    for ch in s:
        if ch in AA_INDEX:
            aa_n += 1
    for i in range(len(s) - 1):
        a, b = s[i], s[i + 1]
        if a in AA_INDEX and b in AA_INDEX:
            mat[AA_INDEX[a], AA_INDEX[b]] += 1.0
            di += 1
    return aa_n, di


def count_proteome(path: Path) -> tuple[np.ndarray, int, int]:
    mat = np.zeros((20, 20), dtype=np.float64)
    total_aa = total_di = 0
    for rec in SeqIO.parse(path, "fasta"):
        aa_n, di = count_dipeptides_from_string(str(rec.seq), mat)
        total_aa += aa_n
        total_di += di
    return mat, total_aa, total_di


def parse_mp_fasta(path: Path) -> list[tuple[str, str]]:
    """List of (GeneSymbol, aa_sequence) from MP FASTA headers (>id|GeneSymbol|...)."""
    out: list[tuple[str, str]] = []
    for rec in SeqIO.parse(path, "fasta"):
        parts = rec.description.split("|")
        gene = parts[1].strip() if len(parts) > 1 else ""
        out.append((gene, str(rec.seq).upper()))
    return out


def count_mp_subset(records: list[tuple[str, str]], gene_filter: set[str] | None) -> tuple[np.ndarray, int, int]:
    mat = np.zeros((20, 20), dtype=np.float64)
    total_aa = total_di = 0
    for gene, seq in records:
        if gene_filter is not None and gene not in gene_filter:
            continue
        aa_n, di = count_dipeptides_from_string(seq, mat)
        total_aa += aa_n
        total_di += di
    return mat, total_aa, total_di


def freq_matrix(counts: np.ndarray, pseudo: float = 0.5) -> np.ndarray:
    """Multinomial-style smoothing over 400 cells."""
    c = counts + pseudo
    return c / c.sum()


def log2_ratio_mat(f_mp: np.ndarray, f_ref: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    return np.log2(np.maximum(f_mp, eps) / np.maximum(f_ref, eps))


def cohens_h_mat(c_mp: np.ndarray, di_mp: float, c_ref: np.ndarray, di_ref: float) -> np.ndarray:
    """
    Cohen's h for two independent proportions (per dipeptide cell).

    Uses Laplace-style smoothing (c+0.5)/(N+200) so each of 400 cells has the same
    prior mass as in freq_matrix; keeps h finite when c=0.
    """
    denom_mp = max(di_mp, 1.0) + 200.0
    denom_ref = max(di_ref, 1.0) + 200.0
    p_mp = (c_mp + 0.5) / denom_mp
    p_ref = (c_ref + 0.5) / denom_ref
    p_mp = np.clip(p_mp, 0.0, 1.0)
    p_ref = np.clip(p_ref, 0.0, 1.0)
    phi_mp = 2.0 * np.arcsin(np.sqrt(p_mp))
    phi_ref = 2.0 * np.arcsin(np.sqrt(p_ref))
    return phi_mp - phi_ref


def plot_heatmap(
    ax,
    data: np.ndarray,
    title: str,
    cbar_label: str,
    *,
    kind: str,
) -> None:
    if kind == "log2":
        vmax = float(np.nanpercentile(np.abs(data), 99))
        vmax = max(vmax, 0.5)
        im = ax.imshow(data, cmap=pal.diverging_log2fc_cmap(), aspect="auto", vmin=-vmax, vmax=vmax)
    elif kind == "cohen_h":
        vmax = float(np.nanpercentile(np.abs(data), 99))
        vmax = max(vmax, 0.02)
        im = ax.imshow(data, cmap=pal.diverging_log2fc_cmap(), aspect="auto", vmin=-vmax, vmax=vmax)
    else:
        raise ValueError(kind)
    ax.set_xticks(np.arange(20))
    ax.set_yticks(np.arange(20))
    ax.set_xticklabels(list(AA20), fontsize=8)
    ax.set_yticklabels(list(AA20), fontsize=8)
    ax.set_xlabel("Second amino acid")
    ax.set_ylabel("First amino acid")
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbar_label)


def main() -> None:
    if not PROTEOME_FA.exists():
        raise FileNotFoundError(f"Missing {PROTEOME_FA}")
    if not MP_FAA.exists():
        raise FileNotFoundError(f"Missing {MP_FAA}. Run export_tcga_filtered_peptides_fasta.py")

    tr_genes = load_tr_genes()
    mp_recs = parse_mp_fasta(MP_FAA)

    print("Counting proteome dipeptides (may take a minute)...")
    c_ref, aa_ref, di_ref = count_proteome(PROTEOME_FA)
    f_ref = freq_matrix(c_ref)

    c_all, aa_all, di_all = count_mp_subset(mp_recs, gene_filter=None)
    c_tr, aa_tr, di_tr = count_mp_subset(mp_recs, gene_filter=tr_genes)

    f_all = freq_matrix(c_all)
    f_tr = freq_matrix(c_tr)

    log_all = log2_ratio_mat(f_all, f_ref)
    log_tr = log2_ratio_mat(f_tr, f_ref)
    h_all = cohens_h_mat(c_all, float(di_all), c_ref, float(di_ref))
    h_tr = cohens_h_mat(c_tr, float(di_tr), c_ref, float(di_ref))

    sig_recs: list[tuple[str, str]] | None = None
    log_smprot: np.ndarray | None = None
    h_smprot: np.ndarray | None = None
    aa_smprot = di_smprot = 0
    if SIG_FAA.exists():
        sig_recs = parse_mp_fasta(SIG_FAA)
        c_smprot, aa_smprot, di_smprot = count_mp_subset(sig_recs, gene_filter=None)
        f_smprot = freq_matrix(c_smprot)
        log_smprot = log2_ratio_mat(f_smprot, f_ref)
        h_smprot = cohens_h_mat(c_smprot, float(di_smprot), c_ref, float(di_ref))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nrows = 3 if sig_recs else 2
    fig, axes = plt.subplots(nrows, 2, figsize=(11, 5.0 * nrows))

    plot_heatmap(
        axes[0, 0],
        log_all,
        f"(A) Log2 ratio: lncRNA-MPs / proteome\n({aa_all:,} aa MPs; {aa_ref:,} aa proteome)",
        "Log2 ratio",
        kind="log2",
    )
    plot_heatmap(
        axes[0, 1],
        log_tr,
        f"(B) Log2 ratio: Tr-lncRNA-MPs / proteome\n({aa_tr:,} aa; {aa_ref:,} aa proteome)",
        "Log2 ratio",
        kind="log2",
    )
    plot_heatmap(
        axes[1, 0],
        h_all,
        "(C) Cohen's h: lncRNA-MPs vs proteome\n(arcsin sqrt proportions)",
        "Cohen's h",
        kind="cohen_h",
    )
    plot_heatmap(
        axes[1, 1],
        h_tr,
        "(D) Cohen's h: Tr-lncRNA-MPs vs proteome\n(arcsin sqrt proportions)",
        "Cohen's h",
        kind="cohen_h",
    )
    if sig_recs is not None and log_smprot is not None and h_smprot is not None:
        n_pep = len(sig_recs)
        plot_heatmap(
            axes[2, 0],
            log_smprot,
            f"(E) Log2 ratio: significant SmProt MPs / proteome\n({aa_smprot:,} aa; {n_pep} peptides; {aa_ref:,} aa proteome)",
            "Log2 ratio",
            kind="log2",
        )
        plot_heatmap(
            axes[2, 1],
            h_smprot,
            "(F) Cohen's h: significant SmProt MPs vs proteome\n(arcsin sqrt proportions)",
            "Cohen's h",
            kind="cohen_h",
        )

    extra = "\n(E–F) significant SmProt peptides (exportable FASTA)." if sig_recs else ""
    fig.suptitle(
        "Fig. 2. Dipeptide aa composition (TCGA-filtered MPs vs canonical Tr subset vs proteome)"
        + extra,
        fontsize=12,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)

    summary = (
        f"Proteome: {aa_ref:,} aa counted (standard AA), {di_ref:,} overlapping dipeptides.\n"
        f"lncRNA-MPs (all filtered FASTA): {aa_all:,} aa, {di_all:,} dipeptides.\n"
        f"Tr-lncRNA-MPs (GeneSymbol in canonical, {len(tr_genes)} genes): {aa_tr:,} aa, {di_tr:,} dipeptides.\n"
    )
    if sig_recs is not None:
        summary += (
            f"Significant SmProt MPs (exportable FASTA): {aa_smprot:,} aa, {di_smprot:,} dipeptides "
            f"({len(sig_recs)} parent sequences).\n"
        )
    summary += f"Saved: {OUT_FIG}"
    print(summary)
    (OUT_DIR / "fig2_dipeptide_summary.txt").write_text(summary, encoding="utf-8")


if __name__ == "__main__":
    main()
