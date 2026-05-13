"""
Bar charts: per cancer and (analysis × transition × direction), the percentage of
genes that **intersect limma and z for that exact stratum** (same cancer type and
transition) with **|z| ≥ 3**, have **≥1 SmProt peptide row** in a chosen peptide table (**default:** TCGA-matrix
``data/smprot_filtered_tcga_expr_genes.tsv``, matching ``smprot_tcga_filtered_peptides.faa``;
**alternate:** ``--peptide-gene-set all_smprot_filtered`` uses full ``data/smprot_filtered.tsv``),
out of all such intersected genes.

- **Stratum:** limma FDR < 0.05 (`limma_*_FDR0.05.csv`) for that cancer × transition
  × direction (**up** / **down**: limma ``logFC`` and z-detail ``log2FC`` both
  positive or both negative; **combined**: limma DE either direction, still merged
  with z rows for that cancer × transition with ``|z| ≥ 3``).

- **Z side:** ``tr_lncrna_output/tr_lncrnas_stage_detail.csv`` or
  ``tr_lncrnas_metastasis_detail.csv`` (same transitions as ``tr_lncrna_de_analysis.py``;
  rows are already z hits; we still enforce ``|z| ≥ 3``).

- **Vermillion dashed line — Overall Tr-lncRNAs:** fraction of **global** canonical
  Tr genes (``canonical_significant_lncRNAs.txt`` = limma gene union ∩ z union,
  ~1608) with ≥1 peptide from the same filtered table.

- **Bluish-green dashed line — Overall TCGA lncRNAs:** fraction of all matrix lncRNA
  columns with ≥1 peptide from the same table.

Outputs (under ``figures/<peptide_gene_set>/`` by default, see ``--figures-dir``):

  - ``peptide_fraction/*.png`` — one bar chart per analysis / transition / direction
  - ``tr_de_peptide_fraction_by_cancer.csv`` — long table for all panels

For ``--peptide-gene-set tcga_matrix``, the **stage E→L combined** panel is also copied to
``figures/fig2b_stage_E_L_combined.png`` (manuscript Fig 2B at repo root).
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

import figure_palettes as pal

DATA = ROOT / "data"
FIGURES = ROOT / "figures"
LIMMA = ROOT / "tr_lncrna_output" / "limma"
CANONICAL = DATA / "canonical_significant_lncRNAs.txt"
LIMMA_Z = LIMMA / "limma_z_intersection_genes.txt"
# Default peptide gene lists (override with --peptide-gene-set).
PEP_ALL_FILTERED_TSV = DATA / "smprot_filtered.tsv"
PEP_TCGA_MATRIX_TSV = DATA / "smprot_filtered_tcga_expr_genes.tsv"
STAGE_CSV = DATA / "primary_exp_stage_lncRNA.csv"
META_CSV = DATA / "primary_exp_metastasis_lncRNA.csv"
STAGE_LIMMA = LIMMA / "limma_stage_FDR0.05.csv"
META_LIMMA = LIMMA / "limma_metastasis_FDR0.05.csv"
TR_Z_STAGE = ROOT / "tr_lncrna_output" / "tr_lncrnas_stage_detail.csv"
TR_Z_META = ROOT / "tr_lncrna_output" / "tr_lncrnas_metastasis_detail.csv"

MIN_SAMPLES_CANCER = 100
FDR_MAX = 0.05
Z_ABS_MIN = 3.0

META_STAGE = ["sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage"]
META_META = ["sample_id", "cancer_type", "stage", "M_stage"]

TRANSITION_TITLES_STAGE: dict[str, str] = {
    "I_II": "Stage I → II",
    "II_III": "Stage II → III",
    "III_IV": "Stage III → IV",
    "E_L": "Early (I–II) → Late (III–IV)",
}
TRANSITION_TITLES_META: dict[str, str] = {
    "M0s_M0l": "M0 short → M0 long",
    "M0l_M1L": "M0 long → M1 lymphovascular",
    "ME_M1L": "M0 (short or long) → M1 lymphovascular",
}


def cancer_types_over_n(meta: pd.DataFrame, n: int = MIN_SAMPLES_CANCER) -> list[str]:
    counts = meta.groupby("cancer_type").size()
    return sorted(counts[counts > n].index.astype(str).tolist())


def matrix_gene_columns(csv_path: Path, meta_cols: list[str]) -> list[str]:
    head = pd.read_csv(csv_path, nrows=0)
    return [c for c in head.columns if c not in meta_cols]


def load_tr_genes() -> set[str]:
    path = CANONICAL if CANONICAL.exists() else LIMMA_Z
    if not path.exists():
        raise FileNotFoundError(
            f"Need {CANONICAL} or {LIMMA_Z}. Run tr_limma_de.R / z pipeline first."
        )
    return {ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()}


def load_peptide_genes_filtered(pep_tsv: Path) -> set[str]:
    if not pep_tsv.exists():
        raise FileNotFoundError(f"Missing {pep_tsv}. Run python build_significant_lncs_smprot.py.")
    pep = pd.read_csv(pep_tsv, sep="\t", usecols=["GeneSymbol"], low_memory=False)
    return set(pep["GeneSymbol"].astype(str).str.strip())


def load_limma_sig(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run tr_limma_de.R.")
    df = pd.read_csv(path, dtype={"gene": str, "cancer_type": str, "transition": str})
    df["cancer_type"] = df["cancer_type"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df["transition"] = df["transition"].astype(str).str.strip()
    return df.loc[df["adj.P.Val"].astype(float) < FDR_MAX].copy()


def filter_direction(df: pd.DataFrame, direction: str) -> pd.DataFrame:
    if direction == "combined":
        return df
    logfc = df["logFC"].astype(float)
    if direction == "up":
        return df.loc[logfc > 0]
    if direction == "down":
        return df.loc[logfc < 0]
    raise ValueError(direction)


def load_z_detail(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run python tr_lncrna_de_analysis.py.")
    df = pd.read_csv(
        path,
        dtype={"gene": str, "cancer_type": str, "transition": str, "analysis": str},
    )
    df["gene"] = df["gene"].astype(str).str.strip()
    df["cancer_type"] = df["cancer_type"].astype(str).str.strip()
    df["transition"] = df["transition"].astype(str).str.strip()
    df["analysis"] = df["analysis"].astype(str).str.strip()
    df["log2FC"] = df["log2FC"].astype(float)
    df["z"] = df["z"].astype(float)
    return df


def genes_limma_z_stratum(
    limma_rows: pd.DataFrame,
    z_frame: pd.DataFrame,
    *,
    analysis: str,
    transition: str,
    cancer_type: str,
    direction: str,
) -> set[str]:
    """Genes DE in limma for this stratum and present in z-detail with |z| >= 3."""
    L = limma_rows.loc[
        limma_rows["cancer_type"] == cancer_type,
        ["gene", "logFC"],
    ].drop_duplicates(subset=["gene"])
    if L.empty:
        return set()
    Z = z_frame.loc[
        (z_frame["analysis"] == analysis)
        & (z_frame["cancer_type"] == cancer_type)
        & (z_frame["transition"] == transition),
        ["gene", "log2FC", "z"],
    ].drop_duplicates(subset=["gene"])
    if Z.empty:
        return set()
    Z = Z.loc[np.isfinite(Z["z"]) & (Z["z"].abs() >= Z_ABS_MIN)]
    if Z.empty:
        return set()
    m = L.merge(Z, on="gene", how="inner")
    if m.empty:
        return set()
    if direction == "up":
        m = m.loc[(m["logFC"] > 0) & (m["log2FC"] > 0)]
    elif direction == "down":
        m = m.loc[(m["logFC"] < 0) & (m["log2FC"] < 0)]
    # combined: limma DE either direction, still require z hit with |z| >= 3
    return set(m["gene"].unique())


def pct_intersection(numer: set[str], denom: set[str]) -> float:
    if not denom:
        return float("nan")
    return 100.0 * len(numer & denom) / len(denom)


def plot_one(
    cancers: list[str],
    n_stratum: list[int],
    n_pep: list[int],
    title: str,
    ylabel: str,
    out_path: Path,
    overall_tr_canonical_pct: float,
    overall_tcga_matrix_pct: float,
) -> None:
    pcts = []
    for d, p in zip(n_stratum, n_pep, strict=True):
        if d > 0:
            pcts.append(100.0 * p / d)
        else:
            pcts.append(0.0)

    x = np.arange(len(cancers))
    fig_w = max(10.0, 0.42 * len(cancers))
    fig, ax = plt.subplots(figsize=(fig_w, 5.2))

    ax.bar(
        x,
        pcts,
        color=pal.OI_SKY_BLUE,
        edgecolor=pal.OI_BLUE,
        linewidth=1.0,
    )
    ax.set_xticks(x)
    ax.set_xticklabels([c.lower() for c in cancers], rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold")
    hi = max(pcts) if pcts else 0.0
    ax.set_ylim(0, max(12.0, hi * 1.15 + 0.01, 1.0))
    ax.yaxis.grid(True, linestyle="-", alpha=0.35)
    ax.set_axisbelow(True)

    for i, (d, p, pct) in enumerate(zip(n_stratum, n_pep, pcts, strict=True)):
        if d > 0:
            lab = f"{pct:.1f}%\n({p}/{d})"
        else:
            lab = "0.0%\n(0/0)"
        ax.text(
            i,
            min(pct + 0.35, ax.get_ylim()[1] * 0.98),
            lab,
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )

    if np.isfinite(overall_tr_canonical_pct):
        ax.axhline(
            overall_tr_canonical_pct,
            color=pal.OI_VERMILLION,
            linestyle="--",
            linewidth=1.2,
            label=f"Overall Tr-lncRNAs: {overall_tr_canonical_pct:.1f}%",
        )
    if np.isfinite(overall_tcga_matrix_pct):
        ax.axhline(
            overall_tcga_matrix_pct,
            color=pal.OI_BLUISH_GREEN,
            linestyle="--",
            linewidth=1.2,
            label=f"Overall TCGA lncRNAs: {overall_tcga_matrix_pct:.1f}%",
        )

    ax.legend(loc="upper right", frameon=True, edgecolor="0.7")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Peptide-fraction bars: limma DE intersect high-|z| genes with >=1 peptide row in a SmProt TSV."
        )
    )
    ap.add_argument(
        "--peptide-gene-set",
        choices=("tcga_matrix", "all_smprot_filtered"),
        default="tcga_matrix",
        help=(
            "Which SmProt TSV defines 'has peptide': tcga_matrix = data/smprot_filtered_tcga_expr_genes.tsv "
            "(default; genes in TCGA lncRNA matrices; matches smprot_tcga_filtered_peptides.faa). "
            "all_smprot_filtered = data/smprot_filtered.tsv (full curated SmProt list)."
        ),
    )
    ap.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help="Parent directory for outputs (default: UNDEFINED/figures). PNGs go in <parent>/<peptide_gene_set>/peptide_fraction/.",
    )
    args = ap.parse_args()

    figures_parent = args.figures_dir if args.figures_dir is not None else FIGURES
    suite_dir = figures_parent / args.peptide_gene_set
    figs = suite_dir / "peptide_fraction"
    out_csv_name = "tr_de_peptide_fraction_by_cancer.csv"

    if args.peptide_gene_set == "tcga_matrix":
        pep_tsv = PEP_TCGA_MATRIX_TSV
    else:
        pep_tsv = PEP_ALL_FILTERED_TSV

    suite_dir.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    tr_set = load_tr_genes()
    pep_set = load_peptide_genes_filtered(pep_tsv)
    genes_stage = set(matrix_gene_columns(STAGE_CSV, META_STAGE))
    genes_meta = set(matrix_gene_columns(META_CSV, META_META))
    if genes_stage != genes_meta:
        print("Note: stage vs metastasis matrix gene columns differ; using stage set for TCGA baseline.")
    all_lnc = genes_stage

    overall_tr_canonical = pct_intersection(pep_set, tr_set)
    overall_tcga_matrix = pct_intersection(pep_set, all_lnc)
    n_tr_pep = len(tr_set & pep_set)
    n_mtx_pep = len(all_lnc & pep_set)
    print(
        f"Canonical Tr-lncRNAs: {len(tr_set)}; with >=1 peptide ({pep_tsv.name}): "
        f"{n_tr_pep} ({overall_tr_canonical:.2f}%)"
    )
    print(
        f"TCGA matrix lncRNA columns: {len(all_lnc)}; with >=1 peptide: "
        f"{n_mtx_pep} ({overall_tcga_matrix:.2f}%)"
    )

    meta_s = pd.read_csv(STAGE_CSV, usecols=["cancer_type"], low_memory=False)
    meta_m = pd.read_csv(META_CSV, usecols=["cancer_type", "M_stage"], low_memory=False)
    meta_m = meta_m.loc[meta_m["cancer_type"].notna() & (meta_m["M_stage"].astype(str) != "M1_s")].copy()
    cancers_stage = cancer_types_over_n(meta_s)
    cancers_meta = cancer_types_over_n(meta_m)

    lim_s = load_limma_sig(STAGE_LIMMA)
    lim_m = load_limma_sig(META_LIMMA)

    summary_rows: list[dict] = []

    z_stage = load_z_detail(TR_Z_STAGE)
    z_meta = load_z_detail(TR_Z_META)

    blocks: list[tuple[str, pd.DataFrame, pd.DataFrame, list[str], dict[str, str]]] = [
        ("stage", lim_s, z_stage, cancers_stage, TRANSITION_TITLES_STAGE),
        ("metastasis", lim_m, z_meta, cancers_meta, TRANSITION_TITLES_META),
    ]

    for analysis, lim_all, z_all, cancers_order, titles in blocks:
        for transition in sorted(lim_all["transition"].unique()):
            base = lim_all.loc[
                (lim_all["analysis"] == analysis) & (lim_all["transition"] == transition)
            ]
            if base.empty:
                continue
            tr_title = titles.get(transition, transition)
            for direction, dir_label in (
                ("up", "Up (limma & z, |z|≥3)"),
                ("down", "Down (limma & z, |z|≥3)"),
                ("combined", "Combined"),
            ):
                lim_rows = base if direction == "combined" else filter_direction(base, direction)

                n_stratum_list: list[int] = []
                n_pep_list: list[int] = []
                for ct in cancers_order:
                    gset = genes_limma_z_stratum(
                        lim_rows,
                        z_all,
                        analysis=analysis,
                        transition=transition,
                        cancer_type=ct,
                        direction=direction,
                    )
                    n = len(gset)
                    w = len(gset & pep_set)
                    n_stratum_list.append(n)
                    n_pep_list.append(w)

                for ct, n, w in zip(cancers_order, n_stratum_list, n_pep_list, strict=True):
                    pct = 100.0 * w / n if n else float("nan")
                    summary_rows.append(
                        {
                            "analysis": analysis,
                            "transition": transition,
                            "direction": direction,
                            "cancer_type": ct,
                            "n_limma_z_stratum": n,
                            "n_limma_z_with_peptide": w,
                            "pct_limma_z_stratum_with_peptide": pct,
                            "overall_tr_canonical_pct": overall_tr_canonical,
                            "n_canonical_tr": len(tr_set),
                            "n_canonical_tr_with_peptide": n_tr_pep,
                            "overall_tcga_matrix_lnc_pct": overall_tcga_matrix,
                            "n_matrix_lnc_genes": len(all_lnc),
                            "n_matrix_lnc_with_peptide": n_mtx_pep,
                            "peptide_gene_set": args.peptide_gene_set,
                            "peptide_genes_tsv": str(pep_tsv),
                        }
                    )
                title = (
                    "lncRNAs with filtered SmProt peptides\n"
                    f"(% of limma DE ∩ |z|≥3 per cancer; {tr_title}; {dir_label})"
                )
                fname = f"{analysis}_{transition}_{direction}.png"
                out_panel = figs / fname
                plot_one(
                    cancers_order,
                    n_stratum_list,
                    n_pep_list,
                    title,
                    "Percentage of peptide-producing lncRNAs",
                    out_panel,
                    overall_tr_canonical,
                    overall_tcga_matrix,
                )
                sum_stratum = sum(n_stratum_list)
                print(f"Wrote {out_panel}  (sum limma-z stratum genes across cancers: {sum_stratum})")
                # Fig 2B (manuscript): stage Early→Late combined panel also at repo figures/ root.
                if (
                    args.peptide_gene_set == "tcga_matrix"
                    and analysis == "stage"
                    and transition == "E_L"
                    and direction == "combined"
                ):
                    root_fig2b = FIGURES / "fig2b_stage_E_L_combined.png"
                    FIGURES.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(out_panel, root_fig2b)
                    print(f"Copied to {root_fig2b} (Fig 2B)")

    summ = pd.DataFrame(summary_rows)
    if not summ.empty:
        out_csv = suite_dir / out_csv_name
        summ.to_csv(out_csv, index=False)
        print(f"Wrote {out_csv} ({len(summ)} rows)")


if __name__ == "__main__":
    main()
