"""
Compare old lists vs canonical lncRNAs (limma∩z) and current SmProt curation outputs.

Filtered-peptide comparison uses only rows whose GeneSymbol is a TCGA lncRNA
expression column (union of stage and metastasis matrices), with smPEP_ID deduped.
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
import json

import pandas as pd

DATA = ROOT / "data"
OUT = DATA / "compare_old_new"
OUT.mkdir(exist_ok=True)

OLD_LNCS = DATA / "significant_lncs_list_old.txt"
OLD_SIG_PEPS = DATA / "significant_peps_old.csv"
OLD_FILT = DATA / "filtered_peps_old.csv"

NEW_LNCS = DATA / "canonical_significant_lncRNAs.txt"
NEW_SIG_PEPS = DATA / "significant_lnc_peptides.tsv"
NEW_FILT_TCGA = DATA / "smprot_filtered_tcga_expr_genes.tsv"
NEW_FILT_ALL = DATA / "smprot_filtered.tsv"

STAGE_EXPR = DATA / "primary_exp_stage_lncRNA.csv"
META_EXPR = DATA / "primary_exp_metastasis_lncRNA.csv"
META_STAGE_COLS = {"sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage"}
META_META_COLS = {"sample_id", "cancer_type", "stage", "M_stage"}


def norm_pep_id(s) -> str:
    if pd.isna(s):
        return ""
    if isinstance(s, float) and s == int(s):
        return str(int(s))
    return str(s).strip()


def load_gene_set(path: Path) -> set[str]:
    if path.suffix == ".csv":
        df = pd.read_csv(path, usecols=["gene_name"])
        return set(df["gene_name"].astype(str).str.strip())
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return {ln.strip() for ln in lines if ln.strip()}


def load_smpep_set(path: Path, sep: str | None = None) -> set[str]:
    kwargs = {"dtype": str} if sep == "\t" else {}
    df = pd.read_csv(path, sep=sep, usecols=["smPEP_ID"], **kwargs)
    return {norm_pep_id(x) for x in df["smPEP_ID"]}


def tcga_lncrna_expr_symbols() -> set[str]:
    """Gene symbols that are expression columns in the TCGA lncRNA matrices (union)."""
    c1 = pd.read_csv(STAGE_EXPR, nrows=0).columns
    c2 = pd.read_csv(META_EXPR, nrows=0).columns
    g1 = {str(c) for c in c1 if c not in META_STAGE_COLS}
    g2 = {str(c) for c in c2 if c not in META_META_COLS}
    return g1 | g2


def load_filtered_smpep_tcga_only(path: Path, sep: str, tcga_genes: set[str]) -> set[str]:
    """smPEP_ID set after restricting rows to GeneSymbol in tcga_genes and deduping smPEP_ID."""
    df = pd.read_csv(path, sep=sep, usecols=["smPEP_ID", "GeneSymbol"], dtype=str, low_memory=False)
    sym = df["GeneSymbol"].astype(str)
    df = df.loc[sym.isin(tcga_genes)].drop_duplicates(subset=["smPEP_ID"], keep="first")
    return {norm_pep_id(x) for x in df["smPEP_ID"]}


def compare_sets(name: str, old: set[str], new: set[str]) -> dict:
    inter = old & new
    old_only = old - new
    new_only = new - old
    rep = {
        "comparison": name,
        "n_old": len(old),
        "n_new": len(new),
        "n_intersection": len(inter),
        "n_old_only": len(old_only),
        "n_new_only": len(new_only),
        "jaccard": len(inter) / len(old | new) if (old | new) else 0.0,
    }
    (OUT / f"{name}_intersection.txt").write_text("\n".join(sorted(inter)), encoding="utf-8")
    (OUT / f"{name}_old_only.txt").write_text("\n".join(sorted(old_only)), encoding="utf-8")
    (OUT / f"{name}_new_only.txt").write_text("\n".join(sorted(new_only)), encoding="utf-8")
    return rep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    old_genes = load_gene_set(OLD_LNCS)
    new_genes = load_gene_set(NEW_LNCS)

    old_sig = load_smpep_set(OLD_SIG_PEPS, sep=",")
    new_sig = load_smpep_set(NEW_SIG_PEPS, sep="\t")

    tcga_syms = tcga_lncrna_expr_symbols()
    old_filt = load_filtered_smpep_tcga_only(OLD_FILT, ",", tcga_syms)
    if NEW_FILT_TCGA.exists():
        new_filt = load_smpep_set(NEW_FILT_TCGA, sep="\t")
    else:
        new_filt = load_filtered_smpep_tcga_only(NEW_FILT_ALL, "\t", tcga_syms)

    reports = [
        compare_sets("significant_lnc_genes_canonical_vs_old", old_genes, new_genes),
        compare_sets("significant_peptides_canonical_vs_old", old_sig, new_sig),
        compare_sets("filtered_peptides_tcga_lncRNA_smPEP_ID", old_filt, new_filt),
    ]

    summary_path = OUT / "comparison_summary.json"
    summary_path.write_text(json.dumps(reports, indent=2), encoding="utf-8")

    print(json.dumps(reports, indent=2))
    print(f"\nWrote lists and {summary_path}")


if __name__ == "__main__":
    main()
