"""
1) Peptide counts from smprot_filtered.tsv (>=10 aa already applied in that file;
   no TCGA expression-column filter).

   Writes:
     - data/smprot_filtered_no_tcga_summary.csv  (aggregate counts)
     - data/smprot_filtered_no_tcga_peptides_by_gene.csv  (per GeneSymbol)

2) Genes (TCGA lncRNA expression columns) tied to significant peptides:
   - New: significant_lnc_peptides.tsv ∩ TCGA gene symbols
   - Old: significant_peps_old.csv ∩ TCGA gene symbols

   Writes:
     - data/compare_old_new/significant_peptide_genes_tcga_old_new.csv
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



import pandas as pd

DATA = ROOT / "data"
OUT_COMPARE = DATA / "compare_old_new"

SMPROT_FILT = DATA / "smprot_filtered.tsv"
SIG_NEW = DATA / "significant_lnc_peptides.tsv"
SIG_OLD = DATA / "significant_peps_old.csv"
STAGE_EXPR = DATA / "primary_exp_stage_lncRNA.csv"
META_EXPR = DATA / "primary_exp_metastasis_lncRNA.csv"
META_STAGE_COLS = {"sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage"}
META_META_COLS = {"sample_id", "cancer_type", "stage", "M_stage"}

SUMMARY_CSV = DATA / "smprot_filtered_no_tcga_summary.csv"
BY_GENE_CSV = DATA / "smprot_filtered_no_tcga_peptides_by_gene.csv"
GENE_COMPARE_CSV = OUT_COMPARE / "significant_peptide_genes_tcga_old_new.csv"


def tcga_lncrna_expr_symbols() -> set[str]:
    c1 = pd.read_csv(STAGE_EXPR, nrows=0).columns
    c2 = pd.read_csv(META_EXPR, nrows=0).columns
    g1 = {str(c) for c in c1 if c not in META_STAGE_COLS}
    g2 = {str(c) for c in c2 if c not in META_META_COLS}
    return g1 | g2


def main() -> None:
    if not SMPROT_FILT.exists():
        raise FileNotFoundError(f"Missing {SMPROT_FILT}. Run build_significant_lncs_smprot.py")

    df = pd.read_csv(
        SMPROT_FILT,
        sep="\t",
        usecols=["smPEP_ID", "GeneSymbol"],
        dtype=str,
        low_memory=False,
    )
    df["GeneSymbol"] = df["GeneSymbol"].astype(str).str.strip()
    df["smPEP_ID"] = df["smPEP_ID"].astype(str).str.strip()
    df = df.drop_duplicates(subset=["smPEP_ID"], keep="first")

    n_pep = len(df)
    n_genes = df["GeneSymbol"].nunique()

    pd.DataFrame(
        [
            {"metric": "unique_smPEP_ID", "value": n_pep, "note": "dedup smPEP_ID; >=10 aa in build pipeline"},
            {"metric": "unique_GeneSymbol", "value": n_genes, "note": "genes with >=1 filtered peptide; no TCGA column filter"},
            {"metric": "source_file", "value": str(SMPROT_FILT), "note": ""},
        ]
    ).to_csv(SUMMARY_CSV, index=False)

    by_gene = (
        df.groupby("GeneSymbol", as_index=False)
        .agg(n_unique_smPEP_ID=("smPEP_ID", "count"))
        .sort_values("n_unique_smPEP_ID", ascending=False)
    )
    by_gene.to_csv(BY_GENE_CSV, index=False)

    print(f"Wrote {SUMMARY_CSV} (peptides={n_pep}, genes={n_genes})")
    print(f"Wrote {BY_GENE_CSV} ({len(by_gene)} genes)")

    tcga = tcga_lncrna_expr_symbols()

    rows_compare: list[dict] = []
    genes_new: set[str] = set()
    genes_old: set[str] = set()

    if SIG_NEW.exists():
        sn = pd.read_csv(SIG_NEW, sep="\t", usecols=["smPEP_ID", "GeneSymbol"], dtype=str, low_memory=False)
        sn["GeneSymbol"] = sn["GeneSymbol"].astype(str).str.strip()
        sn = sn.drop_duplicates(subset=["smPEP_ID"], keep="first")
        rows_compare.append(
            {
                "set_label": "new_significant_lnc_peptides_all_genes",
                "n_unique_smPEP_ID": int(sn["smPEP_ID"].nunique()),
                "n_unique_GeneSymbol": int(sn["GeneSymbol"].nunique()),
                "n_intersection_GeneSymbol": "",
                "n_new_only_GeneSymbol": "",
                "n_old_only_GeneSymbol": "",
                "source": str(SIG_NEW),
                "filter": "no TCGA expression-column filter (matches compare_sig old-vs-new peptide count)",
            }
        )
        sn_tcga = sn.loc[sn["GeneSymbol"].isin(tcga)]
        genes_new = set(sn_tcga["GeneSymbol"].astype(str).unique())
        peps_new = sn_tcga["smPEP_ID"].nunique()
        rows_compare.append(
            {
                "set_label": "new_significant_lnc_peptides_tcga_genes_only",
                "n_unique_smPEP_ID": int(peps_new),
                "n_unique_GeneSymbol": len(genes_new),
                "n_intersection_GeneSymbol": "",
                "n_new_only_GeneSymbol": "",
                "n_old_only_GeneSymbol": "",
                "source": str(SIG_NEW),
                "filter": "GeneSymbol in TCGA lncRNA expression matrix columns",
            }
        )
    else:
        rows_compare.append(
            {
                "set_label": "new_significant_lnc_peptides_missing",
                "n_unique_smPEP_ID": "",
                "n_unique_GeneSymbol": "",
                "n_intersection_GeneSymbol": "",
                "n_new_only_GeneSymbol": "",
                "n_old_only_GeneSymbol": "",
                "source": str(SIG_NEW),
                "filter": "file missing",
            }
        )

    if SIG_OLD.exists():
        so = pd.read_csv(SIG_OLD, usecols=["smPEP_ID", "GeneSymbol"], dtype=str, low_memory=False)
        so["GeneSymbol"] = so["GeneSymbol"].astype(str).str.strip()
        so = so.drop_duplicates(subset=["smPEP_ID"], keep="first")
        so_tcga = so.loc[so["GeneSymbol"].isin(tcga)]
        genes_old = set(so_tcga["GeneSymbol"].astype(str).unique())
        peps_old = so_tcga["smPEP_ID"].nunique()
        rows_compare.append(
            {
                "set_label": "old_significant_peps_tcga_genes_only",
                "n_unique_smPEP_ID": int(peps_old),
                "n_unique_GeneSymbol": len(genes_old),
                "n_intersection_GeneSymbol": "",
                "n_new_only_GeneSymbol": "",
                "n_old_only_GeneSymbol": "",
                "source": str(SIG_OLD),
                "filter": "GeneSymbol in TCGA lncRNA expression matrix columns",
            }
        )
    else:
        rows_compare.append(
            {
                "set_label": "old_significant_peps_missing",
                "n_unique_smPEP_ID": "",
                "n_unique_GeneSymbol": "",
                "n_intersection_GeneSymbol": "",
                "n_new_only_GeneSymbol": "",
                "n_old_only_GeneSymbol": "",
                "source": str(SIG_OLD),
                "filter": "file missing",
            }
        )

    inter_g: set[str] = set()
    new_only_g: set[str] = set()
    old_only_g: set[str] = set()
    if SIG_NEW.exists() and SIG_OLD.exists():
        inter_g = genes_new & genes_old
        new_only_g = genes_new - genes_old
        old_only_g = genes_old - genes_new
        rows_compare.append(
            {
                "set_label": "comparison_GeneSymbol_tcga_only",
                "n_unique_smPEP_ID": "",
                "n_unique_GeneSymbol": "",
                "n_intersection_GeneSymbol": len(inter_g),
                "n_new_only_GeneSymbol": len(new_only_g),
                "n_old_only_GeneSymbol": len(old_only_g),
                "source": "",
                "filter": "GeneSymbol sets: new_tcga_genes_only vs old_tcga_genes_only",
            }
        )

    OUT_COMPARE.mkdir(exist_ok=True)
    pd.DataFrame(rows_compare).to_csv(GENE_COMPARE_CSV, index=False)

    if inter_g or new_only_g or old_only_g:
        (OUT_COMPARE / "significant_peptide_genes_tcga_intersection.txt").write_text(
            "\n".join(sorted(inter_g)), encoding="utf-8"
        )
        (OUT_COMPARE / "significant_peptide_genes_tcga_new_only.txt").write_text(
            "\n".join(sorted(new_only_g)), encoding="utf-8"
        )
        (OUT_COMPARE / "significant_peptide_genes_tcga_old_only.txt").write_text(
            "\n".join(sorted(old_only_g)), encoding="utf-8"
        )

    print(f"Wrote {GENE_COMPARE_CSV}")
    for r in rows_compare:
        print(r)


if __name__ == "__main__":
    main()
