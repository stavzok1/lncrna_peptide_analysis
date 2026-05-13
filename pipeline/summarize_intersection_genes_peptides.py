"""
Per-gene table: limma (FDR) ∩ z-score significant lncRNAs × TCGA-filtered SmProt peptides.

Gene list = canonical_significant_lncRNAs.txt (same as limma_z_intersection_genes.txt
when the limma pipeline has been run). No overlap with the old gene list.

Columns: gene_name, n_unique_filtered_peptides, smPEP_IDs (semicolon-separated).
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
LIMMA_OUT = ROOT / "tr_lncrna_output" / "limma"
CANONICAL = DATA / "canonical_significant_lncRNAs.txt"
LIMMA_Z_TXT = LIMMA_OUT / "limma_z_intersection_genes.txt"
PEPT_TCGA = DATA / "smprot_filtered_tcga_expr_genes.tsv"
OUT = LIMMA_OUT / "limma_z_genes_tcga_filtered_peptide_summary.csv"
OUT_MIN1_SORTED = LIMMA_OUT / "limma_z_genes_tcga_filtered_peptide_summary_min1_by_npeps.csv"


def norm_pep_id(s) -> str:
    if pd.isna(s):
        return ""
    if isinstance(s, float) and s == int(s):
        return str(int(s))
    return str(s).strip()


def load_limma_z_gene_list() -> list[str]:
    if CANONICAL.exists():
        path = CANONICAL
    elif LIMMA_Z_TXT.exists():
        path = LIMMA_Z_TXT
    else:
        raise FileNotFoundError(
            f"Missing limma∩z gene list. Expected {CANONICAL} or {LIMMA_Z_TXT}. "
            "Run: Rscript tr_limma_de.R <project_root> (and/or python build_significant_lncs_smprot.py)."
        )
    return [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def main() -> None:
    if not PEPT_TCGA.exists():
        raise FileNotFoundError(
            f"Missing {PEPT_TCGA}. Run: python filter_peptides_tcga_expr_genes.py"
        )

    genes = load_limma_z_gene_list()
    gene_set = set(genes)

    pep = pd.read_csv(PEPT_TCGA, sep="\t", usecols=["GeneSymbol", "smPEP_ID"], low_memory=False)
    pep["GeneSymbol"] = pep["GeneSymbol"].astype(str)

    rows: list[dict] = []
    for g in sorted(gene_set):
        sub = pep.loc[pep["GeneSymbol"] == g, "smPEP_ID"]
        ids = sorted({norm_pep_id(x) for x in sub if norm_pep_id(x)})
        rows.append(
            {
                "gene_name": g,
                "n_unique_filtered_peptides": len(ids),
                "smPEP_IDs": ";".join(ids),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    n_with = (out["n_unique_filtered_peptides"] > 0).sum()
    min1 = out.loc[out["n_unique_filtered_peptides"] >= 1].copy()
    min1 = min1.sort_values("n_unique_filtered_peptides", ascending=False)
    min1.to_csv(OUT_MIN1_SORTED, index=False)

    print(f"Genes in intersection: {len(out)}")
    print(f"Genes with >=1 TCGA-filtered peptide: {n_with}")
    print(f"Wrote: {OUT}")
    print(f"Wrote: {OUT_MIN1_SORTED} ({len(min1)} rows, sorted by n_unique_filtered_peptides desc)")


if __name__ == "__main__":
    main()
