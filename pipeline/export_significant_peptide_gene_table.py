"""
Build one row per unique SmProt GeneSymbol in **analyzed** ``significant_lnc_peptides.tsv``
(exportable peptides; typically ~124–127 genes depending on FASTA sync),
with SmProt IDs, canonical (limma∩z) gene names/ENSG from significant_lncs.csv via
versionless ENSG match, and TCGA expression-matrix symbol / ENSG from lncrna_genes_small.csv.

Output: data/significant_lnc_peptide_gene_map_127.csv (filename kept for downstream paths).
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
SIG_PEPT = DATA / "significant_lnc_peptides.tsv"
SIG_LNCS = DATA / "significant_lncs.csv"
LNCMAP = DATA / "lncrna_genes_small.csv"
STAGE_EXPR = DATA / "primary_exp_stage_lncRNA.csv"
META_EXPR = DATA / "primary_exp_metastasis_lncRNA.csv"
OUT_CSV = DATA / "significant_lnc_peptide_gene_map_127.csv"

META_STAGE = {"sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage"}
META_META = {"sample_id", "cancer_type", "stage", "M_stage"}


def ens_base(gid: str) -> str:
    gid = str(gid).strip()
    if not gid.startswith("ENSG"):
        return gid
    return gid.split(".", 1)[0]


def tcga_expr_symbols() -> set[str]:
    c1 = pd.read_csv(STAGE_EXPR, nrows=0).columns
    c2 = pd.read_csv(META_EXPR, nrows=0).columns
    g1 = {str(c) for c in c1 if c not in META_STAGE}
    g2 = {str(c) for c in c2 if c not in META_META}
    return g1 | g2


def main() -> None:
    pep = pd.read_csv(SIG_PEPT, sep="\t", dtype=str, low_memory=False)
    for c in ("GeneSymbol", "GeneID", "TranscriptID", "smPEP_ID"):
        pep[c] = pep[c].astype(str).str.strip()

    sig = pd.read_csv(SIG_LNCS, dtype=str)
    sig["gene_name"] = sig["gene_name"].astype(str).str.strip()
    sig["gene_id"] = sig["gene_id"].astype(str).str.strip()
    base_to_canon: dict[str, tuple[str, str]] = {}
    for _, r in sig.iterrows():
        b = ens_base(r["gene_id"])
        if b not in base_to_canon:
            base_to_canon[b] = (r["gene_name"], r["gene_id"])

    lnc = pd.read_csv(LNCMAP, usecols=["gene_id", "gene_name"], dtype=str)
    lnc["gene_name"] = lnc["gene_name"].astype(str).str.strip()
    lnc["gene_id"] = lnc["gene_id"].astype(str).str.strip()
    sym_to_gid = lnc.drop_duplicates(subset=["gene_name"], keep="first").set_index("gene_name")[
        "gene_id"
    ].to_dict()

    tcga = tcga_expr_symbols()

    rows: list[dict] = []
    for sym, g in pep.groupby("GeneSymbol", sort=True):
        gids = sorted({x for x in g["GeneID"].unique() if x and x != "nan"})
        tids = sorted({x for x in g["TranscriptID"].unique() if x and x != "nan"})
        pids = sorted({x for x in g["smPEP_ID"].unique() if x and x != "nan"})

        canon_name = ""
        canon_id = ""
        for gid in gids:
            b = ens_base(gid)
            if b in base_to_canon:
                canon_name, canon_id = base_to_canon[b]
                break

        tcga_sym = ""
        tcga_gid = ""
        note = ""
        if canon_name and canon_name in tcga:
            tcga_sym = canon_name
            tcga_gid = sym_to_gid.get(canon_name, canon_id)
            note = "tcga_symbol_via_canonical_ensg_match"
        elif sym in tcga:
            tcga_sym = sym
            tcga_gid = sym_to_gid.get(sym, "")
            note = "tcga_symbol_smprot_only"
        else:
            note = "no_tcga_matrix_column"

        rows.append(
            {
                "smprot_GeneSymbol": sym,
                "smprot_GeneID_primary": gids[0] if gids else "",
                "smprot_GeneID_all": ";".join(gids),
                "smprot_TranscriptID_all": ";".join(tids),
                "n_significant_lnc_peptides": len(pids),
                "smPEP_IDs": ";".join(pids),
                "canonical_gene_name": canon_name,
                "canonical_gene_id": canon_id,
                "tcga_matrix_gene_symbol": tcga_sym,
                "tcga_matrix_gene_id": tcga_gid,
                "in_tcga_expression_panel": "Y" if tcga_sym else "N",
                "mapping_note": note,
            }
        )

    out = pd.DataFrame(rows).sort_values("smprot_GeneSymbol")
    if len(out) != 127:
        print(
            f"Note: historical run had 127 unique SmProt symbols in the full significant table; "
            f"this export has {len(out)} (reflects analyzed `significant_lnc_peptides.tsv` after FASTA sync)."
        )
    out.to_csv(OUT_CSV, index=False)
    n_panel = (out["in_tcga_expression_panel"] == "Y").sum()
    print(f"Wrote {OUT_CSV} ({len(out)} rows, {n_panel} with TCGA matrix symbol/id).")


if __name__ == "__main__":
    main()
