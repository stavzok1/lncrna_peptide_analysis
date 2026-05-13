"""
Restrict smprot_filtered.tsv to peptides that can be tied to a TCGA expression column.

**Matching (two routes, OR):**

1. **Gene symbol:** SmProt ``GeneSymbol`` equals a column name in
   ``primary_exp_stage_lncRNA.csv`` or ``primary_exp_metastasis_lncRNA.csv`` (metadata
   columns stripped). Same as before.

2. **Ensembl gene ID:** SmProt ``GeneID`` (versioned ENSG… or base ENSG) is looked up in
   ``data/lncrna_genes_small.csv`` → ``gene_name``; if that **gene_name** is a TCGA matrix
   column, the row is kept. This recovers cases where SmProt’s symbol string differs from
   the matrix header but the ENSG maps to the TCGA gene name.

The TCGA gene set is the **union** of stage and metastasis matrix column names.

Also enforces >=10 aa from StartOnTrans/StopOnTrans if smprot_filtered was built
with an older pipeline (redundant when build_significant_lncs_smprot.py is current).

Writes diff reports under ``data/reports/`` when the output TSV already exists.
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



import numpy as np
import pandas as pd

DATA = ROOT / "data"

STAGE = DATA / "primary_exp_stage_lncRNA.csv"
META = DATA / "primary_exp_metastasis_lncRNA.csv"
SMPROT_FILT = DATA / "smprot_filtered.tsv"
OUT = DATA / "smprot_filtered_tcga_expr_genes.tsv"
LNCRNA_MAP = DATA / "lncrna_genes_small.csv"
REPORTS = DATA / "reports"

META_STAGE = {"sample_id", "cancer_type", "ajcc_t", "ajcc_m", "stage"}
META_META = {"sample_id", "cancer_type", "stage", "M_stage"}
MIN_AA_LENGTH = 10


def peptide_aa_length(df: pd.DataFrame) -> pd.Series:
    s = pd.to_numeric(df["StartOnTrans"], errors="coerce")
    e = pd.to_numeric(df["StopOnTrans"], errors="coerce")
    return ((e - s) // 3).astype("Int64")


def expr_genes_from_header(path: Path, meta: set[str]) -> set[str]:
    df = pd.read_csv(path, nrows=0)
    return {c for c in df.columns if c not in meta}


def build_ensg_to_gene_name(path: Path) -> dict[str, str]:
    """Map versioned ENSG, then unversioned base ENSG (first occurrence wins for base)."""
    df = pd.read_csv(path, usecols=["gene_id", "gene_name"], low_memory=False)
    m: dict[str, str] = {}
    for _, row in df.iterrows():
        gid = str(row["gene_id"]).strip()
        name = str(row["gene_name"]).strip()
        if not gid or not name:
            continue
        m[gid] = name
        base = gid.split(".", 1)[0]
        m.setdefault(base, name)
    return m


def resolve_gene_name(gene_id: str, ensg_map: dict[str, str]) -> str | None:
    g = str(gene_id).strip()
    if not g:
        return None
    if g in ensg_map:
        return ensg_map[g]
    base = g.split(".", 1)[0]
    return ensg_map.get(base)


def resolve_series(gene_ids: pd.Series, ensg_map: dict[str, str]) -> pd.Series:
    uniq = gene_ids.astype(str).str.strip().unique()
    cache: dict[str, str | None] = {}
    for u in uniq:
        cache[str(u)] = resolve_gene_name(str(u), ensg_map)
    return gene_ids.astype(str).str.strip().map(cache)


def write_reports(
    old_df: pd.DataFrame | None,
    new_df: pd.DataFrame,
    *,
    n_sym_only: int,
    n_ensg_only: int,
    n_both: int,
) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    new_ids = set(new_df["smPEP_ID"].astype(str).str.strip())
    lines: list[str] = []

    lines.append("TCGA-matrix peptide filter summary")
    lines.append("=" * 60)
    lines.append(f"Output rows (unique smPEP_ID): {len(new_df)}")
    lines.append(f"Matched via GeneSymbol in matrix columns: {n_sym_only + n_both}")
    lines.append(f"Matched only via GeneID -> lncrna_genes_small gene_name: {n_ensg_only}")
    lines.append(f"Matched by both symbol and mapped name: {n_both}")
    lines.append("")

    if old_df is not None:
        old_ids = set(old_df["smPEP_ID"].astype(str).str.strip())
        removed = old_ids - new_ids
        added = new_ids - old_ids
        lines.append(f"Compared to previous {OUT.name}: removed {len(removed)}, added {len(added)}")
        if removed:
            rem_df = old_df.loc[old_df["smPEP_ID"].astype(str).str.strip().isin(removed)].copy()
            rem_path = REPORTS / "tcga_expr_genes_removed_smpep.tsv"
            rem_df.to_csv(rem_path, sep="\t", index=False)
            lines.append(f"Removed rows (full old table columns): {rem_path}")
        if added:
            add_df = new_df.loc[new_df["smPEP_ID"].astype(str).str.strip().isin(added)].copy()
            add_path = REPORTS / "tcga_expr_genes_added_smpep.tsv"
            add_df.to_csv(add_path, sep="\t", index=False)
            lines.append(f"Added rows: {add_path}")
    else:
        lines.append(f"No prior {OUT.name} on disk — skipped removed/added file diff.")
    lines.append("")

    summary_path = REPORTS / "tcga_expr_genes_filter_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {summary_path}")

    import report_smprot_pipeline_stages as smrep

    smrep.main()


def main() -> None:
    g_stage = expr_genes_from_header(STAGE, META_STAGE)
    g_meta = expr_genes_from_header(META, META_META)
    tcga_genes = g_stage | g_meta
    if g_stage != g_meta:
        only_s = sorted(g_stage - g_meta)[:10]
        only_m = sorted(g_meta - g_stage)[:10]
        print("Note: gene column sets differ between matrices.")
        print(f"  only stage (first 10): {only_s}")
        print(f"  only meta (first 10): {only_m}")

    if not LNCRNA_MAP.exists():
        raise FileNotFoundError(f"Need gene map for Ensembl matching: {LNCRNA_MAP}")
    ensg_map = build_ensg_to_gene_name(LNCRNA_MAP)

    old_df: pd.DataFrame | None = None
    if OUT.exists():
        old_df = pd.read_csv(OUT, sep="\t", low_memory=False)

    filt = pd.read_csv(SMPROT_FILT, sep="\t", low_memory=False)
    n0 = len(filt)
    aa_len = peptide_aa_length(filt)
    filt = filt.loc[aa_len.notna() & (aa_len >= MIN_AA_LENGTH)].copy()
    n_after_len = len(filt)

    sym = filt["GeneSymbol"].astype(str).str.strip()
    gid = filt["GeneID"].astype(str).str.strip()
    resolved = resolve_series(gid, ensg_map)
    via_sym = sym.isin(tcga_genes)
    via_ensg = resolved.notna() & resolved.astype(str).isin(tcga_genes)
    in_tcga = via_sym | via_ensg

    n_sym_only = int((via_sym & ~via_ensg).sum())
    n_ensg_only = int((~via_sym & via_ensg).sum())
    n_both = int((via_sym & via_ensg).sum())

    sub = filt.loc[in_tcga].copy()
    mloc = sub.index
    vs = via_sym.loc[mloc].to_numpy()
    ve = via_ensg.loc[mloc].to_numpy()
    sub["TCGA_match_via"] = np.where(
        vs & ve,
        "gene_symbol+ensembl_map",
        np.where(vs, "gene_symbol", "ensembl_gene_id"),
    )

    for col in ("TISPvalue", "RiboPvalue"):
        sub[col] = pd.to_numeric(sub[col], errors="coerce")
    sub = sub.sort_values(
        by=["TISPvalue", "RiboPvalue", "RiboID"],
        ascending=[True, True, True],
        na_position="last",
    )
    sub = sub.drop_duplicates(subset=["smPEP_ID"], keep="first")
    sub.to_csv(OUT, sep="\t", index=False)

    n_sym_in_tcga = int(via_sym.sum())
    n_any = int(in_tcga.sum())
    print(f"TCGA expression gene symbols (union): {len(tcga_genes)}")
    print(f"Ensembl gene_id -> gene_name map entries: {len(ensg_map)}")
    print(
        f"smprot_filtered.tsv rows: {n0}, after >= {MIN_AA_LENGTH} aa: {n_after_len}, "
        f"GeneSymbol in TCGA cols: {n_sym_in_tcga}, symbol OR Ensembl-mapped name in TCGA: {n_any}"
    )
    print(
        f"  match breakdown (pre-dedupe rows): symbol-only {n_sym_only}, "
        f"ensembl_map-only {n_ensg_only}, both {n_both}"
    )
    print(f"After dedupe smPEP_ID: {len(sub)} rows -> {OUT}")

    write_reports(old_df, sub, n_sym_only=n_sym_only, n_ensg_only=n_ensg_only, n_both=n_both)


if __name__ == "__main__":
    main()
