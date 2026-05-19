"""
Match SmProt peptide rows to ``significant_lncs.csv`` (canonical significant lncRNA genes).

Used by ``export_tcga_filtered_peptides_fasta.py`` after ``smprot_filtered.tsv`` is
restricted to FASTA-exportable peptides, to refresh ``significant_lnc_peptides_full.tsv``.
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
SIG_LNCS = DATA / "significant_lncs.csv"
SIG_PEPT_FULL = DATA / "significant_lnc_peptides_full.tsv"
SIG_PEPT = DATA / "significant_lnc_peptides.tsv"
CANONICAL_GENES = DATA / "canonical_significant_lncRNAs.txt"
LNCRNA_MAP = DATA / "lncrna_genes_small.csv"
LIMMA_Z = ROOT / "tr_lncrna_output" / "limma" / "limma_z_intersection_genes.txt"


def _ens_base(gene_id: str) -> str:
    gid = str(gene_id).strip()
    if not gid.startswith("ENSG"):
        return gid
    return gid.split(".", 1)[0]


def gene_matches_significant_lnc(
    gene_id: str,
    gene_sym: str,
    *,
    sig_ids: set[str],
    sig_ids_base: set[str],
    sig_names: set[str],
) -> bool:
    gid = str(gene_id).strip()
    if gid in sig_ids:
        return True
    if _ens_base(gid) in sig_ids_base:
        return True
    if str(gene_sym).strip() in sig_names:
        return True
    return False


def significant_lnc_match_sets() -> tuple[set[str], set[str], set[str]]:
    """(versioned ENSG ids, unversioned ENSG bases, canonical gene symbols) for Tr matching."""
    if SIG_LNCS.exists():
        sig = pd.read_csv(SIG_LNCS, low_memory=False)
        sig_ids = set(sig.loc[sig["gene_id"].notna(), "gene_id"].astype(str).str.strip())
        sig_ids_base = {_ens_base(g) for g in sig_ids}
        sig_names = set(sig["gene_name"].astype(str).str.strip()) - {""}
        return sig_ids, sig_ids_base, sig_names

    path = CANONICAL_GENES if CANONICAL_GENES.exists() else LIMMA_Z
    if not path.exists():
        raise FileNotFoundError(
            f"Need {SIG_LNCS} or {CANONICAL_GENES} (or {LIMMA_Z}) for Ensembl-aware Tr matching."
        )
    sig_names = {ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()}
    sig_ids: set[str] = set()
    sig_ids_base: set[str] = set()
    if LNCRNA_MAP.exists():
        gmap = pd.read_csv(LNCRNA_MAP, usecols=["gene_id", "gene_name"], dtype=str)
        gmap["gene_name"] = gmap["gene_name"].str.strip()
        gmap["gene_id"] = gmap["gene_id"].str.strip()
        sub = gmap.loc[gmap["gene_name"].isin(sig_names)]
        sig_ids = set(sub["gene_id"])
        sig_ids_base = {_ens_base(g) for g in sig_ids}

    # Peptide export tables carry SmProt ENSG ids for Tr genes (symbol aliases included).
    for pep_path in (SIG_PEPT_FULL, SIG_PEPT):
        if not pep_path.exists():
            continue
        pep = pd.read_csv(pep_path, sep="\t", usecols=["GeneID", "GeneSymbol"], low_memory=True)
        for gid in pep["GeneID"].astype(str).str.strip():
            if gid and gid.lower() != "nan":
                sig_ids.add(gid)
                sig_ids_base.add(_ens_base(gid))
        for sym in pep["GeneSymbol"].astype(str).str.strip():
            if sym and sym.lower() != "nan":
                sig_names.add(sym)

    return sig_ids, sig_ids_base, sig_names


def peptide_rows_match_significant_lnc(df: pd.DataFrame) -> pd.Series:
    """True per row: SmProt GeneID and/or GeneSymbol matches canonical Tr lncRNAs."""
    need = {"GeneID", "GeneSymbol"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"peptide table missing columns {sorted(missing)}")
    sig_ids, sig_ids_base, sig_names = significant_lnc_match_sets()
    return df.apply(
        lambda r: gene_matches_significant_lnc(
            r["GeneID"],
            r["GeneSymbol"],
            sig_ids=sig_ids,
            sig_ids_base=sig_ids_base,
            sig_names=sig_names,
        ),
        axis=1,
    )


def write_significant_lnc_peptides_full(smprot_filtered_path: Path) -> int:
    """
    Read ``smprot_filtered.tsv`` (already FASTA-filtered) and write rows whose gene
    matches ``significant_lncs.csv`` to ``significant_lnc_peptides_full.tsv``.
    Returns number of rows written.
    """
    if not SIG_LNCS.exists():
        raise FileNotFoundError(f"Missing {SIG_LNCS}; run build_significant_lncs_smprot.py first.")
    filt = pd.read_csv(smprot_filtered_path, sep="\t", low_memory=False)
    m = peptide_rows_match_significant_lnc(filt)
    matched = filt.loc[m].copy()
    matched.to_csv(SIG_PEPT_FULL, sep="\t", index=False)
    print(
        f"significant_lnc_peptides_full.tsv: {len(matched)} rows "
        f"({matched['smPEP_ID'].nunique()} unique smPEP_ID) -> {SIG_PEPT_FULL}"
    )
    return len(matched)
