"""
Match SmProt peptide rows to ``significant_lncs.csv`` (canonical significant lncRNA genes).

Used by ``export_tcga_filtered_peptides_fasta.py`` after ``smprot_filtered.tsv`` is
restricted to FASTA-exportable peptides, to refresh ``significant_lnc_peptides_full.tsv``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SIG_LNCS = DATA / "significant_lncs.csv"
SIG_PEPT_FULL = DATA / "significant_lnc_peptides_full.tsv"


def write_significant_lnc_peptides_full(smprot_filtered_path: Path) -> int:
    """
    Read ``smprot_filtered.tsv`` (already FASTA-filtered) and write rows whose gene
    matches ``significant_lncs.csv`` to ``significant_lnc_peptides_full.tsv``.
    Returns number of rows written.
    """
    if not SIG_LNCS.exists():
        raise FileNotFoundError(f"Missing {SIG_LNCS}; run build_significant_lncs_smprot.py first.")
    sig = pd.read_csv(SIG_LNCS, sep=",", low_memory=False)
    filt = pd.read_csv(smprot_filtered_path, sep="\t", low_memory=False)

    sig_ids = set(sig.loc[sig["gene_id"].notna(), "gene_id"].astype(str))
    sig_ids_base = {gid.split(".", 1)[0] for gid in sig_ids}
    sig_names = set(sig["gene_name"].astype(str))

    def gene_match(gene_id: str, gene_sym: str) -> bool:
        gid = str(gene_id)
        if gid in sig_ids:
            return True
        base = gid.split(".", 1)[0]
        if base in sig_ids_base:
            return True
        if str(gene_sym) in sig_names:
            return True
        return False

    m = filt.apply(lambda r: gene_match(r["GeneID"], r["GeneSymbol"]), axis=1)
    matched = filt.loc[m].copy()
    matched.to_csv(SIG_PEPT_FULL, sep="\t", index=False)
    print(
        f"significant_lnc_peptides_full.tsv: {len(matched)} rows "
        f"({matched['smPEP_ID'].nunique()} unique smPEP_ID) -> {SIG_PEPT_FULL}"
    )
    return len(matched)
