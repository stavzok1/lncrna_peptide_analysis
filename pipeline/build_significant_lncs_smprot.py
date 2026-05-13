"""
Build significant_lncs.csv from the canonical significant lncRNA gene list,
then SmProt2 curation: lncRNA rows, p-value filter, >=10 aa (from transcript coords), dedupe, peptide match.
``data/significant_lnc_peptides_full.tsv`` is **not** written here. After this script, run
``export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa``:
that export **narrows** ``smprot_filtered.tsv`` to FASTA-exportable peptides only, then writes
``significant_lnc_peptides_full.tsv`` and (by default) re-exports ``significant_lnc_peptides.faa`` and runs
``sync_significant_lnc_peptides_with_fasta.py`` for ``significant_lnc_peptides.tsv``.

Canonical genes = limma (FDR<0.05) ∩ z-score (|z|>=3) intersection, produced by
tr_limma_de.R (written to data/canonical_significant_lncRNAs.txt when z union exists).
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
LIMMA_Z_INTERSECTION = ROOT / "tr_lncrna_output" / "limma" / "limma_z_intersection_genes.txt"
CANONICAL_GENES = DATA / "canonical_significant_lncRNAs.txt"
MAP_PATH = DATA / "lncrna_genes_small.csv"
SMPROT = DATA / "SmProt2.txt"

LNCRNA_PEPT = DATA / "lncRNA_peptides.tsv"
SMPROT_FILT = DATA / "smprot_filtered.tsv"
SIG_LNCS = DATA / "significant_lncs.csv"

CHUNKSIZE = 500_000
# SmProt coordinates are on the transcript in nucleotides (0-based half-open); ORF length in aa.
MIN_AA_LENGTH = 10


def peptide_aa_length(df: pd.DataFrame) -> pd.Series:
    s = pd.to_numeric(df["StartOnTrans"], errors="coerce")
    e = pd.to_numeric(df["StopOnTrans"], errors="coerce")
    return ((e - s) // 3).astype("Int64")


def load_canonical_gene_names() -> pd.DataFrame:
    """Return one-column DataFrame gene_name; sync data/canonical from limma output."""
    if LIMMA_Z_INTERSECTION.exists():
        text = LIMMA_Z_INTERSECTION.read_text(encoding="utf-8", errors="replace")
        CANONICAL_GENES.write_text(text, encoding="utf-8")
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    elif CANONICAL_GENES.exists():
        lines = [ln.strip() for ln in CANONICAL_GENES.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
    else:
        raise FileNotFoundError(
            "Missing canonical gene list. Run Rscript tr_limma_de.R <project_root> first "
            f"(expects {LIMMA_Z_INTERSECTION} or {CANONICAL_GENES})."
        )
    return pd.DataFrame({"gene_name": lines})


def main() -> None:
    names = load_canonical_gene_names()
    gmap = pd.read_csv(MAP_PATH, usecols=["gene_id", "gene_name"])
    # One gene_id per gene_name in reference; if duplicates, keep first
    gmap = gmap.drop_duplicates(subset=["gene_name"], keep="first")
    sig = names.merge(gmap, on="gene_name", how="left")
    sig.to_csv(SIG_LNCS, index=False)
    n_mapped = sig["gene_id"].notna().sum()
    print(f"significant_lncs.csv: {len(sig)} genes, {n_mapped} with gene_id, {len(sig) - n_mapped} unmatched")

    # Pass 1: lncRNA rows only
    first = True
    total_lnc = 0
    for chunk in pd.read_csv(SMPROT, sep="\t", chunksize=CHUNKSIZE, low_memory=False):
        sub = chunk.loc[chunk["GeneType"].astype(str) == "lncRNA"].copy()
        total_lnc += len(sub)
        sub.to_csv(LNCRNA_PEPT, sep="\t", mode="w" if first else "a", header=first, index=False)
        first = False
    print(f"lncRNA_peptides.tsv: {total_lnc} rows written")

    lnc = pd.read_csv(LNCRNA_PEPT, sep="\t", low_memory=False)
    for col in ("TISPvalue", "RiboPvalue"):
        lnc[col] = pd.to_numeric(lnc[col], errors="coerce")
    filt = lnc[(lnc["TISPvalue"] < 0.05) & (lnc["RiboPvalue"] < 0.05)].copy()
    n_before_len = len(filt)
    aa_len = peptide_aa_length(filt)
    filt = filt.loc[aa_len.notna() & (aa_len >= MIN_AA_LENGTH)].copy()
    n_dropped_len = n_before_len - len(filt)
    # Dedupe smPEP_ID: keep the row with strongest TIS then Ribo p-values (not arbitrary file order).
    filt = filt.sort_values(
        by=["TISPvalue", "RiboPvalue", "RiboID"],
        ascending=[True, True, True],
        na_position="last",
    )
    filt = filt.drop_duplicates(subset=["smPEP_ID"], keep="first")
    filt.to_csv(SMPROT_FILT, sep="\t", index=False)
    print(
        f"smprot_filtered.tsv (pre–FASTA filter): {len(filt)} unique smPEP_ID after p<0.05, "
        f">= {MIN_AA_LENGTH} aa (dropped {n_dropped_len} rows on length), dedupe"
    )
    print(
        "Next: python export_tcga_filtered_peptides_fasta.py "
        "--peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa\n"
        "  → restricts smprot_filtered.tsv to FASTA-exportable rows, writes significant_lnc_peptides_full.tsv, "
        "re-exports significant_lnc_peptides.faa, and runs sync_significant_lnc_peptides_with_fasta.py "
        "(use --no-chain-significant-fasta on that export to skip the significant FASTA + sync chain)."
    )


if __name__ == "__main__":
    main()
