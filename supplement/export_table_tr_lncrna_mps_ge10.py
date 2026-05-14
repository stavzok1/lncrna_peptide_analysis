"""
Emit ``table_tr_lncRNA_mps_ge10.csv``: canonical Tr lncRNAs with ≥10 MPs in the
TCGA-expression–filtered table.

``Size_aa`` = **maximum** amino-acid length among MPs of that gene (from
``data/smprot_tcga_filtered_peptides.faa``, keyed by ``smPEP_ID`` in the FASTA
header before the first ``|``).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
DEFAULT_CANON = DATA / "canonical_significant_lncRNAs.txt"
DEFAULT_TCGA = DATA / "smprot_filtered_tcga_expr_genes.tsv"
DEFAULT_FAA = DATA / "smprot_tcga_filtered_peptides.faa"
DEFAULT_OUT = Path(__file__).resolve().parent / "table_tr_lncRNA_mps_ge10.csv"


def load_peptide_lengths(faa: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    for rec in SeqIO.parse(faa, "fasta"):
        smpep = str(rec.id).split("|", 1)[0].strip()
        out[smpep] = len(rec.seq)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--canonical", type=Path, default=DEFAULT_CANON)
    ap.add_argument("--tcga-tsv", type=Path, default=DEFAULT_TCGA)
    ap.add_argument("--peptide-faa", type=Path, default=DEFAULT_FAA)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--min-mps", type=int, default=10)
    args = ap.parse_args()

    canon = {ln.strip() for ln in args.canonical.read_text(encoding="utf-8").splitlines() if ln.strip()}
    lens = load_peptide_lengths(args.peptide_faa)
    tcga = pd.read_csv(args.tcga_tsv, sep="\t")
    tcga = tcga[tcga["GeneSymbol"].astype(str).isin(canon)].copy()
    tcga["smpep_k"] = tcga["smPEP_ID"].astype(str)
    missing = tcga.loc[~tcga["smpep_k"].isin(lens), "smpep_k"]
    if len(missing):
        raise SystemExit(f"{len(missing)} smPEP_ID rows lack FASTA length (e.g. {missing.iloc[0]!r})")
    tcga["aa_len"] = tcga["smpep_k"].map(lens)

    g = (
        tcga.groupby(["GeneSymbol", "GeneID"], as_index=False)
        .agg(nMPs=("smPEP_ID", "count"), Size_aa=("aa_len", "max"), nMPs_gt30_aa=("aa_len", lambda s: int((s > 30).sum())))
    )
    g = g[g["nMPs"] >= args.min_mps].sort_values("nMPs", ascending=False)
    out = g.rename(columns={"GeneSymbol": "Gene_symbol", "GeneID": "Ensembl_ID"})
    out = out[["Gene_symbol", "Ensembl_ID", "nMPs", "Size_aa", "nMPs_gt30_aa"]]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {args.out} ({len(out)} genes)")


if __name__ == "__main__":
    main()
