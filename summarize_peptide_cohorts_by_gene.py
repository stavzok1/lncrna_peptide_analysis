"""
Build gene-level summaries for SmProt-filtered, TCGA-expression-filtered, and
Tr-lncRNA (significant ribo-seq) peptide cohorts.

AA peptide lengths are read from ``data/smprot_all_filtered_peptides.faa`` (all
cohort ``smPEP_ID`` values used here are present in that file).

Optional: longest transcript length per gene (nucleotides) from a GENCODE-style
transcript FASTA (``>ENST...|ENSG...|...``). Auto-detects ``data/gencode*transcript*.fa*``
if present; otherwise pass ``--transcripts-fa`` or columns are left empty.

Writes:
  data/reports/peptide_cohort_gene_summary.tsv
  data/reports/tcga_filtered_peptides_per_gene.tsv
  data/reports/tr_lncrna_filtered_peptides_per_gene.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
REPORTS = DATA / "reports"
SMPEP_FAA = DATA / "smprot_all_filtered_peptides.faa"

DEFAULTS = {
    "smprot_filtered": (
        DATA / "smprot_filtered.tsv",
        "2769 SmProt-filtered peptides (all Ribo-seq evidence)",
    ),
    "tcga_filtered": (
        DATA / "smprot_filtered_tcga_expr_genes.tsv",
        "2606 TCGA-expression-filtered peptides",
    ),
    "tr_lncrna": (
        DATA / "significant_lnc_peptides.tsv",
        "501 Tr-lncRNA filtered peptides (significant ribo-seq, analyzed set)",
    ),
}


def _load_smpep_aa_lengths(faa_path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    for rec in SeqIO.parse(str(faa_path), "fasta"):
        sid = str(rec.id).split("|", 1)[0].strip()
        seq = str(rec.seq).replace(" ", "").upper()
        aa = "".join(c for c in seq if c in "ACDEFGHIKLMNPQRSTVWY")
        out[sid] = len(aa)
    return out


def _lookup_gene_tx(gene_tx: dict[str, int], gene_id: str) -> int | None:
    if gene_id in gene_tx:
        return int(gene_tx[gene_id])
    if isinstance(gene_id, str) and "." in gene_id:
        base = gene_id.rsplit(".", 1)[0]
        if base in gene_tx:
            return int(gene_tx[base])
    return None


def _load_gene_longest_transcript_nt(fa_path: Path) -> dict[str, int]:
    """
    GENCODE transcript FASTA: first '|' field ENST, second ENSG.
    Keep max len(seq) per ENSG (any isoform).
    """
    gene_max: dict[str, int] = {}
    for rec in SeqIO.parse(str(fa_path), "fasta"):
        parts = rec.description.split("|")
        if len(parts) < 2:
            continue
        ensg = parts[1].strip()
        if not ensg.startswith("ENSG"):
            continue
        L = len(rec.seq)
        gene_max[ensg] = max(gene_max.get(ensg, 0), L)
        if "." in ensg:
            base = ensg.rsplit(".", 1)[0]
            gene_max[base] = max(gene_max.get(base, 0), L)
    return gene_max


def _discover_transcript_fasta(data_dir: Path) -> Path | None:
    for p in sorted(data_dir.glob("gencode*transcript*.fa*")) + sorted(
        data_dir.glob("gencode*transcripts*.fa*")
    ):
        if p.is_file():
            return p
    return None


def _attach_aa_len(df: pd.DataFrame, lens: dict[str, int]) -> pd.DataFrame:
    out = df.copy()
    out["smPEP_ID"] = out["smPEP_ID"].astype(str).str.strip()
    out["aa_len"] = out["smPEP_ID"].map(lens)
    miss = int(out["aa_len"].isna().sum())
    if miss:
        print(f"Warning: {miss} smPEP_ID without AA length in FASTA", file=sys.stderr)
    return out


def peptides_per_gene(df: pd.DataFrame, gene_tx: dict[str, int] | None) -> pd.DataFrame:
    g = (
        df.groupby("GeneID", as_index=False)
        .agg(
            GeneSymbol=("GeneSymbol", "first"),
            n_peptides=("smPEP_ID", "count"),
            max_peptide_aa=("aa_len", "max"),
            n_peptides_aa_ge_30=("aa_len", lambda s: int((s >= 30).sum())),
        )
        .sort_values("n_peptides", ascending=False, kind="mergesort")
        .reset_index(drop=True)
    )
    if gene_tx:
        g["max_transcript_nt"] = g["GeneID"].map(lambda gid: _lookup_gene_tx(gene_tx, str(gid)))
        g["max_transcript_nt"] = g["max_transcript_nt"].astype("Int64")
    else:
        g["max_transcript_nt"] = pd.NA
    return g


def _cohort_tx_stat(df: pd.DataFrame, gene_tx: dict[str, int] | None) -> tuple[int | None, int | None]:
    """Max and median of (longest transcript per gene) over genes present in cohort."""
    if not gene_tx:
        return None, None
    vals: list[int] = []
    for gid in df["GeneID"].astype(str).unique():
        v = _lookup_gene_tx(gene_tx, gid)
        if v is not None:
            vals.append(v)
    if not vals:
        return None, None
    vs = sorted(vals)
    return vs[-1], vs[len(vs) // 2]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=REPORTS, help="Directory for output TSVs")
    ap.add_argument(
        "--peptide-fasta",
        type=Path,
        default=SMPEP_FAA,
        help="FASTA of smPROT peptides (header id = smPEP_ID) for AA lengths",
    )
    ap.add_argument(
        "--transcripts-fa",
        type=Path,
        default=None,
        help="GENCODE transcript FASTA (.fa or .gz) for longest transcript per gene; "
        "default: auto-detect under data/ or leave transcript columns empty",
    )
    ap.add_argument(
        "--no-transcript-fasta",
        action="store_true",
        help="Do not auto-detect transcript FASTA (max_transcript_nt columns empty).",
    )
    args = ap.parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.peptide_fasta.is_file():
        raise SystemExit(f"Missing peptide FASTA: {args.peptide_fasta}")

    lens = _load_smpep_aa_lengths(args.peptide_fasta)
    print(f"Loaded {len(lens)} peptide sequences from {args.peptide_fasta.relative_to(ROOT)}")

    tx_path = args.transcripts_fa
    if args.no_transcript_fasta:
        tx_path = None
    elif tx_path is None:
        tx_path = _discover_transcript_fasta(DATA)

    gene_tx: dict[str, int] | None = None
    if tx_path is not None and tx_path.is_file():
        gene_tx = _load_gene_longest_transcript_nt(tx_path)
        n_gene_keys = len({k for k in gene_tx if k.startswith("ENSG") and "." in k})
        print(
            f"Loaded longest-transcript nt length per gene ({n_gene_keys} versioned ENSG keys) "
            f"from {tx_path.relative_to(ROOT)}"
        )
    elif tx_path is not None:
        print(f"Warning: transcript FASTA not found ({tx_path}); max_transcript_nt empty.", file=sys.stderr)
    else:
        print(
            "No transcript FASTA supplied or found under data/ (gencode*transcript*.fa*); "
            "max_transcript_nt columns will be empty. Use --transcripts-fa /path/to/gencode.*.transcripts.fa.gz",
            file=sys.stderr,
        )

    rows = []
    tcga_df: pd.DataFrame | None = None
    tr_df: pd.DataFrame | None = None

    for key, (path, label) in DEFAULTS.items():
        df = pd.read_csv(path, sep="\t")
        df = _attach_aa_len(df, lens)
        n_pep = len(df)
        n_gene = int(df["GeneID"].nunique())
        max_aa_s = df["aa_len"].max()
        max_aa = int(max_aa_s) if pd.notna(max_aa_s) else ""
        n_ge30 = int((df["aa_len"] >= 30).sum()) if df["aa_len"].notna().any() else 0
        tx_max, tx_med = _cohort_tx_stat(df, gene_tx)

        row = {
            "cohort_key": key,
            "cohort_description": label,
            "source_tsv": path.relative_to(ROOT).as_posix(),
            "n_peptides": n_pep,
            "n_unique_genes_GeneID": n_gene,
            "n_unique_GeneSymbol": int(df["GeneSymbol"].nunique()),
            "max_peptide_aa_in_cohort": max_aa,
            "n_peptides_aa_ge_30_in_cohort": n_ge30,
            "max_longest_transcript_nt_among_genes": tx_max if tx_max is not None else "",
            "median_longest_transcript_nt_per_gene": tx_med if tx_med is not None else "",
        }
        rows.append(row)
        if key == "tcga_filtered":
            tcga_df = df
        if key == "tr_lncrna":
            tr_df = df

    summary = pd.DataFrame(rows)
    summary_path = out_dir / "peptide_cohort_gene_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    assert tcga_df is not None and tr_df is not None
    tcga_path = out_dir / "tcga_filtered_peptides_per_gene.tsv"
    peptides_per_gene(tcga_df, gene_tx).to_csv(tcga_path, sep="\t", index=False)

    tr_path = out_dir / "tr_lncrna_filtered_peptides_per_gene.tsv"
    peptides_per_gene(tr_df, gene_tx).to_csv(tr_path, sep="\t", index=False)

    print(summary.to_string(index=False))
    print()
    print(f"Wrote {summary_path.relative_to(ROOT)}")
    print(f"Wrote {tcga_path.relative_to(ROOT)} ({len(tcga_df)} peptides -> {tcga_df['GeneID'].nunique()} genes)")
    print(f"Wrote {tr_path.relative_to(ROOT)} ({len(tr_df)} peptides -> {tr_df['GeneID'].nunique()} genes)")


if __name__ == "__main__":
    main()
