"""
Write ``data/significant_lnc_peptides.tsv`` as the subset of rows from
``data/significant_lnc_peptides_full.tsv`` whose ``smPEP_ID`` appears in
``data/significant_lnc_peptides.faa`` (exportable / analyzed peptides).

Run after ``export_tcga_filtered_peptides_fasta.py --peptides-tsv data/significant_lnc_peptides_full.tsv``
(creates ``significant_lnc_peptides.faa``). When you export from ``data/smprot_filtered.tsv``,
that script normally runs this sync for you after refreshing ``*_full.tsv``.

Downstream figures and NetMHC prep should use ``significant_lnc_peptides.tsv`` so
coordinates without a successful in-silico translation are excluded.
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


import shutil

import pandas as pd
from Bio import SeqIO

DATA = ROOT / "data"
FULL_TSV = DATA / "significant_lnc_peptides_full.tsv"
OUT_TSV = DATA / "significant_lnc_peptides.tsv"
SIG_FAA = DATA / "significant_lnc_peptides.faa"


def fasta_smpep_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    for rec in SeqIO.parse(path, "fasta"):
        part = rec.description.split("|", 1)[0].strip()
        if part.startswith(">"):
            part = part[1:].strip()
        ids.add(part)
    return ids


def main() -> None:
    if not SIG_FAA.exists():
        print(f"Missing {SIG_FAA}; export peptides first.", file=sys.stderr)
        raise SystemExit(1)

    if not FULL_TSV.exists():
        if OUT_TSV.exists():
            shutil.copy2(OUT_TSV, FULL_TSV)
            print(f"Created {FULL_TSV.name} from existing {OUT_TSV.name} (SmProt full list).")
        else:
            print(f"Missing {FULL_TSV} and {OUT_TSV}. Run build_significant_lncs_smprot.py first.", file=sys.stderr)
            raise SystemExit(1)

    have = fasta_smpep_ids(SIG_FAA)
    df = pd.read_csv(FULL_TSV, sep="\t", dtype=str, low_memory=False)
    if "smPEP_ID" not in df.columns:
        raise SystemExit(f"{FULL_TSV} has no smPEP_ID column")

    df["smPEP_ID"] = df["smPEP_ID"].astype(str).str.strip()
    ordered_ids: list[str] = []
    seen_row: set[str] = set()
    for sid in df["smPEP_ID"]:
        if sid in have and sid not in seen_row:
            ordered_ids.append(sid)
            seen_row.add(sid)

    sub = df.loc[df["smPEP_ID"].isin(have)].copy()
    sub["_k"] = pd.Categorical(sub["smPEP_ID"], categories=ordered_ids, ordered=True)
    sub = sub.sort_values("_k").drop(columns=["_k"])

    n_full = df["smPEP_ID"].nunique()
    n_out = sub["smPEP_ID"].nunique()
    missing = sorted(set(df["smPEP_ID"].unique()) - have)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"FASTA smPEP_ID: {len(have)}")
    print(f"Full TSV unique smPEP_ID: {n_full}; wrote analyzed TSV rows: {len(sub)} ({n_out} unique)")
    if missing:
        print(f"Excluded {len(missing)} smPEP_ID not present in FASTA (e.g. {missing[:3]}…).")
    print(f"Saved: {OUT_TSV}")


if __name__ == "__main__":
    main()
