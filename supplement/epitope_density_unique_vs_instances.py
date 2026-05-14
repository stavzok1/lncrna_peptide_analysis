"""
SB epitope **density** for the **same short-parent sampling design as Figure 5**
(``plot_fig5abc_netmhc_sb_triple.py``):

- **Sig lnc:** merged ``netmhcpan_sig_lnc_with_iedb.tsv`` + parent lengths from
  ``sig_parent_micropeptides.csv`` (``length_aa`` = screened ORF / micropeptide parent).

- **Coding proportional whole:** merged ``netmhcpan_coding_proportional_whole_with_iedb.tsv``
  + parent lengths from ``coding_proportional_whole_parent_micropeptides.csv``. These are
  **short proteins sampled by the proportional-whole method** used for Fig 5 — **not**
  random genomic fragments and **not** full-length UniProt chains unless you opt in below.

**SB gate:** canonical Fig 5 stack in ``netmhc_sb_core`` (``--sb-mode``; default ``full``).

**Per Fig 5 parent row (501 sig MPs vs 501 coding controls), two numerators (allele rows collapsed):**

1. **Unique 9-mers:** distinct ``Peptide`` (case-insensitive) among all SB rows pooled for
   that ``smPEP_ID`` (sig) or ``control_id`` (coding).

2. **Positional instances:** distinct ``(input_seq_id_fasta, Pos, Peptide)`` in that pool —
   the same 9-mer at the same ``Pos`` on two different merged parents still counts **twice**;
   under “unique” it counts **once**.

**Two reporting modes (both printed):**

- **Per-parent densities** (+ Mann–Whitney): **exactly 501 vs 501** — one row per entry in
  ``sig_parent_micropeptides.csv`` vs ``coding_proportional_whole_parent_micropeptides.csv``.
  All merged NetMHC SB rows whose ``input_seq_id_fasta`` maps to that MP (``smPEP_ID``)
  or coding ``control_id`` are aggregated; **positional** = distinct
  ``(input_seq_id_fasta, Pos, Peptide)`` (allele rows collapsed); **unique** = distinct
  ``Peptide``. Density = count × 100 / ``length_aa`` from the same parent row. MPs with
  no SB rows contribute **0**.

- **Cohort-level** (printed after): one numerator over **total aa screened** — sum of
  ``length_aa`` across **all** distinct ``input_seq_id_fasta`` in the merged TSV, with
  global dedup of SB peptides or of ``(input_seq_id_fasta, Pos, Peptide)`` slots.
"""
from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DATA_NET = REPO_ROOT / "data" / "netmhc"
DATA = REPO_ROOT / "data"
KNOWN_FA = DATA / "known_proteins.fasta"
DEFAULT_UNIPROT_CACHE = DATA_NET / "uniprot_protein_length_cache.tsv"
sys.path.insert(0, str(HERE.parent / "scripts"))
from netmhc_sb_core import (  # noqa: E402
    FIG5_IEDB_EL_RANK_MAX_DEFAULT,
    FIG5_IEDB_IC50_MAX_NM_DEFAULT,
    ba_score_min_for_ic50_lt,
    pick_iedb_ic50_column,
    sb_mask_fig5_defaults,
)


def cod_control_key(parent: str) -> str:
    p = str(parent)
    if p.startswith("COD|"):
        p = p[4:]
    parts = p.split("|")
    if len(parts) >= 7:
        return "|".join(parts[:6])
    return p


def uniprot_acc_from_parent(parent: str) -> str | None:
    m = re.search(r"sp\|([A-Z0-9]+)\|", str(parent))
    return m.group(1) if m else None


def load_fasta_lengths(fa: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    for rec in SeqIO.parse(fa, "fasta"):
        hid = rec.id
        if hid.startswith("sp|"):
            parts = hid.split("|")
            if len(parts) >= 2:
                out[parts[1]] = len(rec.seq)
    return out


def load_parent_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
    sig = pd.read_csv(DATA_NET / "sig_parent_micropeptides.csv", dtype=str)
    cod = pd.read_csv(DATA_NET / "coding_proportional_whole_parent_micropeptides.csv", dtype=str)
    sm_len = sig.set_index(sig["smPEP_ID"].astype(str))["length_aa"].astype(float)
    cod_len = cod.set_index(cod["control_id"].astype(str))["length_aa"].astype(float)
    return sig, cod, sm_len, cod_len


def read_length_cache(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype={"accession": str, "length": int})
    return dict(zip(df["accession"].str.strip(), df["length"]))


def write_length_cache(path: Path, lengths: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(lengths.items(), key=lambda kv: kv[0])
    pd.DataFrame(rows, columns=["accession", "length"]).to_csv(path, sep="\t", index=False)


def fetch_uniprot_lengths(
    accessions: set[str],
    cache_path: Path,
    *,
    batch_size: int = 80,
    sleep_s: float = 0.2,
) -> dict[str, int]:
    """Merge cache + UniProt REST (Swiss-Prot stream) for missing accessions."""
    lengths = read_length_cache(cache_path)
    missing = sorted(accessions - set(lengths.keys()))
    if not missing:
        return {a: lengths[a] for a in accessions if a in lengths}

    for i in range(0, len(missing), batch_size):
        chunk = missing[i : i + batch_size]
        q = "accession:(" + " OR ".join(chunk) + ")"
        url = (
            "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=tsv"
            "&fields=accession,length&query=" + urllib.parse.quote(q)
        )
        try:
            raw = urllib.request.urlopen(url, timeout=120).read().decode("utf-8")
        except urllib.error.HTTPError as e:
            raise SystemExit(f"UniProt HTTP error for batch starting {chunk[0]!r}: {e}") from e
        lines = raw.strip().splitlines()
        if len(lines) < 2:
            time.sleep(sleep_s)
            continue
        for ln in lines[1:]:
            parts = ln.split("\t")
            if len(parts) >= 2 and parts[0] and parts[1].isdigit():
                lengths[parts[0].strip()] = int(parts[1].strip())
        time.sleep(sleep_s)

    write_length_cache(cache_path, lengths)
    return {a: lengths[a] for a in accessions if a in lengths}


def read_sb_usecols(path: Path) -> set[str]:
    hdr = set(pd.read_csv(path, sep="\t", nrows=0).columns)
    need = {
        "input_seq_id_fasta",
        "Pos",
        "Peptide",
        "EL_rank",
        "iedb_score",
        "iedb_processing_score",
    }
    iedb_ic50 = pick_iedb_ic50_column(sorted(hdr))
    if iedb_ic50:
        need.add(iedb_ic50)
    else:
        need.add("BA_score")
    missing = need - hdr
    if missing:
        raise SystemExit(f"{path}: missing columns {sorted(missing)}")
    return need


def apply_sb(path: Path, sb_mode: str) -> pd.DataFrame:
    usecols = read_sb_usecols(path)
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    hdr = set(df.columns)
    iedb_ic50 = pick_iedb_ic50_column(sorted(hdr))
    ba_min = ba_score_min_for_ic50_lt(FIG5_IEDB_IC50_MAX_NM_DEFAULT)
    mask = sb_mask_fig5_defaults(
        df,
        el_max=FIG5_IEDB_EL_RANK_MAX_DEFAULT,
        el_lte=False,
        ba_min=ba_min,
        ic50_max_nm=FIG5_IEDB_IC50_MAX_NM_DEFAULT,
        iedb_ic50_col=iedb_ic50,
        sb_mode=sb_mode,
    )
    return df.loc[mask].copy()


def sig_metrics_501(
    df_sb: pd.DataFrame,
    sig_parent: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, int, dict[str, int]]:
    """One density per ``sig_parent_micropeptides`` row: pool SB rows by ``smPEP_ID`` (second ``|`` field)."""
    u_out: list[float] = []
    w_out: list[float] = []
    skipped = 0
    n_with_sb = 0
    grouped: dict[str, pd.DataFrame] = {}
    if len(df_sb):
        gsb = df_sb.copy()
        gsb["smpep"] = gsb["input_seq_id_fasta"].astype(str).str.split("|").str[1]
        gsb["pep"] = gsb["Peptide"].astype(str).str.upper()
        grouped = {str(k): g for k, g in gsb.groupby("smpep", sort=False)}
    for _, row in sig_parent.iterrows():
        smpep = str(row["smPEP_ID"])
        L = float(row["length_aa"])
        if not np.isfinite(L) or L <= 0:
            skipped += 1
            u_out.append(np.nan)
            w_out.append(np.nan)
            continue
        g = grouped.get(smpep)
        if g is None or len(g) == 0:
            nu, nw = 0, 0
        else:
            n_with_sb += 1
            nu = int(g["pep"].nunique())
            nw = int(g.groupby(["input_seq_id_fasta", "Pos", "pep"], sort=False).ngroups)
        u_out.append(100.0 * nu / L)
        w_out.append(100.0 * nw / L)
    n_par = len(u_out)
    stats = {"n_parents": n_par, "n_with_sb": n_with_sb, "n_zero_sb": n_par - n_with_sb}
    return np.asarray(u_out), np.asarray(w_out), skipped, stats


def cod_metrics_501(
    df_sb: pd.DataFrame,
    cod_parent: pd.DataFrame,
    len_by_acc: dict[str, int],
    *,
    length_source: str,
) -> tuple[np.ndarray, np.ndarray, int, dict[str, int]]:
    """One density per ``coding_proportional_whole_parent_micropeptides`` row (``control_id``)."""
    u_out: list[float] = []
    w_out: list[float] = []
    skipped = 0
    n_with_sb = 0
    grouped: dict[str, pd.DataFrame] = {}
    if len(df_sb):
        gsb = df_sb.copy()
        gsb["cid"] = gsb["input_seq_id_fasta"].map(cod_control_key)
        gsb["pep"] = gsb["Peptide"].astype(str).str.upper()
        grouped = {str(k): g for k, g in gsb.groupby("cid", sort=False)}
    for _, row in cod_parent.iterrows():
        cid = str(row["control_id"])
        if length_source == "fragment":
            L = float(row["length_aa"])
        else:
            acc = uniprot_acc_from_parent(cid)
            if not acc:
                skipped += 1
                u_out.append(np.nan)
                w_out.append(np.nan)
                continue
            L = float(len_by_acc.get(acc, np.nan))
        if not np.isfinite(L) or L <= 0:
            skipped += 1
            u_out.append(np.nan)
            w_out.append(np.nan)
            continue
        g = grouped.get(cid)
        if g is None or len(g) == 0:
            nu, nw = 0, 0
        else:
            n_with_sb += 1
            nu = int(g["pep"].nunique())
            nw = int(g.groupby(["input_seq_id_fasta", "Pos", "pep"], sort=False).ngroups)
        u_out.append(100.0 * nu / L)
        w_out.append(100.0 * nw / L)
    n_par = len(u_out)
    stats = {"n_parents": n_par, "n_with_sb": n_with_sb, "n_zero_sb": n_par - n_with_sb}
    return np.asarray(u_out), np.asarray(w_out), skipped, stats


def summarize(label: str, x: np.ndarray, y: np.ndarray, alt: str) -> None:
    xf = x[np.isfinite(x)]
    yf = y[np.isfinite(y)]
    print(f"\n=== {label} ===")
    print(
        f"sig lnc: n={len(xf)} (finite) mean={float(np.mean(xf)):.6g} median={float(np.median(xf)):.6g}"
    )
    print(
        f"coding:  n={len(yf)} (finite) mean={float(np.mean(yf)):.6g} median={float(np.median(yf)):.6g}"
    )
    if len(xf) < 2 or len(yf) < 2:
        print("Mann–Whitney skipped (need at least 2 finite values per cohort).")
        return
    res = mannwhitneyu(yf, xf, alternative=alt)
    print(f"Mann–Whitney coding vs sig, alternative={alt!r}: p={float(res.pvalue):.6g}")
    res2 = mannwhitneyu(xf, yf, alternative="two-sided")
    print(f"Mann–Whitney two-sided: p={float(res2.pvalue):.6g}")


def cohort_scanned_aa_sig(sig_path: Path, sm_len: pd.Series) -> float:
    """Sum ``length_aa`` over every distinct ``input_seq_id_fasta`` in the merged sig TSV."""
    parents = pd.read_csv(sig_path, sep="\t", usecols=["input_seq_id_fasta"])["input_seq_id_fasta"].unique()
    tot = 0.0
    for parent in parents:
        smpep = str(parent).split("|")[1]
        L = float(sm_len.get(str(smpep), np.nan))
        if np.isfinite(L) and L > 0:
            tot += L
    return tot


def cohort_scanned_aa_cod(
    cod_path: Path,
    cod_len: pd.Series,
    len_by_acc: dict[str, int],
    *,
    length_source: str,
) -> float:
    parents = pd.read_csv(cod_path, sep="\t", usecols=["input_seq_id_fasta"])["input_seq_id_fasta"].unique()
    tot = 0.0
    for parent in parents:
        if length_source == "fragment":
            key = cod_control_key(parent)
            L = float(cod_len.get(key, np.nan))
        else:
            acc = uniprot_acc_from_parent(parent)
            if not acc:
                continue
            L = float(len_by_acc.get(acc, np.nan))
        if np.isfinite(L) and L > 0:
            tot += L
    return tot


def print_cohort_summary(
    label: str,
    df_sb: pd.DataFrame,
    sum_aa_scanned: float,
) -> None:
    df_sb = df_sb.copy()
    df_sb["pep"] = df_sb["Peptide"].astype(str).str.upper()
    u = int(df_sb["pep"].nunique())
    pos = int(df_sb.groupby(["input_seq_id_fasta", "Pos", "pep"], sort=False).ngroups)
    print(f"\n--- Cohort-level ({label}) ---")
    print(f"Total aa scanned (sum of parent length_aa over merged parents): {sum_aa_scanned:.6g}")
    print(f"Global unique SB 9-mers (dedup Peptide): {u}")
    print(f"Global positional SB slots (dedup input_seq_id_fasta, Pos, Peptide): {pos}")
    if sum_aa_scanned > 0:
        print(f"Unique SB 9-mers per 100 aa (cohort): {100.0 * u / sum_aa_scanned:.6g}")
        print(f"Positional SB slots per 100 aa (cohort): {100.0 * pos / sum_aa_scanned:.6g}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--sb-mode",
        choices=("full", "no_ic50", "ic50_only"),
        default="full",
        help="SB composition (default: full Fig 5 stack).",
    )
    ap.add_argument(
        "--coding-length-source",
        choices=("fragment", "fasta", "uniprot"),
        default="fragment",
        help="Coding denominator: Fig 5 short parent length_aa (default), bundled FASTA, or UniProt.",
    )
    ap.add_argument(
        "--uniprot-cache",
        type=Path,
        default=DEFAULT_UNIPROT_CACHE,
        help="TSV cache of accession\\tlength for UniProt mode.",
    )
    ap.add_argument(
        "--known-fasta",
        type=Path,
        default=KNOWN_FA,
        help="FASTA used when --coding-length-source fasta.",
    )
    ap.add_argument(
        "--mw-alt",
        choices=("less", "greater", "two-sided"),
        default="less",
        help="Primary Mann–Whitney alternative for coding vs sig (default: less).",
    )
    args = ap.parse_args()

    sig_parent, cod_parent, sm_len, cod_len = load_parent_tables()
    sig_path = DATA_NET / "netmhcpan_sig_lnc_with_iedb.tsv"
    cod_path = DATA_NET / "netmhcpan_coding_proportional_whole_with_iedb.tsv"

    print(f"SB mode: {args.sb_mode!r}")
    print(
        "Cohorts: Fig 5 merged TSVs - sig lnc vs coding proportional whole "
        "(short sampled parents; not random genomic fragments)."
    )
    print(f"Coding length source: {args.coding_length_source!r}")

    sig_sb = apply_sb(sig_path, args.sb_mode)
    cod_sb = apply_sb(cod_path, args.sb_mode)
    print(f"SB rows retained: sig {len(sig_sb)} / coding {len(cod_sb)}")

    len_by_acc: dict[str, int] = {}
    if args.coding_length_source == "uniprot":
        accs: set[str] = set()
        for p in pd.read_csv(cod_path, sep="\t", usecols=["input_seq_id_fasta"])["input_seq_id_fasta"].unique():
            a = uniprot_acc_from_parent(str(p))
            if a:
                accs.add(a)
        len_by_acc = fetch_uniprot_lengths(accs, args.uniprot_cache)
        got = len(len_by_acc)
        print(f"UniProt lengths resolved for {got} / {len(accs)} distinct accessions (see {args.uniprot_cache})")
    elif args.coding_length_source == "fasta":
        len_by_acc = load_fasta_lengths(args.known_fasta)
        print(f"FASTA lengths from {args.known_fasta}")
    else:
        print(
            "Denominator (Fig 5 default): parent length_aa from "
            "sig_parent_micropeptides + coding_proportional_whole_parent_micropeptides."
        )

    su, sw, skip_sig, st_sig = sig_metrics_501(sig_sb, sig_parent)
    cu, cw, skip_cod, st_cod = cod_metrics_501(
        cod_sb, cod_parent, len_by_acc, length_source=args.coding_length_source
    )
    print(f"Sig parents skipped (missing length_aa): {skip_sig}")
    print(f"Coding parents skipped (missing length): {skip_cod}")
    print(
        "Per-parent universe (Mann–Whitney): sig n="
        f"{st_sig['n_parents']} (SB>0: {st_sig['n_with_sb']}, zero SB: {st_sig['n_zero_sb']}); "
        f"coding n={st_cod['n_parents']} (SB>0: {st_cod['n_with_sb']}, zero SB: {st_cod['n_zero_sb']})"
    )

    sum_sig_aa = cohort_scanned_aa_sig(sig_path, sm_len)
    sum_cod_aa = cohort_scanned_aa_cod(cod_path, cod_len, len_by_acc, length_source=args.coding_length_source)
    print(
        "\nCohort totals: denominator = sum of parent length_aa over **all** distinct "
        "input_seq_id_fasta rows in each merged TSV (total aa screened), "
        "not restricted to parents with SB."
    )
    print_cohort_summary("sig lnc", sig_sb, sum_sig_aa)
    print_cohort_summary("coding proportional whole", cod_sb, sum_cod_aa)

    m_sig = np.isfinite(sw) & np.isfinite(su)
    m_cod = np.isfinite(cw) & np.isfinite(cu)
    if np.any((sw > su) & m_sig):
        n = int(np.sum((sw > su) & m_sig))
        print(f"\nNote: {n} sig MPs (501 table) have positional density > unique.")
    else:
        print(
            "\nNote (501 sig MPs): under this SB filter, positional density did not exceed "
            "unique density on any MP (often true when each screened string yields at most "
            "one row per (Pos, peptide) under the gate)."
        )

    if np.any((cw > cu) & m_cod):
        n2 = int(np.sum((cw > cu) & m_cod))
        print(f"Note: {n2} coding controls (501 table) have positional density > unique.")

    summarize(
        "Unique SB 9-mers per 100 aa (501 vs 501 parent tables; zero if no SB)",
        su,
        cu,
        args.mw_alt,
    )
    summarize(
        "Positional SB per 100 aa (501 vs 501 parent tables; zero if no SB)",
        sw,
        cw,
        args.mw_alt,
    )


if __name__ == "__main__":
    main()
