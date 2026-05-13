#!/usr/bin/env python3
"""
Merge wide NetMHCpan-4.2 ``*.xls`` (tab-separated) outputs with IEDB Next-Gen CSV
(``iedb_*_el_tap_immuno.csv``) by adding IEDB columns to a long per-allele table.

Join strategies
---------------
``stable`` (default): Build the same ``stable_key`` as the IEDB runner
(``input_seq_id:start:end:peptide:allele``). Wide rows are aligned to FASTA
records **in file order** (``--ninemers-fasta``): NetMHCpan often writes
``n_wide / n_fasta`` rows per record (for example a 9-mer plus an ``X``-padded
frame).

If IEDB was run on the **full parent** peptide FASTA (e.g. ``significant_lnc_peptides.faa``)
but NetMHC was run on ``SIG|smPEP|gene|k|...`` ninemers, pass
``--parent-peptide-fasta`` pointing at that parent FASTA. The script then maps
each ninemer header to the parent ``input_seq_id`` (same numeric smPEP id as
the parent's first ``|`` field) and 1-based ``start/end`` from the ninemer
offset ``k``, and merges on that IEDB ``stable_key``. NetMHC ``X``-padded rows
are skipped. ``COD|coding|sp|...`` ninemers (proteome control) map to parent FASTA
ids ``coding|sp|...`` built from header fields 1–7, with window index in field 8.

``smpep_window``: Join on
``(smPEP_id, start, end, peptide, allele)`` where NetMHC ``smPEP_id`` is parsed
from ``ID`` (``SIG_<digits>_...``) and ``start`` is the ``Pos`` column. IEDB
``smPEP_id`` is the first ``|``-separated field of ``input_seq_id``. This only
aligns rows when both pipelines use the **same** numeric smPEP ids and the same
1-based window coordinates.

IEDB columns prefixed with ``iedb_`` are merged in. With ``--skip-iedb-binding``,
IEDB ``netmhcpan_el_*`` and most ``netmhcpan_ba_*`` fields are dropped so they do
not duplicate local NetMHC wide EL/BA — **``netmhcpan_ba_ic50`` (nM) is kept**
when present so downstream filters can use IEDB-reported IC50.

**``--synthetic-iedb-pass``** (omit ``--iedb-csv``): skip the IEDB join; add
``iedb_score`` / ``iedb_processing_score`` constants that pass default Fig 5 gates
and ``iedb_netmhcpan_ba_ic50`` from NetMHC ``BA_score`` (TTN-companion-style mapping).
For proportional-whole coding ninemers (``COD|coding_pwhole|…``) without an IEDB run.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def normalize_hla_netmhc_to_iedb(name: str) -> str:
    """HLA-A01:01 -> HLA-A*01:01 (IEDB style). Other tokens unchanged."""
    s = name.strip()
    if not s:
        return s
    up = s.upper()
    m = re.fullmatch(r"HLA-([ABCEFG])(\d{2}:\d{2})", up)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return s


def parse_netmhc_sig_smpep(seq_id: str) -> str | None:
    m = re.match(r"^SIG_(\d+)_", seq_id.strip())
    return m.group(1) if m else None


def parse_fasta_records(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    cur_id: str | None = None
    cur_seq: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                seq = "".join(cur_seq).strip()
                if seq:
                    records.append((cur_id, seq))
            cur_id = line[1:].split()[0] if len(line) > 1 else "seq"
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id is not None:
        seq = "".join(cur_seq).strip()
        if seq:
            records.append((cur_id, seq))
    return records


@dataclass(frozen=True)
class ParentPeptideIndex:
    """SIG parents: lookup by numeric smPEP (first header field). COD parents: full header id."""

    by_numeric_smpep: dict[str, tuple[str, str]]
    by_full_id: dict[str, tuple[str, str]]


def load_parent_peptide_index(parent_fasta: Path) -> ParentPeptideIndex:
    """
    Index parent FASTA: full header id for all records; additionally ``by_numeric_smpep``
    when the first ``|`` field is all digits (lnc smPEP ids).
    """
    by_full: dict[str, tuple[str, str]] = {}
    by_num: dict[str, tuple[str, str]] = {}
    for rid, seq in parse_fasta_records(_read_text(parent_fasta)):
        by_full[rid] = (rid, seq)
        head = rid.split("|", 1)[0]
        if re.fullmatch(r"\d+", head):
            if head in by_num and by_num[head][1] != seq:
                raise SystemExit(
                    f"Duplicate numeric parent id {head!r} with different sequences in {parent_fasta}"
                )
            by_num[head] = (rid, seq)
    return ParentPeptideIndex(by_numeric_smpep=by_num, by_full_id=by_full)


def crosswalk_stable_key_from_ninemer_row(
    ninemer_header_id: str,
    peptide: str,
    allele_iedb: str,
    parent_index: ParentPeptideIndex,
) -> str | None:
    """
    Map ninemer header + 9-mer to IEDB ``stable_key`` (SIG lnc or COD proteome).
    """
    parts = ninemer_header_id.split("|")
    pep = str(peptide).strip()
    if len(pep) != 9 or "X" in pep.upper():
        return None
    al = str(allele_iedb).strip()

    if len(parts) >= 5 and parts[0].upper() == "SIG":
        sm = parts[1]
        try:
            k0 = int(parts[3])
        except ValueError:
            return None
        row = parent_index.by_numeric_smpep.get(sm)
        if not row:
            return None
        sid, pseq = row
        if k0 < 0 or k0 + 9 > len(pseq) or pseq[k0 : k0 + 9] != pep:
            return None
        return f"{sid}:{k0 + 1}:{k0 + 9}:{pep}:{al}"

    if len(parts) >= 10 and parts[0].upper() == "COD":
        # COD|coding|sp|UNIPROT|GENE|...|pool|k|ref_peptide — parent id = fields 1..7
        parent_sid = "|".join(parts[1:8])
        try:
            k0 = int(parts[8])
        except ValueError:
            return None
        row = parent_index.by_full_id.get(parent_sid)
        if not row:
            return None
        sid, pseq = row
        if k0 < 0 or k0 + 9 > len(pseq) or pseq[k0 : k0 + 9] != pep:
            return None
        return f"{sid}:{k0 + 1}:{k0 + 9}:{pep}:{al}"

    return None


def add_crosswalk_stable_key_column(
    df: pd.DataFrame,
    parent_index: ParentPeptideIndex,
) -> pd.DataFrame:
    h = df["input_seq_id_fasta"].astype(str).to_numpy()
    p = df["Peptide"].astype(str).to_numpy()
    a = df["allele"].astype(str).to_numpy()
    keys: list[str | None] = [
        crosswalk_stable_key_from_ninemer_row(h[i], p[i], a[i], parent_index)
        for i in range(len(df))
    ]
    out = df.copy()
    out["crosswalk_stable_key"] = keys
    return out


def parse_allele_header_row(line: str) -> list[str]:
    toks = line.split("\t")
    alleles: list[str] = []
    for t in toks:
        t = t.strip()
        if t.startswith("HLA-"):
            if not alleles or alleles[-1] != t:
                alleles.append(t)
    return alleles


def wide_xls_to_long(path: Path) -> pd.DataFrame:
    lines = _read_text(path).splitlines()
    hi = next(i for i, ln in enumerate(lines) if ln.startswith("Pos\t"))
    allele_row = lines[hi - 1] if hi > 0 else ""
    alleles = parse_allele_header_row(allele_row)
    if not alleles:
        raise SystemExit(f"No HLA alleles parsed from line before header in {path}")
    n_block = 10  # core,icore,EL_score,EL_rank,BA_score,BA_rank,Pathogen,Pathogen_rank,Neo,Neo_rank
    block_cols = [
        "core",
        "icore",
        "EL_score",
        "EL_rank",
        "BA_score",
        "BA_rank",
        "Pathogen_score",
        "Pathogen_rank",
        "Neo_score",
        "Neo_rank",
    ]
    rows_out: list[dict[str, Any]] = []
    wide_row_ix = 0
    warned_trunc = False
    for ln in lines[hi + 1 :]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 3 + n_block:
            continue
        n_blocks = (len(parts) - 3) // n_block
        if n_blocks < 1:
            continue
        if n_blocks < len(alleles) and not warned_trunc:
            print(
                f"[merge_netmhcpan_iedb] warning: {path.name} rows have {n_blocks} allele block(s) "
                f"but header lists {len(alleles)}; trailing alleles may be missing (e.g. column cap).",
                file=sys.stderr,
            )
            warned_trunc = True
        try:
            pos = int(float(parts[0]))
        except ValueError:
            continue
        peptide = parts[1]
        seq_id = parts[2]
        for ai in range(min(len(alleles), n_blocks)):
            allele = alleles[ai]
            j = 3 + ai * n_block
            chunk = parts[j : j + n_block]
            if len(chunk) < n_block:
                continue
            row: dict[str, Any] = {
                "wide_row_ix": wide_row_ix,
                "Pos": pos,
                "Peptide": peptide,
                "ID": seq_id,
                "allele_netmhc": allele,
                "allele": normalize_hla_netmhc_to_iedb(allele),
            }
            for name, val in zip(block_cols, chunk, strict=True):
                row[name] = val
            rows_out.append(row)
        wide_row_ix += 1
    return pd.DataFrame(rows_out)


def attach_stable_keys_from_fasta_order(df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    NetMHCpan may emit multiple wide rows per FASTA record (e.g. real 9-mer + X-padded
    frame). Map ``wide_row_ix // (n_wide / n_fasta)`` to FASTA record order and build
    ``stable_key`` like the IEDB runner: ``input_seq_id:start:end:peptide:allele``.
    """
    recs = parse_fasta_records(_read_text(fasta_path))
    n_fasta = len(recs)
    if n_fasta == 0:
        raise SystemExit(f"No sequences in {fasta_path}")
    n_wide = int(df["wide_row_ix"].max()) + 1 if len(df) else 0
    if n_wide % n_fasta != 0:
        raise SystemExit(
            f"Wide row count {n_wide} is not divisible by FASTA records {n_fasta} "
            f"({fasta_path}); cannot align by order."
        )
    rpf = n_wide // n_fasta
    rids = [rid for rid, _ in recs]

    def stable_for_row(i: int) -> str:
        wi = int(df["wide_row_ix"].iloc[i])
        rid = rids[wi // rpf]
        pos = int(df["Pos"].iloc[i])
        pep = str(df["Peptide"].iloc[i])
        end = pos + len(pep) - 1
        al = str(df["allele"].iloc[i])
        return f"{rid}:{pos}:{end}:{pep}:{al}"

    out = df.copy()
    out["input_seq_id_fasta"] = [rids[int(df["wide_row_ix"].iloc[i]) // rpf] for i in range(len(df))]
    out["stable_key"] = [stable_for_row(i) for i in range(len(df))]
    return out


def smpep_window_join_key(df: pd.DataFrame) -> pd.Series:
    def row_key(i: int) -> str | None:
        sp = parse_netmhc_sig_smpep(str(df["ID"].iloc[i]))
        if not sp:
            return None
        s = int(df["Pos"].iloc[i])
        pep = str(df["Peptide"].iloc[i])
        e = s + len(pep) - 1
        al = str(df["allele"].iloc[i])
        return f"{sp}:{s}:{e}:{pep}:{al}"

    return pd.Series([row_key(i) for i in range(len(df))], index=df.index)


def iedb_smpep_window_key_from_row(r: pd.Series) -> str:
    sid = str(r["iedb_input_seq_id"]).split("|", 1)[0]
    return (
        f"{sid}:{int(r['iedb_start'])}:{int(r['iedb_end'])}:"
        f"{r['iedb_peptide']}:{r['iedb_allele']}"
    )


def load_iedb_csv(path: Path, skip_binding: bool) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    drop = {"sequence_number", "peptide_index", "median_percentile"}
    if skip_binding:
        for n in df.columns:
            if n.startswith("netmhcpan_el_"):
                drop.add(n)
            elif n.startswith("netmhcpan_ba_"):
                # Keep IEDB BA IC50 (nM) for binding filters; drop other BA fields (core/icore/percentile).
                if "ic50" not in n.lower():
                    drop.add(n)
    use = [c for c in df.columns if c not in drop]
    out = df[use].copy()
    rename = {c: f"iedb_{c}" for c in out.columns if not c.startswith("iedb_")}
    out = out.rename(columns=rename)
    return out


def ba_score_to_ic50_nm(ba: float) -> float:
    x = float(ba)
    if not math.isfinite(x):
        return float("nan")
    x = min(1.0 - 1e-12, max(1e-12, x))
    return float(50000.0 ** (1.0 - x))


def add_synthetic_iedb_pass_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Minimal IEDB-like columns so ``sb_mask_spec`` full mode can use BA-derived IC50."""
    out = df.copy()
    out["iedb_score"] = 0.25
    out["iedb_processing_score"] = 2.0
    ba = pd.to_numeric(out["BA_score"], errors="coerce")
    out["iedb_netmhcpan_ba_ic50"] = [ba_score_to_ic50_nm(float(x)) if pd.notna(x) else float("nan") for x in ba]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--netmhc-xls", type=Path, required=True, help="Wide NetMHCpan .xls (tab-separated).")
    ap.add_argument(
        "--iedb-csv",
        type=Path,
        default=None,
        help="IEDB peptide_table CSV with stable_key (omit with --synthetic-iedb-pass).",
    )
    ap.add_argument("--out-tsv", type=Path, required=True, help="Output long TSV path.")
    ap.add_argument(
        "--ninemers-fasta",
        type=Path,
        help="Ninemers FASTA used for NetMHC (required for --join stable).",
    )
    ap.add_argument(
        "--parent-peptide-fasta",
        type=Path,
        help="Full peptide FASTA used for IEDB (e.g. significant_lnc_peptides.faa). "
        "When set with --join stable, maps SIG ninemers to parent stable_key for the merge.",
    )
    ap.add_argument(
        "--join",
        choices=("stable", "smpep_window"),
        default="stable",
        help="stable: merge on stable_key (see --parent-peptide-fasta). "
        "smpep_window: merge on first field of input_seq_id + start/end/peptide/allele.",
    )
    ap.add_argument(
        "--skip-iedb-binding",
        action="store_true",
        help="Drop IEDB netmhcpan_el_* and netmhcpan_ba_* except IC50 (keep netmhcpan_ba_ic50 + TAP/immuno/processing).",
    )
    ap.add_argument(
        "--synthetic-iedb-pass",
        action="store_true",
        help="Skip IEDB CSV; add synthetic imm/proc + BA→IC50 columns (see module doc).",
    )
    args = ap.parse_args()

    if args.synthetic_iedb_pass and args.iedb_csv is not None:
        ap.error("Do not pass --iedb-csv with --synthetic-iedb-pass")
    if not args.synthetic_iedb_pass and args.iedb_csv is None:
        ap.error("Provide --iedb-csv or use --synthetic-iedb-pass")
    if args.synthetic_iedb_pass and args.join != "stable":
        ap.error("--synthetic-iedb-pass requires --join stable")
    if args.synthetic_iedb_pass and not args.ninemers_fasta:
        ap.error("--synthetic-iedb-pass requires --ninemers-fasta")

    long_df = wide_xls_to_long(args.netmhc_xls)
    unmatched = 0

    if args.synthetic_iedb_pass:
        long_df = attach_stable_keys_from_fasta_order(long_df, args.ninemers_fasta)
        merged = add_synthetic_iedb_pass_columns(long_df)
        print(
            f"[merge_netmhcpan_iedb] synthetic_iedb_pass rows={len(merged)} (no IEDB join)",
            file=sys.stderr,
        )
    else:
        iedb = load_iedb_csv(args.iedb_csv, skip_binding=args.skip_iedb_binding)

        if args.join == "stable":
            if not args.ninemers_fasta:
                ap.error("--join stable requires --ninemers-fasta")
            long_df = attach_stable_keys_from_fasta_order(long_df, args.ninemers_fasta)
            if args.parent_peptide_fasta:
                parent_ix = load_parent_peptide_index(args.parent_peptide_fasta)
                long_df = add_crosswalk_stable_key_column(long_df, parent_ix)
                merged = long_df.merge(
                    iedb,
                    left_on="crosswalk_stable_key",
                    right_on="iedb_stable_key",
                    how="left",
                )
                merged = merged.drop(columns=["iedb_stable_key"], errors="ignore")
                n_xw = int(long_df["crosswalk_stable_key"].notna().sum())
                print(
                    f"[merge_netmhcpan_iedb] crosswalk_stable_key non-null rows={n_xw} / {len(long_df)}",
                    file=sys.stderr,
                )
            else:
                merged = long_df.merge(
                    iedb,
                    left_on="stable_key",
                    right_on="iedb_stable_key",
                    how="left",
                )
                merged = merged.drop(columns=["iedb_stable_key"], errors="ignore")
        else:
            long_df["_jw"] = smpep_window_join_key(long_df)
            iedb = iedb.copy()
            iedb["_jw"] = iedb.apply(iedb_smpep_window_key_from_row, axis=1)
            unmatched = int(long_df["_jw"].isna().sum())
            merged = long_df.merge(iedb, on="_jw", how="left")
            merged = merged.drop(columns=["_jw"], errors="ignore")

        n = len(merged)
        if "iedb_total_score" in merged.columns:
            got = int(merged["iedb_total_score"].notna().sum())
        else:
            icols = [c for c in merged.columns if c.startswith("iedb_")]
            got = int(merged[icols].notna().any(axis=1).sum()) if icols else 0
        print(
            f"[merge_netmhcpan_iedb] rows={n} join={args.join!r} "
            f"rows_with_iedb_total_score={int(got)} netmhc_missing_join_key={int(unmatched)}",
            file=sys.stderr,
        )

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out_tsv, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Wrote {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
