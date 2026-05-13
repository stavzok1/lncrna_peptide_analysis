"""
Incremental NetMHCpan helpers for the significant-lncRNA 9-mer track.

Compare current ``data/significant_lnc_peptides.tsv`` + ``significant_lnc_peptides.faa``
to the last full NetMHC prep snapshot ``data/netmhc/sig_parent_micropeptides.csv``.

Commands
~~~~~~~~

**strip** — Write ``data/netmhc/netmhcpan_sig_lnc_stripped.xls``: same first three lines
as ``netmhcpan_sig_lnc.xls``, but data rows whose NetMHC ``ID`` encodes a **removed**
``smPEP_ID`` (regex ``SIG_<digits>_`` on column 3) are dropped. Also writes
``data/netmhc/ninemers_sig_lnc_added_only.fasta`` for parents in TSV+FASTA today but
absent from ``sig_parent_micropeptides.csv`` (run NetMHCpan on this in WSL).

**merge** — After ``netMHCpan ... -xlsfile data/netmhc/netmhcpan_sig_lnc_added_only.xls``,
concatenate ``netmhcpan_sig_lnc_stripped.xls`` with **data** lines from the added-only XLS
(skip its first three lines) into ``netmhcpan_sig_lnc_merged.xls``.

**Caveat:** coding controls stay from the old prep unless you rerun
``prepare_netmhc_tr_vs_coding_epitopes.py`` end-to-end (preferred for strict paired design).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
NET = DATA / "netmhc"
SIG_TSV = DATA / "significant_lnc_peptides.tsv"
SIG_FAA = DATA / "significant_lnc_peptides.faa"
OLD_PARENTS = NET / "sig_parent_micropeptides.csv"
SIG_XLS = NET / "netmhcpan_sig_lnc.xls"
OUT_STRIPPED = NET / "netmhcpan_sig_lnc_stripped.xls"
OUT_ADDED_FA = NET / "ninemers_sig_lnc_added_only.fasta"
OUT_MERGED = NET / "netmhcpan_sig_lnc_merged.xls"

_ID_RE = re.compile(r"^SIG_(\d+)_")


def smpep_from_netmhc_id(field: str) -> str | None:
    m = _ID_RE.match(str(field).strip())
    return m.group(1) if m else None


def diff_sets() -> tuple[set[str], set[str], set[str], set[str]]:
    cur = set(pd.read_csv(SIG_TSV, sep="\t", dtype=str)["smPEP_ID"].str.strip())
    if not OLD_PARENTS.exists():
        raise SystemExit(f"Missing prior parent table: {OLD_PARENTS}")
    old = set(pd.read_csv(OLD_PARENTS, dtype=str)["smPEP_ID"].astype(str).str.strip())
    return cur, old, cur - old, old - cur


def cmd_strip(args: argparse.Namespace) -> None:
    import prepare_netmhc_tr_vs_coding_epitopes as prep

    xls_in = Path(args.xls_in)
    if not xls_in.exists():
        raise SystemExit(f"Missing NetMHCpan XLS: {xls_in}")
    if not SIG_FAA.exists():
        raise SystemExit(f"Missing {SIG_FAA}")

    cur, old, added, removed = diff_sets()
    removed = set(removed)
    added = set(added)

    lines = xls_in.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 4:
        raise SystemExit(f"Unexpected short XLS: {xls_in}")
    header = lines[:3]
    kept: list[str] = []
    dropped = 0
    no_id = 0
    for ln in lines[3:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 3:
            kept.append(ln)
            continue
        sid = smpep_from_netmhc_id(parts[2])
        if sid is None:
            no_id += 1
            kept.append(ln)
            continue
        if sid in removed:
            dropped += 1
            continue
        kept.append(ln)

    out_xls = Path(args.out_stripped)
    out_xls.parent.mkdir(parents=True, exist_ok=True)
    out_xls.write_text("\n".join(header + kept) + "\n", encoding="utf-8")

    by_id = {s: (sym, seq) for s, sym, seq in prep.parse_sig_peptide_fasta(SIG_FAA)}
    sig_9: list[tuple[str, str]] = []
    for pid in sorted(added, key=lambda x: int(x)):
        if pid not in by_id:
            print(f"Warning: added smPEP_ID {pid} not in FASTA — no 9-mers emitted.", flush=True)
            continue
        sym, seq = by_id[pid]
        mer_id = f"{pid}|{sym}"
        sig_9.extend(prep.sliding_ninemers(seq, mer_id, "SIG"))
    prep.write_fasta(Path(args.out_added_fa), sig_9)

    print(
        f"Current TSV: {len(cur)} | Old parent CSV: {len(old)} | "
        f"Added (need new NetMHCpan): {len(added)} | Removed (stripped from XLS): {len(removed)}"
    )
    print(f"XLS rows dropped (matched removed smPEP_ID): {dropped}")
    print(f"XLS rows kept without parseable SIG_<id>_: {no_id}")
    print(f"Wrote {out_xls}")
    print(f"Wrote {args.out_added_fa} ({len(sig_9)} 9-mer records)")


def cmd_merge(args: argparse.Namespace) -> None:
    stripped = Path(args.stripped_xls)
    added_xls = Path(args.added_xls)
    if not stripped.exists():
        raise SystemExit(f"Missing {stripped} (run strip first)")
    if not added_xls.exists():
        raise SystemExit(f"Missing {added_xls}")

    base_lines = stripped.read_text(encoding="utf-8", errors="replace").splitlines()
    add_lines = added_xls.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(add_lines) < 4:
        raise SystemExit(f"Unexpected short added XLS: {added_xls}")

    merged = base_lines + add_lines[3:]
    out = Path(args.out_merged)
    out.write_text("\n".join(merged) + "\n", encoding="utf-8")
    print(f"Wrote {out} ({len(merged)} lines)")


def main() -> None:
    ap = argparse.ArgumentParser(description="Strip / merge NetMHCpan sig lnc XLS incrementally.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("strip", help="Drop removed-parent rows; write FASTA for added parents.")
    sp.add_argument("--xls-in", type=Path, default=SIG_XLS)
    sp.add_argument("--out-stripped", type=Path, default=OUT_STRIPPED)
    sp.add_argument("--out-added-fa", type=Path, default=OUT_ADDED_FA)
    sp.set_defaults(func=cmd_strip)

    mg = sub.add_parser("merge", help="Append added-only NetMHCpan XLS data to stripped base.")
    mg.add_argument("--stripped-xls", type=Path, default=OUT_STRIPPED)
    mg.add_argument("--added-xls", type=Path, default=NET / "netmhcpan_sig_lnc_added_only.xls")
    mg.add_argument("--out-merged", type=Path, default=OUT_MERGED)
    mg.set_defaults(func=cmd_merge)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
