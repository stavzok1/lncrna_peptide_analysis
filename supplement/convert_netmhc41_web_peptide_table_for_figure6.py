#!/usr/bin/env python3
"""
Convert **NetMHCpan 4.1 web** ``peptide_table*.tsv`` (long format) into inputs for
``plot_figure6_ttn_as1_allele_coverage.py`` with ``--gating iedb_sb``:

1. **Synthetic wide XLS** (same 10-column-per-allele layout as local ``-xls 1`` exports).
2. **Companion CSV** with ``stable_key`` plus columns that merge to ``iedb_score``,
   ``iedb_processing_score``, and ``iedb_netmhcpan_ba_ic50``. Processing values prefer the
   web column **``processing score``** (falls back to ``processing total score``), matching
   strict web-style filters such as
   ``immunogenicity score > 0.1`` & ``processing score > 1.5`` & ``netmhcpan_el percentile <= 2``.

``stable_key`` format matches ``plot_figure6_ttn_as1_allele_coverage.build_ttn_long_for_iedb_merge``::

    {input_seq_id}:{start1}:{end1}:{PEPTIDE}:{HLA-A*01:01}

where ``start1`` / ``end1`` are **1-based inclusive** coordinates on the TTN-AS1 parent,
``PEPTIDE`` is the 9-mer, and the allele token uses an asterisk (e.g. ``HLA-A*01:01``).

**Defaults:** allele order from ``data/netmhc/hla_european27_class1.txt``;
``--input-seq-id`` default ``108065`` (must match ``--iedb-parent-input-seq-id`` on the plot).

Example::

    python scripts/convert_netmhc41_web_peptide_table_for_figure6.py \\
        --web-tsv data/netmhc/predictions/ttn_as1_web_netmhc41/peptide_table_a08440bd.tsv \\
        --out-dir data/netmhc/predictions/ttn_as1_web_netmhc41

Then Fig 6 (see ``figures/manuscript_netmhc/fig6/ttn_web_netmhc41/README.txt`` for ``sb_full`` vs
``sb_ic50_lt150nm`` folders). Example (full SB)::

    python plot_figure6_ttn_as1_allele_coverage.py --gating iedb_sb --sb-mode full --ic50-max-nm 150 \\
        --netmhc-xls data/netmhc/predictions/ttn_as1_web_netmhc41/netmhcpan_ttn_as1_web_netmhc41.xls \\
        --iedb-csv data/netmhc/predictions/ttn_as1_web_netmhc41/iedb_companion_for_fig6.csv \\
        --iedb-parent-input-seq-id 108065 --split-panels --also-write-unique \\
        -o figures/manuscript_netmhc/fig6/ttn_web_netmhc41/sb_full/fig6_split.png
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


import argparse
import math
import re

import pandas as pd

def xls_allele_token(line: str) -> str:
    s = line.strip().upper().replace("*", "")
    if not s.startswith("HLA-"):
        s = "HLA-" + s.lstrip("-")
    return s


def iedb_stable_allele(xls_token: str) -> str:
    """HLA-A01:01 -> HLA-A*01:01 (matches plot_figure6 ``normalize_hla_netmhc_to_iedb``)."""
    s = xls_token.strip().upper()
    m = re.fullmatch(r"HLA-([ABCEFG])(\d{2}:\d{2})", s)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return s


def read_fasta_peptides(path: Path) -> list[tuple[int, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    out: list[tuple[int, str]] = []
    tag: str | None = None
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            tag = ln[1:].strip()
            continue
        if tag is None:
            continue
        if not tag.isdigit():
            raise SystemExit(f"Unexpected FASTA header (want digits): >{tag}")
        out.append((int(tag), ln.strip().upper()))
        tag = None
    if not out:
        raise SystemExit(f"No peptides parsed from {path}")
    return out


def ic50_nm_to_ba_score(ic50_nm: float) -> float:
    if not math.isfinite(ic50_nm) or ic50_nm <= 0:
        return float("nan")
    return float(max(1e-12, min(1.0 - 1e-12, 1.0 - math.log10(ic50_nm) / math.log10(50000.0))))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--web-tsv", type=Path, required=True, help="NetMHCpan web peptide_table TSV export.")
    ap.add_argument("--out-dir", type=Path, default=None, help="Directory for XLS + companion CSV (default: web-tsv parent).")
    ap.add_argument(
        "--fasta",
        type=Path,
        default=ROOT / "data/netmhc/ttn_as1_108065_ninemers.fasta",
        help="Sliding 9-mer FASTA (defines row order and 0-based IDs for wide XLS).",
    )
    ap.add_argument(
        "--alleles",
        type=Path,
        default=ROOT / "data/netmhc/hla_european27_class1.txt",
        help="Allele order for wide XLS columns (HLA-A01:01 style, one per line).",
    )
    ap.add_argument(
        "--input-seq-id",
        type=str,
        default="108065",
        help="Token for stable_key prefix (must match --iedb-parent-input-seq-id when plotting).",
    )
    args = ap.parse_args()

    out_dir = args.out_dir if args.out_dir is not None else args.web_tsv.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    peptides = read_fasta_peptides(args.fasta)
    alleles_xls = [xls_allele_token(ln) for ln in args.alleles.read_text(encoding="utf-8").splitlines() if ln.strip()]

    df = pd.read_csv(args.web_tsv, sep="\t", low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]
    need = {"peptide", "allele", "peptide index", "netmhcpan_el percentile", "netmhcpan_ba ic50", "netmhcpan_ba percentile"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"Web TSV missing columns {sorted(miss)}; have: {list(df.columns)}")

    imm_col = "immunogenicity score" if "immunogenicity score" in df.columns else None
    # Match common strict web filters: ``processing score`` (not ``processing total score``).
    if "processing score" in df.columns:
        proc_col = "processing score"
    elif "processing total score" in df.columns:
        proc_col = "processing total score"
    else:
        proc_col = None
    if imm_col is None or proc_col is None:
        raise SystemExit(
            "Expected immunogenicity / processing columns not found "
            f"(want 'immunogenicity score' and 'processing score' or 'processing total score'); "
            f"columns={list(df.columns)}"
        )

    el_score_col = "netmhcpan_el score" if "netmhcpan_el score" in df.columns else None

    def norm_allele_web(a: object) -> str:
        return xls_allele_token(str(a))

    df["_al_x"] = df["allele"].map(norm_allele_web)
    df["peptide"] = df["peptide"].astype(str).str.strip().str.upper()
    df["peptide index"] = pd.to_numeric(df["peptide index"], errors="coerce")

    web_alleles = sorted(df["_al_x"].unique())
    if set(web_alleles) != set(alleles_xls):
        raise SystemExit(
            "Allele set in web TSV does not match --alleles file.\n"
            f"  web only: {sorted(set(web_alleles) - set(alleles_xls))}\n"
            f"  panel only: {sorted(set(alleles_xls) - set(web_alleles))}"
        )

    # Lookup: (0-based offset, peptide, allele_xls) -> metrics (web uses 1-based start in "peptide index")
    idx_map: dict[tuple[int, str, str], dict[str, float]] = {}
    for _, r in df.iterrows():
        s1 = int(r["peptide index"])
        if not (1 <= s1 <= 10_000):
            continue
        pos0 = s1 - 1
        pep = str(r["peptide"])
        al = str(r["_al_x"])
        idx_map[(pos0, pep, al)] = {
            "el_pct": float(r["netmhcpan_el percentile"]),
            "ba_pct": float(r["netmhcpan_ba percentile"]),
            "ic50": float(r["netmhcpan_ba ic50"]),
            "imm": float(r[imm_col]),
            "proc": float(r[proc_col]),
            "el_score": float(r[el_score_col]) if el_score_col and pd.notna(r.get(el_score_col)) else float("nan"),
        }

    merged: dict[tuple[str, str], dict[str, float]] = {}
    for pos0, pep in peptides:
        for al in alleles_xls:
            hit = idx_map.get((pos0, pep, al))
            if hit is None:
                continue
            ic50 = hit["ic50"]
            merged[(pep, al)] = {
                "ic50": ic50,
                "ba_rank": hit["ba_pct"],
                "el_rank": hit["el_pct"],
                "el_score": hit["el_score"],
                "_imm": hit["imm"],
                "_proc": hit["proc"],
            }

    xls_path = out_dir / "netmhcpan_ttn_as1_web_netmhc41.xls"
    a_arg = ",".join(alleles_xls)
    line0 = f"# NetMHCpan 4.1 web peptide_table → synthetic -a {a_arg} -l 9 -xls 1 (for Fig 6)"
    line1 = ""
    n_per = 10
    cols = ["Pos", "Peptide", "ID"] + [
        x
        for _ in alleles_xls
        for x in (
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
        )
    ]
    line2 = "\t".join(cols)
    lines = [line0, line1, line2]

    csv_rows: list[dict[str, object]] = []
    input_seq_id = str(args.input_seq_id).strip()

    for pid, pep in peptides:
        pos0 = int(pid)
        s1 = pos0 + 1
        e1 = pos0 + len(pep)
        row: list[str] = [str(pid + 1), pep, str(pos0)]
        for al in alleles_xls:
            m = merged.get((pep, al))
            if m is None:
                row.extend(["", "", "", "", "", "", "0.0", "50.0", "0.0", "50.0"])
                continue
            ic50 = float(m["ic50"])
            ba_score = ic50_nm_to_ba_score(ic50)
            ba_rank = float(m["ba_rank"])
            el_score = float(m["el_score"])
            el_rank = float(m["el_rank"])
            core = icore = pep
            row.extend(
                [
                    core,
                    icore,
                    f"{el_score:.6f}" if math.isfinite(el_score) else "",
                    f"{el_rank:.6f}" if math.isfinite(el_rank) else "",
                    f"{ba_score:.6f}" if math.isfinite(ba_score) else "",
                    f"{ba_rank:.6f}" if math.isfinite(ba_rank) else "",
                    "0.0",
                    "50.0",
                    "0.0",
                    "50.0",
                ]
            )
            al_i = iedb_stable_allele(al)
            stable = f"{input_seq_id}:{s1}:{e1}:{pep}:{al_i}"
            csv_rows.append(
                {
                    "stable_key": stable,
                    "score": m["_imm"],
                    "processing_score": m["_proc"],
                    "netmhcpan_ba_ic50": ic50,
                }
            )
        lines.append("\t".join(row))

    xls_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    csv_path = out_dir / "iedb_companion_for_fig6.csv"
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)

    meta = out_dir / "FIG6_WEB_CONVERSION.txt"
    meta.write_text(
        "\n".join(
            [
                f"source_tsv={args.web_tsv}",
                f"input_seq_id={input_seq_id}",
                f"wide_xls={xls_path}",
                f"iedb_csv={csv_path}",
                "Figure outputs: figures/manuscript_netmhc/fig6/ttn_web_netmhc41/sb_strict_filtered/ (strict web SB), "
                "sb_full/, sb_ic50_lt150nm/",
                f"Run notes: {out_dir / 'figure6'}",
                "Plot (strict web SB: imm>0.1, proc>1.5, EL<=2%, no IC50): plot_figure6_ttn_as1_allele_coverage.py --gating iedb_sb --sb-mode no_ic50 --imm-min 0.1 --proc-min 1.5 --el-rank-max 2 --el-rank-lte \\",
                f"  --netmhc-xls {xls_path.as_posix()} \\",
                f"  --iedb-csv {csv_path.as_posix()} \\",
                f"  --iedb-parent-input-seq-id {input_seq_id} \\",
                "  --split-panels --also-write-unique -o figures/manuscript_netmhc/fig6/ttn_web_netmhc41/sb_strict_filtered/fig6_split.png",
                "Plot (full SB + IC50): .../sb_full/fig6_split.png (see script docstring).",
                "Plot (IC50<150 only): .../sb_ic50_lt150nm/fig6_split.png",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {xls_path}", flush=True)
    print(f"Wrote {csv_path} ({len(csv_rows)} rows)", flush=True)
    print(f"Wrote {meta}", flush=True)


if __name__ == "__main__":
    main()
