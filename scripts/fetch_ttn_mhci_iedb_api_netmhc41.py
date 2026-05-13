#!/usr/bin/env python3
"""
Fetch **NetMHCpan-4.1** class I predictions for TTN-AS1 sliding 9-mers via the **IEDB Tools API**
(REST POST), combining **binding** (`netmhcpan_ba-4.1`) and **eluted ligand** (`netmhcpan_el-4.1`).

This does **not** use a local NetMHCpan binary. It calls::

    https://tools-cluster-interface.iedb.org/tools_api/mhci/

and writes a **synthetic wide XLS** compatible with ``plot_figure6_ttn_as1_allele_coverage.py``
(same 10-column-per-allele block: core, icore, EL_score, EL_rank, BA_score, BA_rank, …).

**BA_score** is reconstructed from IEDB IC50 (nM) using the usual NetMHCpan mapping
``IC50 = 50000 ** (1 - BA_score)``. **BA_rank** and **EL_rank** use IEDB ``percentile_rank``;
**EL_score** uses IEDB ``score`` from the EL method response.

Example::

    python scripts/fetch_ttn_mhci_iedb_api_netmhc41.py --out-dir data/netmhc/predictions/ttn_as1_smpep108065_iedb_api_netmhc41

Then plot Fig 6::

    python plot_figure6_ttn_as1_allele_coverage.py --netmhc-xls data/netmhc/predictions/ttn_as1_smpep108065_iedb_api_netmhc41/netmhcpan_ttn_as1_iedb_api_netmhc41.xls --split-panels --also-write-unique -o figures/manuscript_netmhc/fig6/ttn_iedb_api_netmhc41/fig6_split.png

On Windows, if TLS revocation checks fail, add ``--insecure`` (same idea as ``curl -k``).

Large FASTA payloads can trigger **HTTP 500** from the public cluster; the script defaults to
**batched** POSTs (``--batch-size``) and **retries** on transient 5xx.
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

IEDB_MHCI_URL_DEFAULT = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"


def ic50_nm_to_ba_score(ic50_nm: float) -> float:
    """Inverse of NetMHCpan BA_score → IC50 mapping (same as repo-wide scripts)."""
    if not math.isfinite(ic50_nm) or ic50_nm <= 0:
        return float("nan")
    return float(max(1e-12, min(1.0 - 1e-12, 1.0 - math.log10(ic50_nm) / math.log10(50000.0))))


def xls_allele_token(line: str) -> str:
    """HLA-A01:01 style (no asterisk), matching ``hla_european27_class1.txt``."""
    s = line.strip().upper().replace("*", "")
    if not s.startswith("HLA-"):
        s = "HLA-" + s.lstrip("-")
    return s


def iedb_allele_token(xls_token: str) -> str:
    """HLA-A*01:01 style for IEDB POST ``allele`` parameter."""
    s = xls_token.strip().upper()
    m = re.fullmatch(r"HLA-([ABCEFG])(\d{2}:\d{2})", s)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return s


def read_fasta_peptides(path: Path) -> list[tuple[int, str]]:
    """Return (numeric_id, peptide) in file order; id from header ``>000``."""
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


def read_allele_list(path: Path) -> list[str]:
    rows = [xls_allele_token(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not rows:
        raise SystemExit(f"No alleles in {path}")
    return rows


def post_form(
    url: str,
    fields: dict[str, str],
    *,
    insecure: bool,
    max_retries: int = 3,
    retry_backoff_s: float = 2.0,
) -> str:
    data = urllib.parse.urlencode(fields).encode("utf-8")
    ctx = ssl.create_default_context()
    if insecure:
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    retry_codes = {500, 502, 503, 504}
    backoff = float(retry_backoff_s)
    last_code = 0
    last_body = ""
    for attempt in range(max_retries + 1):
        req = urllib.request.Request(
            url,
            data=data,
            method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        try:
            with urllib.request.urlopen(req, timeout=600, context=ctx if insecure else None) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            last_body = e.read().decode("utf-8", errors="replace") if e.fp else ""
            last_code = e.code
            if e.code in retry_codes and attempt < max_retries:
                time.sleep(backoff)
                backoff = min(backoff * 2.0, 60.0)
                continue
            preview = last_body[:2000] if last_body else "(empty body)"
            raise SystemExit(
                f"IEDB HTTP {last_code} (attempt {attempt + 1}/{max_retries + 1}); "
                f"body_len={len(last_body)} preview:\n{preview}"
            ) from e
    raise SystemExit(f"IEDB HTTP {last_code}: unreachable retry loop exit")


def parse_iedb_tsv(tsv_text: str) -> list[dict[str, str]]:
    t = tsv_text.strip()
    if not t:
        return []
    f = io.StringIO(t)
    r = csv.DictReader(f, delimiter="\t")
    return [dict(row) for row in r]


def write_dict_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=keys, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for row in rows:
        w.writerow({k: str(row.get(k, "")) for k in keys})
    path.write_text(buf.getvalue(), encoding="utf-8")


def fetch_method_for_allele(
    *,
    url: str,
    method: str,
    allele_iedb: str,
    fasta: str,
    insecure: bool,
    max_retries: int,
    retry_backoff_s: float,
) -> list[dict[str, str]]:
    fields = {
        "method": method,
        "allele": allele_iedb,
        "length": "9",
        "sequence_text": fasta,
    }
    raw = post_form(
        url,
        fields,
        insecure=insecure,
        max_retries=max_retries,
        retry_backoff_s=retry_backoff_s,
    )
    return parse_iedb_tsv(raw)


def peptide_batches(peptides: list[tuple[int, str]], batch_size: int) -> list[list[tuple[int, str]]]:
    if batch_size <= 0:
        return [peptides]
    return [peptides[i : i + batch_size] for i in range(0, len(peptides), batch_size)]


def fetch_method_for_allele_batched(
    *,
    url: str,
    method: str,
    allele_iedb: str,
    peptides: list[tuple[int, str]],
    batch_size: int,
    insecure: bool,
    max_retries: int,
    retry_backoff_s: float,
    sleep_between_batches_s: float,
) -> list[dict[str, str]]:
    chunks = peptide_batches(peptides, batch_size)
    all_rows: list[dict[str, str]] = []
    for ci, chunk in enumerate(chunks):
        if len(chunks) > 1:
            print(f"  batch {ci + 1}/{len(chunks)} ({len(chunk)} peptides)", flush=True)
        fasta = build_fasta(chunk)
        rows = fetch_method_for_allele(
            url=url,
            method=method,
            allele_iedb=allele_iedb,
            fasta=fasta,
            insecure=insecure,
            max_retries=max_retries,
            retry_backoff_s=retry_backoff_s,
        )
        all_rows.extend(rows)
        if ci < len(chunks) - 1:
            time.sleep(float(sleep_between_batches_s))
    return all_rows


def build_fasta(peptides: list[tuple[int, str]]) -> str:
    parts: list[str] = []
    for pid, pep in peptides:
        parts.append(f">{pid:03d}")
        parts.append(pep)
    return "\n".join(parts) + "\n"


def merge_predictions(
    _alleles_xls: list[str],
    peptides: list[tuple[int, str]],
    ba_by_allele: dict[str, list[dict[str, str]]],
    el_by_allele: dict[str, list[dict[str, str]]],
) -> dict[tuple[str, str], dict[str, float]]:
    """
    Key (peptide, allele_xls). Values: ic50, ba_rank_pct, el_score, el_rank_pct.
    """
    out: dict[tuple[str, str], dict[str, float]] = {}
    pset = {p for _, p in peptides}

    def ingest(rows: list[dict[str, str]], *, kind: str) -> None:
        for row in rows:
            pep = str(row.get("peptide", "")).strip().upper()
            if pep not in pset:
                continue
            al = str(row.get("allele", "")).strip().upper()
            al_x = xls_allele_token(al.replace("*", ""))
            key = (pep, al_x)
            if key not in out:
                out[key] = {}
            if kind == "ba":
                ic50 = float(row.get("ic50", "nan"))
                out[key]["ic50"] = ic50
                out[key]["ba_rank"] = float(row.get("percentile_rank", "nan"))
            else:
                out[key]["el_score"] = float(row.get("score", "nan"))
                out[key]["el_rank"] = float(row.get("percentile_rank", "nan"))

    for rows in ba_by_allele.values():
        ingest(rows, kind="ba")
    for rows in el_by_allele.values():
        ingest(rows, kind="el")
    return out


def write_synthetic_wide_xls(
    path: Path,
    *,
    alleles_xls: list[str],
    peptides: list[tuple[int, str]],
    merged: dict[tuple[str, str], dict[str, float]],
) -> None:
    a_arg = ",".join(alleles_xls)
    line0 = f"# IEDB tools_api/mhci netmhcpan_ba-4.1 + netmhcpan_el-4.1 (merged) -a {a_arg} -l 9 -xls 1 (synthetic header for Fig 6)"
    # Line 1 is not parsed by ``plot_figure6`` (only line 0 supplies alleles).
    line1 = ""
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
    for pid, pep in peptides:
        row: list[str] = [str(pid + 1), pep, f"{pid:03d}"]
        for al in alleles_xls:
            m = merged.get((pep, al), {})
            ic50 = float(m.get("ic50", float("nan")))
            ba_score = ic50_nm_to_ba_score(ic50)
            ba_rank = float(m.get("ba_rank", float("nan")))
            el_score = float(m.get("el_score", float("nan")))
            el_rank = float(m.get("el_rank", float("nan")))
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
        lines.append("\t".join(row))

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "data/netmhc/predictions/ttn_as1_smpep108065_iedb_api_netmhc41",
        help="Directory for TSV logs + synthetic wide XLS.",
    )
    ap.add_argument(
        "--fasta",
        type=Path,
        default=ROOT / "data/netmhc/ttn_as1_108065_ninemers.fasta",
        help="Sliding 9-mer FASTA (same as local NetMHCpan runs).",
    )
    ap.add_argument(
        "--alleles",
        type=Path,
        default=ROOT / "data/netmhc/hla_european27_class1.txt",
        help="One allele per line (HLA-A01:01 style).",
    )
    ap.add_argument(
        "--url",
        type=str,
        default=IEDB_MHCI_URL_DEFAULT,
        help="IEDB MHC-I binding API endpoint.",
    )
    ap.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Seconds between API calls and between peptide batches (be polite to the public server).",
    )
    ap.add_argument(
        "--batch-size",
        type=int,
        default=20,
        help="Peptides per POST (smaller avoids IEDB HTTP 500 on large FASTA). Use 0 for one POST with all peptides.",
    )
    ap.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Retries per POST on HTTP 500/502/503/504 (exponential backoff).",
    )
    ap.add_argument(
        "--retry-backoff",
        type=float,
        default=2.0,
        help="Initial seconds before first retry on 5xx.",
    )
    ap.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification (use if Windows revocation check fails).",
    )
    args = ap.parse_args()

    peptides = read_fasta_peptides(args.fasta)
    alleles_xls = read_allele_list(args.alleles)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    ba_by_allele: dict[str, list[dict[str, str]]] = {}
    el_by_allele: dict[str, list[dict[str, str]]] = {}

    for al in alleles_xls:
        ai = iedb_allele_token(al)
        print(f"[IEDB] BA netmhcpan_ba-4.1 {ai} …", flush=True)
        rows_ba = fetch_method_for_allele_batched(
            url=args.url,
            method="netmhcpan_ba-4.1",
            allele_iedb=ai,
            peptides=peptides,
            batch_size=int(args.batch_size),
            insecure=bool(args.insecure),
            max_retries=int(args.max_retries),
            retry_backoff_s=float(args.retry_backoff),
            sleep_between_batches_s=float(args.sleep),
        )
        write_dict_tsv(out_dir / f"raw_ba_{ai.replace('*', '_')}.tsv", rows_ba)
        ba_by_allele[ai] = rows_ba
        time.sleep(float(args.sleep))

        print(f"[IEDB] EL netmhcpan_el-4.1 {ai} …", flush=True)
        rows_el = fetch_method_for_allele_batched(
            url=args.url,
            method="netmhcpan_el-4.1",
            allele_iedb=ai,
            peptides=peptides,
            batch_size=int(args.batch_size),
            insecure=bool(args.insecure),
            max_retries=int(args.max_retries),
            retry_backoff_s=float(args.retry_backoff),
            sleep_between_batches_s=float(args.sleep),
        )
        write_dict_tsv(out_dir / f"raw_el_{ai.replace('*', '_')}.tsv", rows_el)
        el_by_allele[ai] = rows_el
        time.sleep(float(args.sleep))

    merged = merge_predictions(alleles_xls, peptides, ba_by_allele, el_by_allele)
    out_xls = out_dir / "netmhcpan_ttn_as1_iedb_api_netmhc41.xls"
    write_synthetic_wide_xls(out_xls, alleles_xls=alleles_xls, peptides=peptides, merged=merged)

    meta = out_dir / "RUN_INFO.txt"
    meta.write_text(
        "\n".join(
            [
                f"url={args.url}",
                "methods=netmhcpan_ba-4.1, netmhcpan_el-4.1",
                f"batch_size={args.batch_size}",
                f"max_retries={args.max_retries}",
                f"n_peptides={len(peptides)}",
                f"n_alleles={len(alleles_xls)}",
                f"out_xls={out_xls}",
                "Note: BA_score reconstructed from IEDB IC50 (nM); ranks from IEDB percentile_rank.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {out_xls}", flush=True)


if __name__ == "__main__":
    main()
