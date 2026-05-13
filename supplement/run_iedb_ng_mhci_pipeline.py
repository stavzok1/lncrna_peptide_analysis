#!/usr/bin/env python3
"""
Submit IEDB Next-Generation Class I (mhci) jobs and poll for results.

Uses only the standard library (urllib). Docs:
https://nextgen-tools.iedb.org/docs/api/endpoints/api_references.html
Tool parameters: GET https://api-nextgen-tools.iedb.org/api/v1/mhci

Note: "basic_processing" (proteasome + TAP + MHC in pathway) typically emits
tap_score / proteasome_score / processing_score on *protein* inputs with a
sliding peptide_length_range. For isolated short peptides (peptide_length_range
null), the API may return binding columns only for that processing head.

IEDB rejects FASTA records **longer than 14 aa** when peptide_length_range is omitted
(peptide-list mode). This script then auto-uses sliding **9-mers** (9-9) so long
parents are processed instead of dropped.

Examples:
  python scripts/run_iedb_ng_mhci_pipeline.py \\
    --fasta epitopes.faa --alleles "HLA-A*02:01" \\
    --out-json out/iedb_mhci.json --out-csv out/iedb_mhci_peptide_table.csv

  python scripts/run_iedb_ng_mhci_pipeline.py \\
    --fasta parent_proteins.faa --alleles-file data/netmhc/hla_european27_class1.txt \\
    --peptide-length-range 9-9 \\
    --out-json out/iedb_processing.json --out-csv out/iedb_processing.csv

  # Whole peptide library: chunk + stable IDs (input_seq_id, stable_key) + EL + TAP + immunogenicity
  python scripts/run_iedb_ng_mhci_pipeline.py --basic-processing --immunogenicity \\
    --fasta your_peptides.faa --alleles-file data/netmhc/hla_european27_class1.txt \\
    --chunk-size 40 --out-json out/iedb_batch.json --out-csv out/iedb_merged.csv

Parameter notes (see epilog from --help):
  - EL: default --binding / --binding-method netmhcpan_el (omit --no-binding).
  - basic_processing uses BA internally (--basic-processing-binding); independent of EL head.
  - Immunogenicity: use --mask-choice by_allele for multi-allele runs unless you need custom positions.
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
import csv
import json
import re
import time
import urllib.error
import urllib.request
from typing import Any, Iterable

API_PIPELINE = "https://api-nextgen-tools.iedb.org/api/v1/pipeline"
API_MHCI = "https://api-nextgen-tools.iedb.org/api/v1/mhci"
TERMINAL_OK = {"done", "complete"}
TERMINAL_FAIL = {"failed", "error"}


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _load_fasta_or_raw(path: Path) -> str:
    """Return file contents as-is (IEDB expects FASTA or line-oriented peptide blocks)."""
    return _read_text(path).strip() + "\n"


def _parse_fasta_records(text: str) -> list[tuple[str, str]]:
    """Return list of (record_id, sequence) from FASTA text. Skips empty sequences."""
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


def _fasta_from_records(chunk: list[tuple[str, str]]) -> str:
    lines: list[str] = []
    for rid, seq in chunk:
        lines.append(f">{rid}\n{seq}")
    return "\n".join(lines) + "\n"


def _normalize_iedb_hla_allele(name: str) -> str:
    """
    NetMHCpan allele files often use HLA-A02:01 (no star). IEDB Next-Gen expects
    HLA-A*02:01 for classical class I. Non-matching tokens (e.g. H2-Kb) are unchanged.
    """
    s = name.strip()
    if not s:
        return s
    up = s.upper()
    if re.fullmatch(r"HLA-([ABCEFG])\*(\d{2}:\d{2})", up):
        return up
    m = re.fullmatch(r"HLA-([ABCEFG])(\d{2}:\d{2})", up)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    m = re.fullmatch(r"HLA-([ABCEFG])(\d{2})(\d{2})", up)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    return s


def _normalize_allele_list(raw: str, *, enabled: bool) -> str:
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if not parts:
        return ""
    if enabled:
        parts = [_normalize_iedb_hla_allele(p) for p in parts]
    return ",".join(parts)


def _load_alleles(arg: str | None, path: Path | None, *, normalize: bool) -> str:
    if path:
        lines = []
        for line in _read_text(path).splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            lines.append(s)
        if not lines:
            raise SystemExit(f"No alleles found in {path}")
        joined = ",".join(lines)
        return _normalize_allele_list(joined, enabled=normalize)
    if not arg:
        raise SystemExit("Provide --alleles and/or --alleles-file")
    return _normalize_allele_list(arg.strip(), enabled=normalize)


def _parse_length_range(s: str | None) -> list[int] | None:
    if s is None or s.lower() in ("null", "none", ""):
        return None
    s = s.replace(" ", "")
    if "-" in s:
        a, b = s.split("-", 1)
        return [int(a), int(b)]
    if "," in s:
        a, b = s.split(",", 1)
        return [int(a), int(b)]
    n = int(s)
    return [n, n]


def _http_json(
    url: str,
    method: str = "GET",
    payload: dict[str, Any] | None = None,
    timeout: int = 120,
) -> Any:
    data = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        raw = json.dumps(payload).encode("utf-8")
        data = raw
        headers["Content-Type"] = "application/json; charset=utf-8"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        err_body = e.read().decode("utf-8", errors="replace")
        raise SystemExit(f"HTTP {e.code} {url}\n{err_body}") from e
    except urllib.error.URLError as e:
        raise SystemExit(f"Request failed: {url}\n{e}") from e
    if not body.strip():
        return None
    return json.loads(body)


def _poll_results(
    results_uri: str,
    *,
    poll_interval: float,
    timeout: int,
    progress_prefix: str,
    verbose: bool,
    poll_log_every: int,
) -> dict[str, Any]:
    deadline = time.monotonic() + timeout
    last: dict[str, Any] | None = None
    poll_n = 0
    t0 = time.monotonic()
    log_every = max(1, poll_log_every)
    while time.monotonic() < deadline:
        poll_n += 1
        last = _http_json(results_uri, "GET")  # type: ignore[assignment]
        status = (last or {}).get("status")
        if verbose and (poll_n == 1 or poll_n % log_every == 0):
            elapsed = int(time.monotonic() - t0)
            print(
                f"{progress_prefix}[poll {poll_n}] status={status!r} elapsed={elapsed}s",
                file=sys.stderr,
                flush=True,
            )
        if status in TERMINAL_OK:
            if verbose and poll_n > 1 and poll_n % log_every != 0:
                elapsed = int(time.monotonic() - t0)
                print(
                    f"{progress_prefix}[poll {poll_n}] status={status!r} elapsed={elapsed}s (final)",
                    file=sys.stderr,
                    flush=True,
                )
            return last  # type: ignore[return-value]
        if status in TERMINAL_FAIL:
            raise SystemExit(f"Pipeline status={status!r}\n{json.dumps(last, indent=2)}")
        time.sleep(poll_interval)
    raise SystemExit(
        f"{progress_prefix}Timed out after {timeout}s. Last status: {(last or {}).get('status')!r}"
    )


def _submit_pipeline(pipeline: dict[str, Any], post_timeout: int) -> str:
    submitted = _http_json(API_PIPELINE, "POST", pipeline, timeout=post_timeout)
    results_uri = submitted.get("results_uri")
    if not results_uri:
        hint = ""
        errs = submitted.get("errors") or []
        warns = submitted.get("warnings") or []
        blob = json.dumps(submitted, indent=2)
        if any("maximum peptide length (14)" in str(w) for w in warns):
            hint = (
                "\n\nIEDB: with peptide_length_range omitted, each FASTA sequence must be ≤14 aa. "
                "Re-run with --peptide-length-range 9-9 (sliding 9-mers on parents), or omit only "
                "short peptides. This script auto-switches to 9-9 when the longest input exceeds 14 aa."
            )
        if errs or warns:
            raise SystemExit(f"IEDB pipeline POST rejected or produced no job.{hint}\n{blob}")
        raise SystemExit(f"Unexpected POST response (no results_uri):\n{blob}")
    return results_uri


def _build_predictors(args: argparse.Namespace) -> list[dict[str, Any]]:
    preds: list[dict[str, Any]] = []
    if args.binding:
        preds.append({"type": "binding", "method": args.binding_method})
    if args.basic_processing:
        preds.append(
            {
                "type": "processing",
                "method": "basic_processing",
                "mhc_binding_method": args.basic_processing_binding,
                "proteasome": args.proteasome,
                "tap_precursor": int(args.tap_precursor),
                "tap_alpha": float(args.tap_alpha),
            }
        )
    if args.netchop:
        preds.append(
            {
                "type": "processing",
                "method": "netchop",
                "network_method": args.netchop_network,
                "threshold": float(args.netchop_threshold),
            }
        )
    if args.immunogenicity:
        p: dict[str, Any] = {
            "type": "immunogenicity",
            "method": "immunogenicity",
            "mask_choice": args.mask_choice,
        }
        if args.mask_choice == "custom" and args.position_to_mask:
            p["position_to_mask"] = args.position_to_mask
        preds.append(p)
    if not preds:
        raise SystemExit("No predictors selected. Enable binding and/or processing and/or immunogenicity.")
    return preds


def _first_peptide_table(results_obj: dict[str, Any]) -> dict[str, Any] | None:
    data = results_obj.get("data") or {}
    for block in data.get("results") or []:
        if block.get("type") == "peptide_table":
            return block
    return None


def _peptide_table_to_rows(table: dict[str, Any]) -> tuple[list[str], list[list[Any]]]:
    cols = [c["name"] for c in table.get("table_columns") or []]
    rows = table.get("table_data") or []
    return cols, rows


def _write_csv(path: Path, columns: list[str], rows: Iterable[list[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for row in rows:
            w.writerow(row)


def _append_csv_rows(path: Path, columns: list[str], rows: Iterable[list[Any]], *, write_header: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        if write_header:
            w.writerow(columns)
        for row in rows:
            w.writerow(row)


def _row_indices_for_columns(have: list[str], need: list[str]) -> list[int] | None:
    idx = []
    for name in need:
        try:
            idx.append(have.index(name))
        except ValueError:
            return None
    return idx


def _reorder_rows(rows: list[list[Any]], col_names: list[str], master_cols: list[str]) -> list[list[Any]]:
    mapping = _row_indices_for_columns(col_names, master_cols)
    if mapping is None:
        raise SystemExit(
            f"Column set mismatch between chunks. Expected columns {master_cols!r} "
            f"but got {col_names!r}"
        )
    return [[row[i] for i in mapping] for row in rows]


STABLE_ID_COLS = ["iedb_chunk_index", "input_seq_id", "stable_key"]


def _enrich_rows_stable_ids(
    cols: list[str],
    rows: list[list[Any]],
    *,
    chunk_index: int,
    ordered_seq_ids: list[str],
) -> tuple[list[str], list[list[Any]]]:
    """
    IEDB ``sequence_number`` restarts every job chunk. Map it to FASTA header order in that chunk
    and append stable join keys. ``input_seq_id`` = first FASTA header token; ``stable_key``
    encodes seq id + window + peptide + allele.
    """
    try:
        i_sn = cols.index("sequence_number")
        i_start = cols.index("start")
        i_end = cols.index("end")
        i_pep = cols.index("peptide")
        i_allele = cols.index("allele")
    except ValueError as e:
        raise SystemExit(f"peptide_table missing expected column: {e}") from e
    ext_cols = cols + STABLE_ID_COLS
    out: list[list[Any]] = []
    for row in rows:
        sn = int(float(row[i_sn]))
        if sn < 1 or sn > len(ordered_seq_ids):
            raise SystemExit(
                f"Chunk {chunk_index}: sequence_number={sn} but chunk has {len(ordered_seq_ids)} "
                f"FASTA records (check FASTA order vs IEDB)."
            )
        sid = ordered_seq_ids[sn - 1]
        start = row[i_start]
        end = row[i_end]
        pep = row[i_pep]
        alle = row[i_allele]
        stable = f"{sid}:{start}:{end}:{pep}:{alle}"
        out.append(list(row) + [chunk_index, sid, stable])
    return ext_cols, out


def main() -> None:
    ap = argparse.ArgumentParser(
        description="IEDB Next-Gen MHCI pipeline (CLI via urllib).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Predictor / parameter hints
---------------------------
  EL (presentation): IEDB ``binding`` + ``netmhcpan_el``. Default ON; use ``--no-binding`` only if
  you truly want processing/immunogenicity without an extra EL head (you already have local NetMHC).

  Processing + TAP: ``basic_processing`` always uses an MHC IC50 method internally
  (``--basic-processing-binding``, default netmhcpan_ba). That is separate from the EL head.

  TAP knobs (``--tap-precursor``, ``--tap-alpha``): IEDB defaults (1, 0.2). Tune only if you are
  reproducing a specific published pathway setting.

  Proteasome: ``--proteasome immuno`` (default) vs constitutive - pick what matches your biology.

  Immunogenicity: for many alleles, ``--mask-choice by_allele`` is often more appropriate than
  ``custom`` + fixed ``--position-to-mask``. Custom 2,5,9 is the IEDB doc example for 9-mers.

  Chunked CSV: rows get ``iedb_chunk_index``, ``input_seq_id`` (FASTA id), ``stable_key`` so
  ``sequence_number`` need not be unique across the merged file.
""",
    )
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--fasta", type=Path, help="Path to FASTA (proteins or one peptide per record).")
    src.add_argument("--sequence-text", help="Raw input_sequence_text (for quick tests).")

    ap.add_argument("--alleles", help="Comma-separated allele names, e.g. HLA-A*02:01,HLA-B*07:02")
    ap.add_argument(
        "--alleles-file",
        type=Path,
        help="One allele per line. NetMHC-style HLA-A02:01 is rewritten to HLA-A*02:01 unless --no-normalize-alleles.",
    )
    ap.add_argument(
        "--no-normalize-alleles",
        action="store_true",
        help="Send allele strings exactly as given (skip HLA-A02:01 -> HLA-A*02:01 conversion).",
    )
    ap.add_argument(
        "--peptide-length-range",
        metavar="RANGE",
        help='Sliding window on each FASTA entry, e.g. 9-9. Omit (null) = each record is one peptide '
        "(IEDB max length **14 aa**; longer parents are dropped). Default null; see auto-switch in --help.",
    )
    ap.add_argument("--out-json", type=Path, help="Write full results JSON here (required unless --list-predictors).")
    ap.add_argument("--out-csv", type=Path, help="Write first peptide_table to CSV (optional).")
    ap.add_argument("--email", help="Optional email for IEDB completion notices.")
    ap.add_argument("--poll-interval", type=float, default=3.0, help="Seconds between GET results_uri polls.")
    ap.add_argument("--timeout", type=int, default=900, help="Max total seconds waiting for results (per chunk).")
    ap.add_argument(
        "--http-post-timeout",
        type=int,
        default=300,
        help="Timeout seconds for the pipeline POST (large FASTA payloads).",
    )
    ap.add_argument(
        "--chunk-size",
        type=int,
        default=0,
        metavar="N",
        help="If >0 and input is --fasta, split into chunks of N sequences (one IEDB job per chunk). "
        "Use for whole libraries; merge rows into --out-csv. Default 0 = single job for entire file.",
    )
    ap.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable stderr progress lines (poll status and chunk counters).",
    )
    ap.add_argument(
        "--poll-log-every",
        type=int,
        default=1,
        metavar="N",
        help="Print a poll line every N polls (default 1). Use e.g. 5 to reduce stderr volume.",
    )
    ap.add_argument("--dry-run", action="store_true", help="Print pipeline JSON and exit (no HTTP).")

    ap.add_argument("--binding", dest="binding", action="store_true", default=True)
    ap.add_argument("--no-binding", dest="binding", action="store_false")
    ap.add_argument(
        "--binding-method",
        default="netmhcpan_el",
        help="binding predictor short_name (default: netmhcpan_el).",
    )
    ap.add_argument("--basic-processing", action="store_true", help="Add basic_processing (proteasome + TAP pathway).")
    ap.add_argument(
        "--basic-processing-binding",
        default="netmhcpan_ba",
        help="mhc_binding_method for basic_processing (default: netmhcpan_ba).",
    )
    ap.add_argument("--proteasome", choices=("immuno", "constitutive"), default="immuno")
    ap.add_argument("--tap-precursor", type=int, default=1, help="basic_processing tap_precursor (default 1).")
    ap.add_argument("--tap-alpha", type=float, default=0.2, help="basic_processing tap_alpha (default 0.2).")
    ap.add_argument("--netchop", action="store_true", help="Add NetChop processing head.")
    ap.add_argument("--netchop-network", choices=("c_term", "20s"), default="c_term")
    ap.add_argument("--netchop-threshold", type=float, default=0.5)
    ap.add_argument("--immunogenicity", action="store_true", help="Add Class I immunogenicity score.")
    ap.add_argument(
        "--mask-choice",
        choices=("custom", "display", "by_allele"),
        default="custom",
        help="immunogenicity mask_choice (default custom 2,5,9). Try by_allele for multi-allele panels.",
    )
    ap.add_argument(
        "--position-to-mask",
        default="2,5,9",
        help="Only if --mask-choice custom: comma-separated 1-based positions (default 2,5,9).",
    )
    ap.add_argument(
        "--defaults-tap-immuno",
        action="store_true",
        help="Shorthand: --basic-processing --immunogenicity (EL binding stays on unless --no-binding).",
    )
    ap.add_argument(
        "--panel-el-tap-immune",
        action="store_true",
        help="Shorthand: --binding --basic-processing --immunogenicity (EL + pathway + immunogenicity).",
    )
    ap.add_argument(
        "--list-predictors",
        action="store_true",
        help="Print GET /mhci predictor catalog as JSON and exit.",
    )
    args = ap.parse_args()
    if args.list_predictors:
        spec = _http_json(API_MHCI, "GET")
        json.dump(spec, sys.stdout, indent=2)
        print()
        return

    if not args.out_json:
        ap.error("--out-json is required (unless using --list-predictors).")

    if not args.fasta and not args.sequence_text:
        ap.error("Provide --fasta and/or --sequence-text.")

    if args.defaults_tap_immuno:
        args.basic_processing = True
        args.immunogenicity = True

    if args.panel_el_tap_immune:
        args.binding = True
        args.basic_processing = True
        args.immunogenicity = True

    if not (args.binding or args.basic_processing or args.netchop or args.immunogenicity):
        args.binding = True
        args.basic_processing = True
        args.immunogenicity = True

    verbose = not args.no_progress
    alleles = _load_alleles(args.alleles, args.alleles_file, normalize=not args.no_normalize_alleles)
    pep_range = _parse_length_range(args.peptide_length_range)
    predictors = _build_predictors(args)

    def build_pipeline(seq_text: str) -> dict[str, Any]:
        pl: dict[str, Any] = {
            "pipeline_id": "",
            "run_stage_range": [1, 1],
            "stages": [
                {
                    "stage_number": 1,
                    "tool_group": "mhci",
                    "input_sequence_text": seq_text,
                    "input_parameters": {
                        "alleles": alleles,
                        "peptide_length_range": pep_range,
                        "predictors": predictors,
                    },
                }
            ],
        }
        if args.email:
            pl["email"] = args.email
        return pl

    chunk_seq_ids: list[list[str]]

    if args.sequence_text:
        txt = args.sequence_text.strip() + "\n"
        job_texts = [txt]
        tr_inline = _parse_fasta_records(txt)
        if tr_inline:
            chunk_seq_ids = [[rid for rid, _ in tr_inline]]
            n_seq = len(tr_inline)
        else:
            chunk_seq_ids = [["1"]]
            n_seq = 1
    else:
        raw = _read_text(args.fasta)  # type: ignore[union-attr]
        recs = _parse_fasta_records(raw)
        if not recs:
            raise SystemExit(f"No FASTA sequences found in {args.fasta}")
        n_seq = len(recs)
        if args.chunk_size > 0:
            step = max(1, args.chunk_size)
            job_texts = [_fasta_from_records(recs[i : i + step]) for i in range(0, len(recs), step)]
            chunk_seq_ids = [[rid for rid, _ in recs[i : i + step]] for i in range(0, len(recs), step)]
        else:
            job_texts = [raw.strip() + "\n"]
            chunk_seq_ids = [[rid for rid, _ in recs]]

    n_chunks = len(job_texts)

    all_lens = [len(s) for block in job_texts for _, s in _parse_fasta_records(block)]
    if pep_range is None and all_lens and max(all_lens) > 14:
        pep_range = [9, 9]
        print(
            "[iedb] IEDB Class I with peptide_length_range omitted treats each FASTA record as one "
            "peptide and accepts **length ≤14 aa** only (longer records are excluded). "
            f"Longest sequence here is **{max(all_lens)} aa** — using sliding **9-mers** "
            "(same as --peptide-length-range 9-9). Output rows are 9-mers; columns `start`/`end` "
            "refer to positions within each parent record.",
            file=sys.stderr,
            flush=True,
        )
        if min(all_lens) < 9:
            print(
                f"[iedb] WARNING: some parents are shorter than 9 aa (min={min(all_lens)}); "
                "those records may yield no peptides.",
                file=sys.stderr,
                flush=True,
            )

    if args.dry_run:
        json.dump(build_pipeline(job_texts[0]), sys.stdout, indent=2)
        print()
        return

    if args.chunk_size > 0 and args.sequence_text:
        print("Warning: --chunk-size applies to --fasta only; ignoring for --sequence-text.", file=sys.stderr)

    if args.basic_processing and pep_range is None:
        lens = []
        for block in job_texts:
            for line in block.splitlines():
                if line.startswith(">"):
                    continue
                s = line.strip()
                if s:
                    lens.append(len(s))
        max_len = max(lens, default=0)
        if max_len and max_len <= 15:
            print(
                "Warning: basic_processing with peptide_length_range null and short sequences "
                "may omit tap_score/processing_score; use parent protein FASTA and --peptide-length-range 9-9 "
                "for TAP columns.",
                file=sys.stderr,
            )

    n_allele_warn = len([a for a in alleles.split(",") if a.strip()])
    if args.chunk_size == 0 and n_seq * n_allele_warn > 8000:
        print(
            f"Warning: single job with ~{n_seq * n_allele_warn} rows may be slow or hit size limits; "
            "consider --chunk-size 30–80.",
            file=sys.stderr,
            flush=True,
        )

    if verbose:
        n_allele = len([a for a in alleles.split(",") if a.strip()])
        if pep_range == [9, 9] and all_lens:
            approx_rows = sum(max(0, L - 8) for L in all_lens) * n_allele
        else:
            approx_rows = n_seq * n_allele
        mode = f"sliding {pep_range[0]}-{pep_range[1]}" if pep_range else "peptide-list (≤14 aa/record)"
        print(
            f"[iedb] {n_seq} FASTA record(s), {n_chunks} job(s), {n_allele} allele(s), mode={mode}, "
            f"~{approx_rows} peptide×allele rows (estimate); poll every {args.poll_interval}s, "
            f"timeout {args.timeout}s/job",
            file=sys.stderr,
            flush=True,
        )

    all_results: list[dict[str, Any]] = []
    master_cols: list[str] | None = None
    total_rows = 0

    if args.out_csv and args.out_csv.exists():
        args.out_csv.unlink(missing_ok=True)

    for ci, block in enumerate(job_texts, start=1):
        prefix = f"[iedb chunk {ci}/{n_chunks}] " if n_chunks > 1 else "[iedb] "
        ordered_ids = chunk_seq_ids[ci - 1]
        if verbose:
            nrec = block.count(">")
            print(f"{prefix}POST ({nrec} sequence(s)) …", file=sys.stderr, flush=True)
        pipeline = build_pipeline(block)
        results_uri = _submit_pipeline(pipeline, args.http_post_timeout)
        if verbose:
            print(f"{prefix}submitted results_uri={results_uri}", file=sys.stderr, flush=True)
        last = _poll_results(
            results_uri,
            poll_interval=args.poll_interval,
            timeout=args.timeout,
            progress_prefix=prefix if verbose else "",
            verbose=verbose,
            poll_log_every=args.poll_log_every,
        )
        all_results.append(last)
        if args.out_csv:
            pt = _first_peptide_table(last)
            if not pt:
                print(f"{prefix}No peptide_table in results; skipping CSV rows for this chunk.", file=sys.stderr)
            else:
                cols, rows = _peptide_table_to_rows(pt)
                cols_e, rows_e = _enrich_rows_stable_ids(
                    cols, rows, chunk_index=ci, ordered_seq_ids=ordered_ids
                )
                if master_cols is None:
                    master_cols = cols_e
                    _append_csv_rows(args.out_csv, master_cols, rows_e, write_header=True)
                else:
                    rows_e = _reorder_rows(rows_e, cols_e, master_cols)
                    _append_csv_rows(args.out_csv, master_cols, rows_e, write_header=False)
                total_rows += len(rows_e)
                if verbose:
                    print(
                        f"{prefix}done peptide_table rows={len(rows_e)} cumulative={total_rows}",
                        file=sys.stderr,
                        flush=True,
                    )

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    if n_chunks == 1:
        args.out_json.write_text(json.dumps(all_results[0], indent=2), encoding="utf-8")
    else:
        bundle = {
            "format": "iedb_ng_mhci_batch_v1",
            "n_chunks": n_chunks,
            "chunk_size": args.chunk_size,
            "results": all_results,
        }
        args.out_json.write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    print(f"Wrote {args.out_json}", file=sys.stderr)

    if args.out_csv and master_cols is not None:
        print(f"Wrote {args.out_csv} ({total_rows} rows)", file=sys.stderr)
    elif args.out_csv:
        print("No peptide_table in any chunk; no CSV written.", file=sys.stderr)


if __name__ == "__main__":
    main()
