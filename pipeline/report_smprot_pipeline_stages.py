"""
Summarize SmProt-derived peptide processing: counts kept / excluded / FASTA failures.

Writes:

  - ``data/reports/smprot_pipeline_stages.json`` — machine-readable
  - ``data/reports/smprot_pipeline_stages.md`` — human-readable tables (regenerate anytime)

Run after ``build_significant_lncs_smprot.py``, ``filter_peptides_tcga_expr_genes.py``,
and FASTA export + sync scripts so on-disk artifacts match the report.

**Canonical** appears only for the *significant-lncRNA* branch: it means peptides whose
``GeneID`` / ``GeneSymbol`` matched ``significant_lncs.csv`` (see ``build_significant_lncs_smprot.py``).
It is unrelated to the TCGA expression-matrix column filter.
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


import json
from collections import Counter
from datetime import datetime, timezone

import pandas as pd

DATA = ROOT / "data"
REPORTS = DATA / "reports"
LNCRNA_PEP = DATA / "lncRNA_peptides.tsv"
SMPROT_FILT = DATA / "smprot_filtered.tsv"
TCGA_TSV = DATA / "smprot_filtered_tcga_expr_genes.tsv"
ALL_FAA = DATA / "smprot_all_filtered_peptides.faa"
TCGA_FAA = DATA / "smprot_tcga_filtered_peptides.faa"
SIG_FULL = DATA / "significant_lnc_peptides_full.tsv"
SIG_FAA = DATA / "significant_lnc_peptides.faa"
SIG_TSV = DATA / "significant_lnc_peptides.tsv"
SIG_LNCS = DATA / "significant_lncs.csv"


def count_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    n = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for _ in fh:
            n += 1
    return n


def fasta_record_count(path: Path) -> int | None:
    if not path.exists():
        return None
    c = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                c += 1
    return c


def fasta_smpep_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    if not path.exists():
        return ids
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split("|", 1)[0].strip())
    return ids


def tsv_smpep_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", usecols=["smPEP_ID"], dtype=str, low_memory=False)
    return set(df["smPEP_ID"].astype(str).str.strip())


def tsv_row_count(path: Path) -> int | None:
    if not path.exists():
        return None
    return len(pd.read_csv(path, sep="\t", usecols=["smPEP_ID"], dtype=str, low_memory=False))


def build_report() -> dict:
    out: dict = {"root": str(ROOT), "stages": []}

    # --- Stage 0: raw lncRNA peptide dump (optional) ---
    lnc_lines = count_lines(LNCRNA_PEP)
    lnc_data = (lnc_lines - 1) if lnc_lines is not None and lnc_lines > 0 else None
    out["stages"].append(
        {
            "id": "0_lncrna_peptides_tsv",
            "title": "Intermediate: all lncRNA rows from SmProt2",
            "artifacts": {"tsv": str(LNCRNA_PEP.relative_to(ROOT))},
            "script": "build_significant_lncs_smprot.py (pass 1)",
            "data_rows": lnc_data,
            "notes": "Unfiltered lncRNA GeneType rows copied from SmProt2.txt for downstream filtering.",
        }
    )

    # --- Stage 1: smprot_filtered.tsv ---
    n_filt = tsv_row_count(SMPROT_FILT)
    filt_ids = tsv_smpep_ids(SMPROT_FILT)
    out["stages"].append(
        {
            "id": "1_smprot_filtered",
            "title": "SmProt-filtered lncRNA MPs (all genes, FASTA-synced)",
            "artifacts": {"tsv": str(SMPROT_FILT.relative_to(ROOT))},
            "script": "build_significant_lncs_smprot.py then export_tcga_filtered_peptides_fasta.py on smprot_filtered.tsv",
            "data_rows": n_filt,
            "unique_smPEP_ID": len(filt_ids),
            "rules": (
                "Build: lncRNA; TIS & Ribo p < 0.05; ORF >= 10 aa; dedupe smPEP_ID. "
                "Export (all-filtered FASTA): rewrites this TSV to rows with a successful FASTA translation only."
            ),
        }
    )

    # --- Stage 2: all-filtered FASTA ---
    n_all_faa = fasta_record_count(ALL_FAA)
    all_faa_ids = fasta_smpep_ids(ALL_FAA)
    failed_all = filt_ids - all_faa_ids if filt_ids else set()
    out["stages"].append(
        {
            "id": "2_all_filtered_fasta",
            "title": "Translated FASTA for all filtered MPs",
            "artifacts": {"faa": str(ALL_FAA.relative_to(ROOT))},
            "script": "export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa",
            "fasta_records": n_all_faa,
            "in_tsv_not_in_fasta": len(failed_all),
            "notes": (
                "After FASTA-synced smprot_filtered.tsv, this count should match the TSV row count "
                "(in_tsv_not_in_fasta → 0). Non-zero means the FASTA is stale vs the TSV."
            ),
        }
    )

    # --- Stage 3: TCGA matrix–mapped TSV ---
    n_tcga = tsv_row_count(TCGA_TSV)
    tcga_ids = tsv_smpep_ids(TCGA_TSV)
    excluded_tcga = filt_ids - tcga_ids if filt_ids and tcga_ids else set()
    match_via: dict[str, int] = {}
    if TCGA_TSV.exists() and "TCGA_match_via" in pd.read_csv(TCGA_TSV, sep="\t", nrows=0).columns:
        dfm = pd.read_csv(TCGA_TSV, sep="\t", usecols=["TCGA_match_via"], dtype=str, low_memory=False)
        match_via = dict(Counter(dfm["TCGA_match_via"].fillna("unknown").astype(str)))

    out["stages"].append(
        {
            "id": "3_tcga_matrix_tsv",
            "title": "TCGA expression–matrix mapped MPs",
            "artifacts": {"tsv": str(TCGA_TSV.relative_to(ROOT))},
            "script": "filter_peptides_tcga_expr_genes.py",
            "data_rows": n_tcga,
            "unique_smPEP_ID": len(tcga_ids),
            "excluded_not_in_matrix": len(excluded_tcga),
            "TCGA_match_via_counts": match_via,
            "rules": "GeneSymbol in matrix column names OR GeneID→lncrna_genes_small gene_name in columns; dedupe smPEP_ID after sort.",
        }
    )

    # --- Stage 4: TCGA FASTA ---
    n_tcga_faa = fasta_record_count(TCGA_FAA)
    tcga_faa_ids = fasta_smpep_ids(TCGA_FAA)
    failed_tcga = tcga_ids - tcga_faa_ids if tcga_ids else set()
    out["stages"].append(
        {
            "id": "4_tcga_filtered_fasta",
            "title": "Translated FASTA for TCGA-mapped MPs",
            "artifacts": {"faa": str(TCGA_FAA.relative_to(ROOT))},
            "script": "export_tcga_filtered_peptides_fasta.py (defaults)",
            "fasta_records": n_tcga_faa,
            "in_tsv_not_in_fasta": len(failed_tcga),
        }
    )

    # --- Stage 5: significant full (canonical *gene* list match) ---
    n_sig_full = tsv_row_count(SIG_FULL)
    sig_full_ids = tsv_smpep_ids(SIG_FULL)
    overlap_sig_in_all = len(sig_full_ids & filt_ids) if sig_full_ids and filt_ids else None
    out["stages"].append(
        {
            "id": "5_significant_full_tsv",
            "title": "Peptides on canonical significant lncRNA *genes* (SmProt table)",
            "artifacts": {"tsv": str(SIG_FULL.relative_to(ROOT))},
            "script": "build_significant_lncs_smprot.py",
            "data_rows": n_sig_full,
            "unique_smPEP_ID": len(sig_full_ids),
            "smPEP_also_in_smprot_filtered": overlap_sig_in_all,
            "gene_list": str(SIG_LNCS.relative_to(ROOT)) if SIG_LNCS.exists() else None,
            "rules": "Subset of smprot_filtered rows whose GeneID/GeneSymbol matches significant_lncs.csv (limma DE ∩ z union genes + ENSG mapping). Not filtered by TCGA matrix columns.",
        }
    )

    # --- Stage 6: significant FASTA + analyzed TSV ---
    n_sig_faa = fasta_record_count(SIG_FAA)
    sig_faa_ids = fasta_smpep_ids(SIG_FAA)
    n_sig_tsv = tsv_row_count(SIG_TSV)
    sig_tsv_ids = tsv_smpep_ids(SIG_TSV)
    failed_sig_faa = sig_full_ids - sig_faa_ids if sig_full_ids and sig_faa_ids else set()
    excluded_from_analyzed = sig_full_ids - sig_tsv_ids if sig_full_ids and sig_tsv_ids else set()
    out["stages"].append(
        {
            "id": "6_significant_fasta_and_analyzed",
            "title": "Significant: exportable FASTA + FASTA-synced analyzed TSV",
            "artifacts": {
                "faa": str(SIG_FAA.relative_to(ROOT)),
                "tsv": str(SIG_TSV.relative_to(ROOT)),
            },
            "script": "export_tcga_filtered_peptides_fasta.py (--peptides-tsv significant_lnc_peptides_full.tsv) then sync_significant_lnc_peptides_with_fasta.py",
            "fasta_records": n_sig_faa,
            "analyzed_tsv_rows": n_sig_tsv,
            "in_full_not_in_fasta": len(failed_sig_faa),
            "in_full_not_in_analyzed_tsv": len(excluded_from_analyzed),
            "notes": "analyzed TSV keeps only rows whose smPEP_ID appears in the FASTA (translatable ORFs).",
        }
    )

    return out


def _md_table_cell(s: str) -> str:
    """Escape characters that break pipe tables."""
    return str(s).replace("|", " / ").replace("\n", " ")


def snapshot_table_lines(data: dict) -> list[str]:
    """Flatten key counts into one markdown table (actual numbers from current files)."""
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "## Snapshot (counts from current files)",
        "",
        f"*Generated: {ts}*",
        "",
        "| Stage | Metric | Value |",
        "|-------|--------|-------|",
    ]
    for st in data.get("stages", []):
        sid = str(st.get("id", ""))
        title = _md_table_cell(st.get("title", sid))
        for key in (
            "data_rows",
            "unique_smPEP_ID",
            "fasta_records",
            "in_tsv_not_in_fasta",
            "excluded_not_in_matrix",
            "analyzed_tsv_rows",
            "in_full_not_in_fasta",
            "in_full_not_in_analyzed_tsv",
            "smPEP_also_in_smprot_filtered",
        ):
            if key in st and st[key] is not None:
                lines.append(f"| {title} | `{key}` | {st[key]} |")
        m = st.get("TCGA_match_via_counts")
        if isinstance(m, dict):
            for k, v in sorted(m.items(), key=lambda kv: (-int(kv[1]), kv[0])):
                lines.append(f"| TCGA matrix TSV | `TCGA_match_via` = `{k}` | {v} |")
    lines.append("")
    return lines


def to_markdown(data: dict) -> str:
    lines: list[str] = [
        "# SmProt peptide pipeline — stage counts",
        "",
        "Regenerate: `python report_smprot_pipeline_stages.py`",
        "",
        "**Canonical** (significant branch only): gene membership in `significant_lncs.csv` from `build_significant_lncs_smprot.py` — *not* the TCGA matrix filter.",
        "",
    ]
    lines.extend(snapshot_table_lines(data))
    for st in data.get("stages", []):
        lines.append(f"## {st.get('title', st.get('id'))}")
        lines.append("")
        lines.append(f"- **Stage id:** `{st.get('id')}`")
        if st.get("script"):
            lines.append(f"- **Script(s):** {st['script']}")
        for k, v in st.items():
            if k in ("id", "title", "script", "notes", "rules"):
                continue
            if v is None:
                continue
            if isinstance(v, dict):
                if k == "artifacts":
                    lines.append("- **Artifacts:**")
                elif k == "TCGA_match_via_counts":
                    lines.append("- **TCGA_match_via (how rows matched the matrix):**")
                else:
                    lines.append(f"- **{k}:**")
                for sk, sv in v.items():
                    lines.append(f"  - `{sk}`: {sv}")
            else:
                lines.append(f"- **{k}:** {v}")
        if st.get("rules"):
            lines.append(f"- **Rules:** {st['rules']}")
        if st.get("notes"):
            lines.append(f"- **Notes:** {st['notes']}")
        lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    data = build_report()
    REPORTS.mkdir(parents=True, exist_ok=True)
    jpath = REPORTS / "smprot_pipeline_stages.json"
    mpath = REPORTS / "smprot_pipeline_stages.md"
    jpath.write_text(json.dumps(data, indent=2), encoding="utf-8")
    mpath.write_text(to_markdown(data), encoding="utf-8")
    print(f"Wrote {jpath}")
    print(f"Wrote {mpath}")


if __name__ == "__main__":
    main()
