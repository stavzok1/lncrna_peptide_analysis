"""
Export protein (and optionally DNA) FASTA for rows in a SmProt peptide TSV.

Default input is ``data/smprot_filtered_tcga_expr_genes.tsv`` (TCGA-matrix genes).
Use ``--peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa``
to export **all** curated SmProt-filtered peptides (for Fig. 2 / 3B ``--peptide-set all_smprot_filtered``).
That run **also** rewrites ``smprot_filtered.tsv`` to **only** rows that produced a FASTA record
(FASTA success is part of the filter), then refreshes ``significant_lnc_peptides_full.tsv`` via
``smprot_gene_match.py``, and by default re-exports ``significant_lnc_peptides.faa`` and runs
``sync_significant_lnc_peptides_with_fasta.py`` (disable with ``--no-chain-significant-fasta``).

Use ``--peptides-tsv data/significant_lnc_peptides_full.tsv`` to export only the significant list
(writes ``data/significant_lnc_peptides.faa`` when ``--out-aa`` follows the default naming rule).
Then run ``sync_significant_lnc_peptides_with_fasta.py`` if you did not use the all-filtered chain above.

Sequences come only from a **local** GENCODE transcript FASTA (spliced cDNA),
by slicing ``full[StartOnTrans:StopOnTrans]`` on the matching ENST sequence from
GENCODE **transcript** FASTA. Indexing is **0-based, half-open** (Python slice):
nucleotide positions ``StartOnTrans`` inclusive through ``StopOnTrans - 1``
inclusive are translated; ``StopOnTrans`` is exclusive. The first codon uses
bases ``full[s], full[s+1], full[s+2]``. SmProt coordinates are interpreted in
that same frame on the sequence **as stored** in the FASTA (transcript 5'→3';
no reverse-complement is applied in this script).

Why not genome FASTA alone?
  Peptide coordinates are on the **transcript** (spliced RNA). Rebuilding that
  from genomic DNA requires exon structure (GTF) and splicing logic. The
  transcript FASTA already matches SmProt's TranscriptID + coordinates.

Why Ensembl REST failed before?
  - Versioned ENST IDs (e.g. *.6) are often absent from the current REST API.
  - urllib quoting turned '.' into %2E (400 errors) in an earlier version.

Default transcript path points at your WSL GENCODE v49 file; override with
  --transcripts-fa

Optional: --ensembl-fallback uses network only for rows still missing after
  the local scan (off by default).
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
import gzip
import json
import subprocess
import time
import urllib.error
import urllib.request
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

DATA = ROOT / "data"
TCGA_PEP = DATA / "smprot_filtered_tcga_expr_genes.tsv"
SMPROT_FILT = DATA / "smprot_filtered.tsv"
SIG_PEPT_FULL = DATA / "significant_lnc_peptides_full.tsv"
MIN_AA_LENGTH = 10
OUT_AA_DEFAULT = DATA / "smprot_tcga_filtered_peptides.faa"
OUT_DNA_DEFAULT = DATA / "smprot_tcga_filtered_peptides_nt.fna"
CDNA_CACHE = DATA / "ensembl_cdna_cache.json"

# Local GENCODE transcriptome (UNC path to WSL file — override with --transcripts-fa).
DEFAULT_TRANSCRIPTS_FA = Path(
    r"\\wsl.localhost\Ubuntu\home\stavz\masters\gdc\APM\data\gencode.v49.transcripts.fa.gz"
)


def header_to_transcript_id(desc: str) -> str:
    """First field of FASTA header (e.g. ENST...|ENSG...|...)."""
    line = desc.splitlines()[0].lstrip(">")
    return line.split("|", 1)[0].split()[0]


def load_disk_cache() -> dict[str, str]:
    if not CDNA_CACHE.exists():
        return {}
    try:
        raw = json.loads(CDNA_CACHE.read_text(encoding="utf-8"))
        return {str(k): str(v) for k, v in raw.items() if isinstance(v, str) and v}
    except (json.JSONDecodeError, OSError):
        return {}


def save_disk_cache(cache: dict[str, str]) -> None:
    CDNA_CACHE.write_text(json.dumps(cache, indent=0), encoding="utf-8")


def fetch_cdna_ensembl(
    ensembl_id: str,
    sleep_s: float,
    disk_cache: dict[str, str],
) -> Optional[str]:
    candidates = [ensembl_id]
    if "." in ensembl_id:
        b = ensembl_id.rsplit(".", 1)[0]
        if b not in candidates:
            candidates.append(b)
    for tid in candidates:
        if tid in disk_cache:
            return disk_cache[tid]
        url = f"https://rest.ensembl.org/sequence/id/{tid}?type=cdna"
        req = urllib.request.Request(
            url,
            headers={
                "Accept": "text/plain",
                "User-Agent": "UNDEFINED-tcga-peptide-fasta-export (local research)",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                body = resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                continue
            raise
        lines = body.strip().splitlines()
        if not lines:
            continue
        seq = (
            "".join(lines[1:]) if lines[0].startswith(">") else "".join(lines)
        ).replace(" ", "").upper()
        if seq:
            disk_cache[tid] = seq
            time.sleep(sleep_s)
            return seq
    return None


def translate_dna(dna: str) -> str:
    dna = dna.upper().replace("U", "T")
    if len(dna) % 3 != 0:
        dna = dna[: len(dna) - (len(dna) % 3)]
    aa = str(Seq(dna).translate())
    return aa.split("*")[0]


def stream_transcript_fasta(
    fa_path: Path,
    needed_ids: set[str],
) -> tuple[dict[str, str], dict[str, list[str]]]:
    """
    One pass over GENCODE transcripts.fa(.gz).

    Stores every transcript whose id is in needed_ids OR whose versionless
    ENST prefix matches any needed id (so we can resolve single-isoform
    fallbacks when SmProt uses an older patch version not in this release).
    """
    needed_bases = {t.rsplit(".", 1)[0] for t in needed_ids if "." in t}

    if not fa_path.exists():
        raise FileNotFoundError(
            f"Transcript FASTA not found: {fa_path}\n"
            "Pass --transcripts-fa /path/to/gencode.*.transcripts.fa.gz"
        )

    tx: dict[str, str] = {}
    by_base: dict[str, list[str]] = {}

    opener = gzip.open if fa_path.suffix == ".gz" else open
    mode = "rt"
    encoding = "utf-8"
    with opener(fa_path, mode, encoding=encoding) as fh:  # type: ignore[arg-type]
        for rec in SeqIO.parse(fh, "fasta"):
            rid = header_to_transcript_id(rec.description)
            seq = str(rec.seq).upper().replace("U", "T")
            if "." in rid:
                base = rid.rsplit(".", 1)[0]
                if rid in needed_ids or base in needed_bases:
                    tx[rid] = seq
                    by_base.setdefault(base, []).append(rid)
            elif rid in needed_ids:
                tx[rid] = seq

    return tx, by_base


def get_full_sequence(tx: dict[str, str], by_base: dict[str, list[str]], tid: str) -> Optional[str]:
    if tid in tx:
        return tx[tid]
    if "." not in tid:
        return None
    base = tid.rsplit(".", 1)[0]
    alts = sorted(set(by_base.get(base, [])))
    if len(alts) == 1:
        return tx.get(alts[0])
    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description="FASTA for SmProt peptide rows (GENCODE transcript slice + translate)."
    )
    ap.add_argument("--dna", action="store_true", help="Also write nucleotide FASTA.")
    ap.add_argument("--out-aa", type=Path, default=OUT_AA_DEFAULT)
    ap.add_argument("--out-dna", type=Path, default=OUT_DNA_DEFAULT)
    ap.add_argument(
        "--peptides-tsv",
        type=Path,
        default=None,
        help=(
            "Peptide table TSV (columns smPEP_ID, GeneSymbol, GeneID, TranscriptID, "
            "StartOnTrans, StopOnTrans). Default: smprot_filtered_tcga_expr_genes.tsv."
        ),
    )
    ap.add_argument(
        "--transcripts-fa",
        type=Path,
        default=DEFAULT_TRANSCRIPTS_FA,
        help="GENCODE transcripts FASTA (.fa or .fa.gz).",
    )
    ap.add_argument(
        "--ensembl-fallback",
        action="store_true",
        help="If local FASTA has no transcript, try Ensembl REST (needs internet).",
    )
    ap.add_argument("--ensembl-sleep", type=float, default=0.08)
    ap.add_argument(
        "--min-aa-length",
        type=int,
        default=MIN_AA_LENGTH,
        help=f"Minimum translated peptide length to emit (default: {MIN_AA_LENGTH}). "
        "Use 9 for NetMHC 9-mer parents if you accept shorter ORFs than the SmProt ≥10 aa filter.",
    )
    ap.add_argument(
        "--no-restrict-smprot-filtered-to-fasta",
        action="store_true",
        help="When --peptides-tsv is data/smprot_filtered.tsv, do not rewrite that TSV to FASTA-only rows.",
    )
    ap.add_argument(
        "--no-chain-significant-fasta",
        action="store_true",
        help="When restricting smprot_filtered.tsv, skip re-export of significant_lnc_peptides.faa and sync.",
    )
    args = ap.parse_args()

    pep_path = args.peptides_tsv if args.peptides_tsv is not None else TCGA_PEP
    out_aa = args.out_aa
    out_dna = args.out_dna
    if args.peptides_tsv is not None:
        if out_aa == OUT_AA_DEFAULT:
            p = args.peptides_tsv
            if p.stem == "significant_lnc_peptides_full":
                out_aa = p.with_name("significant_lnc_peptides.faa")
            else:
                out_aa = p.with_suffix(".faa")
        if out_dna == OUT_DNA_DEFAULT:
            p = args.peptides_tsv
            if p.stem == "significant_lnc_peptides_full":
                out_dna = p.with_name("significant_lnc_peptides_nt.fna")
            else:
                out_dna = p.with_name(f"{p.stem}_nt.fna")

    tc = pd.read_csv(pep_path, sep="\t", dtype=str)
    need_cols = {"smPEP_ID", "GeneSymbol", "GeneID", "TranscriptID", "StartOnTrans", "StopOnTrans"}
    if not need_cols <= set(tc.columns):
        raise SystemExit(f"TCGA TSV missing columns: {need_cols - set(tc.columns)}")

    needed_ids = set(tc["TranscriptID"].astype(str).str.strip())
    print(f"Scanning transcript FASTA (need {len(needed_ids)} ENST labels)...")
    tx, by_base = stream_transcript_fasta(args.transcripts_fa, needed_ids)
    print(f"Loaded transcript sequences in memory: {len(tx)} FASTA records")

    disk_cache = load_disk_cache() if args.ensembl_fallback else {}
    n_local = n_fallback = n_fail = 0

    aa_lines: list[str] = []
    dna_lines: list[str] = []
    min_aa = max(1, int(args.min_aa_length))
    written_smpep: set[str] = set()

    for _, row in tc.iterrows():
        sid = str(row["smPEP_ID"])
        sym = str(row["GeneSymbol"])
        gid = str(row["GeneID"])
        tid = str(row["TranscriptID"])
        try:
            s = int(row["StartOnTrans"])
            e = int(row["StopOnTrans"])
        except (TypeError, ValueError):
            n_fail += 1
            continue

        def slice_translate(cdna: str) -> Optional[str]:
            if s < 0 or e > len(cdna) or s >= e:
                return None
            if (e - s) // 3 < min_aa:
                return None
            dna = cdna[s:e].upper()
            try:
                aa = translate_dna(dna)
            except Exception:
                return None
            if len(aa) < min_aa:
                return None
            return aa

        full = get_full_sequence(tx, by_base, tid)
        src = "gencode_transcript"
        aa: Optional[str] = slice_translate(full) if full is not None else None

        # Local transcript can be wrong isoform / too short; still try Ensembl when enabled.
        if aa is None and args.ensembl_fallback:
            if tid not in disk_cache:
                disk_cache[tid] = fetch_cdna_ensembl(tid, args.ensembl_sleep, disk_cache)  # type: ignore[assignment]
            full_e = disk_cache.get(tid)
            if full_e is not None:
                aa2 = slice_translate(full_e)
                if aa2:
                    aa = aa2
                    src = "ensembl_cdna"

        if aa is None:
            n_fail += 1
            continue

        if src == "gencode_transcript":
            n_local += 1
        else:
            n_fallback += 1

        hdr = f">{sid}|{sym}|{gid}|{tid}|{s}-{e}|{src}"
        aa_lines.append(hdr)
        aa_lines.append(aa)
        written_smpep.add(sid)
        if args.dna:
            dna_lines.append(hdr)
            dna_lines.append(dna)

    out_aa.write_text("\n".join(aa_lines) + "\n", encoding="utf-8")
    if args.dna and dna_lines:
        out_dna.write_text("\n".join(dna_lines) + "\n", encoding="utf-8")

    if args.ensembl_fallback:
        save_disk_cache(disk_cache)

    n_written = len(aa_lines) // 2
    print(f"Rows in {pep_path.name}: {len(tc)}")
    print(f"FASTA records written: {n_written} -> {out_aa}")
    if args.dna:
        print(f"DNA FASTA records: {len(dna_lines)//2} -> {out_dna}")
    print(f"Resolved from local GENCODE (exact or 1-isoform fallback): {n_local}")
    if args.ensembl_fallback:
        print(f"Resolved via Ensembl fallback: {n_fallback}")
    print(f"Skipped / failed: {n_fail}")
    if n_fail and not args.ensembl_fallback:
        print(
            "Tip: use a GENCODE release that matches SmProt ENST versions, or "
            "--ensembl-fallback for missing IDs."
        )

    # --- FASTA as filtration for all-filtered SmProt table + downstream significant chain ---
    if (
        pep_path.resolve() == SMPROT_FILT.resolve()
        and not args.no_restrict_smprot_filtered_to_fasta
    ):
        n_before = len(tc)
        tc_ok = tc.loc[tc["smPEP_ID"].astype(str).str.strip().isin(written_smpep)].copy()
        tc_ok = tc_ok[list(tc.columns)]
        tc_ok.to_csv(SMPROT_FILT, sep="\t", index=False)
        print(
            f"Rewrote {SMPROT_FILT.name} to {len(tc_ok)} FASTA-exportable rows "
            f"(removed {n_before - len(tc_ok)} without a successful translation)."
        )
        from smprot_gene_match import write_significant_lnc_peptides_full

        write_significant_lnc_peptides_full(SMPROT_FILT)

        if not args.no_chain_significant_fasta and SIG_PEPT_FULL.exists():
            sig_cmd = [
                sys.executable,
                str(ROOT / "export_tcga_filtered_peptides_fasta.py"),
                "--peptides-tsv",
                str(SIG_PEPT_FULL),
                "--transcripts-fa",
                str(args.transcripts_fa),
                "--no-restrict-smprot-filtered-to-fasta",
                "--no-chain-significant-fasta",
            ]
            if args.ensembl_fallback:
                sig_cmd.append("--ensembl-fallback")
            if args.dna:
                sig_cmd.append("--dna")
            print("+", " ".join(sig_cmd))
            r1 = subprocess.run(sig_cmd, cwd=str(ROOT))
            if r1.returncode != 0:
                print(f"Warning: significant FASTA export exited {r1.returncode}", file=sys.stderr)
            sync_cmd = [sys.executable, str(ROOT / "sync_significant_lnc_peptides_with_fasta.py")]
            print("+", " ".join(sync_cmd))
            r2 = subprocess.run(sync_cmd, cwd=str(ROOT))
            if r2.returncode != 0:
                print(f"Warning: sync_significant_lnc_peptides_with_fasta.py exited {r2.returncode}", file=sys.stderr)


if __name__ == "__main__":
    main()
