"""
Prepare NetMHCpan-4.1 inputs: sliding-window 9-mers from (1) all significant
lncRNA SmProt peptides (``data/significant_lnc_peptides.tsv``, 549 ``smPEP_ID``)
and (2) length-matched coding **controls** sampled from a human coding proteome FASTA.

By default, controls are **contiguous fragments** of that length cut from longer proteins
(common for short micropeptides: full UniProt entries rarely match e.g. 12 aa exactly).
Use **``--coding-control-mode whole_protein``** to instead use **whole FASTA entries**
whose sequence length **equals** each required length (needs enough exact-length proteins
in the FASTA; often fails for very short lengths).

Use **``--coding-control-mode proportional_whole``** to sample **whole** proteome proteins
into **length bins** (default equal-width over the reference length range) so the
**fraction of parents per bin** matches the Tr-lncRNA (reference) micropeptide set as
closely as possible for a target parent count **``--proportional-target-n``** (default =
number of reference parents). Phase **1** uses largest-remainder targets per bin; an optional
**top-up** pass (default **on**) draws additional proteins from bins that still have unused
candidates, weighted by reference bin fraction, until ``target_n`` is reached or pools are
exhausted (``--no-proportional-top-up`` to disable). Shortfalls are reported in
``coding_proportional_bin_summary.csv``.

Primary sequences must exist in ``data/significant_lnc_peptides.faa``. Export from
the **full** SmProt table, then sync the analyzed TSV::

  python export_tcga_filtered_peptides_fasta.py --peptides-tsv data/significant_lnc_peptides_full.tsv \\
    --transcripts-fa /path/to/gencode.*.transcripts.fa.gz [--ensembl-fallback] [--min-aa-length 9]
  python sync_significant_lnc_peptides_with_fasta.py

Does NOT run NetMHCpan (requires local Linux/WSL install). Writes FASTA of
9-mers and metadata under ``data/netmhc/``.

Usage:
  python prepare_netmhc_tr_vs_coding_epitopes.py
  python prepare_netmhc_tr_vs_coding_epitopes.py --seed 42 --coding-fa data/known_proteins.fasta
  python prepare_netmhc_tr_vs_coding_epitopes.py --coding-control-mode proportional_whole --max-proteins 0
    (reference length bins default to ~501 significant MPs: significant_lnc_peptides.faa + .tsv)
  python prepare_netmhc_tr_vs_coding_epitopes.py --coding-control-mode proportional_whole \\
    --peptide-faa data/smprot_tcga_filtered_peptides.faa --peptide-tsv data/smprot_filtered_tcga_expr_genes.tsv \\
    --max-proteins 0
    (optional: ~2.6k TCGA-matrix-filtered parents as reference — not the 501 significant cohort)
"""
from __future__ import annotations

import argparse
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = DATA / "netmhc"
SIG_TSV = DATA / "significant_lnc_peptides.tsv"
SIG_FAA = DATA / "significant_lnc_peptides.faa"
CODING_FA_DEFAULT = DATA / "known_proteins.fasta"

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def clean_aa(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in AA20)


def load_expected_smpep_ids(tsv_path: Path) -> set[str]:
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    if "smPEP_ID" not in df.columns:
        raise SystemExit(f"{tsv_path} has no smPEP_ID column")
    return set(df["smPEP_ID"].astype(str).str.strip())


def parse_sig_peptide_fasta(path: Path) -> list[tuple[str, str, str]]:
    """List of (smPEP_ID, GeneSymbol, aa_sequence) from export header >id|sym|..."""
    out: list[tuple[str, str, str]] = []
    for rec in SeqIO.parse(path, "fasta"):
        parts = rec.description.split("|")
        sid = parts[0].strip() if parts else rec.id
        sym = parts[1].strip() if len(parts) > 1 else ""
        seq = clean_aa(str(rec.seq))
        if seq:
            out.append((sid, sym, seq))
    return out


def sliding_ninemers(seq: str, parent_id: str, source: str) -> list[tuple[str, str]]:
    if len(seq) < 9:
        return []
    rows: list[tuple[str, str]] = []
    for i in range(0, len(seq) - 8):
        mer = seq[i : i + 9]
        uid = f"{source}|{parent_id}|{i}|{mer}"
        rows.append((uid, mer))
    return rows


def sample_whole_protein_matched_coding_mps(
    length_hist: Counter[int],
    coding_fa: Path,
    rng: random.Random,
    max_proteins_scan: int,
) -> list[tuple[str, str, str]]:
    """
    For each length L needed, pick ``n`` distinct FASTA records whose **full**
    cleaned sequence length is exactly ``L`` (whole entry, not a substring).

    Stops scanning early once every length bucket has at least ``n`` candidates
    (unless ``max_proteins_scan`` caps the scan first).
    """
    need: dict[int, int] = {L: int(c) for L, c in length_hist.items() if c > 0}
    lengths_needed = frozenset(need.keys())
    bucket: dict[int, list[tuple[str, str]]] = {L: [] for L in need}

    def sufficient() -> bool:
        return all(len(bucket[L]) >= need[L] for L in need)

    idx = SeqIO.index(str(coding_fa), "fasta")
    scanned = 0
    try:
        ids = list(idx.keys())
        rng.shuffle(ids)
        for raw_id in ids:
            scanned += 1
            if max_proteins_scan > 0 and scanned > max_proteins_scan:
                break
            rec = idx[raw_id]
            seq = clean_aa(str(rec.seq))
            L = len(seq)
            if L not in lengths_needed:
                continue
            pid = raw_id.replace(" ", "_")[:120]
            bucket[L].append((pid, seq))
            if sufficient():
                break
        if not sufficient():
            miss = {L: max(0, need[L] - len(bucket[L])) for L in need}
            hint = (
                "Very short micropeptide lengths rarely occur as full-length entries in a "
                "canonical human proteome FASTA. Options: use default fragment controls "
                "(omit --coding-control-mode whole_protein), add a peptide/sORF database, "
                "or use --max-proteins 0 and a larger FASTA."
            )
            raise RuntimeError(
                f"Whole-protein control sampling: insufficient exact-length proteins after "
                f"scanning {scanned} FASTA record(s). Shortfall per length (need - found): {miss}. {hint}"
            )
        out: list[tuple[str, str, str]] = []
        for L in sorted(need.keys()):
            n = need[L]
            chosen = rng.sample(bucket[L], n)
            for i, (pid, seq) in enumerate(chosen):
                oid = f"coding_whole|{pid}|L{L}|{i}"
                out.append((oid, "CODING", seq))
        return out
    finally:
        idx.close()


def sample_length_matched_coding_mps(
    length_hist: Counter[int],
    coding_fa: Path,
    rng: random.Random,
    max_proteins_scan: int,
) -> list[tuple[str, str, str]]:
    """
    Same number of parent MPs per length as the significant set, each a
    contiguous substring drawn from a longer proteome sequence. Uses Bio.SeqIO.index
    so the full FASTA is not loaded into RAM.
    """
    need: dict[int, int] = {L: int(c) for L, c in length_hist.items() if c > 0}
    out: list[tuple[str, str, str]] = []
    idx = SeqIO.index(str(coding_fa), "fasta")
    try:
        ids = list(idx.keys())
        rng.shuffle(ids)
        scanned = 0
        for raw_id in ids:
            scanned += 1
            if max_proteins_scan > 0 and scanned > max_proteins_scan and sum(need.values()) > 0:
                break
            rec = idx[raw_id]
            seq = clean_aa(str(rec.seq))
            min_need_L = min(need.keys(), default=0)
            if not need or len(seq) < min_need_L:
                continue
            pid = raw_id.replace(" ", "_")[:120]
            for L in sorted(need.keys(), reverse=True):
                if need.get(L, 0) <= 0 or len(seq) < L:
                    continue
                starts = list(range(len(seq) - L + 1))
                rng.shuffle(starts)
                for st in starts:
                    if need[L] <= 0:
                        break
                    frag = seq[st : st + L]
                    if len(frag) != L:
                        continue
                    oid = f"coding|{pid}|{st}|L{L}|{len(out)}"
                    out.append((oid, "CODING", frag))
                    need[L] -= 1
            if sum(need.values()) == 0:
                break
        if sum(need.values()) > 0:
            hint = (
                "Increase --max-proteins (0 = no cap), use a larger proteome FASTA "
                "(e.g. UniProt human reference), or check rare peptide lengths."
            )
            raise RuntimeError(
                f"Could not fill all length buckets (scanned {scanned} proteome records). "
                f"Remaining counts: {need}. {hint}"
            )
    finally:
        idx.close()
    return out


def _bin_key(lo: int, hi: int) -> str:
    return f"{lo}-{hi}"


def length_to_bin_edges(L: int, l_min: int, bin_width: int) -> tuple[int, int]:
    """Inclusive [lo, hi] bin for length L using bins aligned from l_min."""
    if L < l_min:
        raise ValueError(f"length {L} < l_min {l_min}")
    k = (L - l_min) // bin_width
    lo = l_min + k * bin_width
    hi = lo + bin_width - 1
    return lo, hi


def ref_lengths_to_bins(lengths_counter: Counter[int], bin_width: int) -> tuple[dict[str, int], int, int]:
    """
    Aggregate reference parent MP counts by equal-width bins spanning
    [min(length), max(length)] among lengths with count > 0.
    """
    active = [L for L, c in lengths_counter.items() if c > 0]
    if not active:
        raise ValueError("empty length histogram")
    l_min, l_max = min(active), max(active)
    bin_counts: dict[str, int] = defaultdict(int)
    for L, c in lengths_counter.items():
        if c <= 0:
            continue
        lo, hi = length_to_bin_edges(L, l_min, bin_width)
        bin_counts[_bin_key(lo, hi)] += int(c)
    return dict(bin_counts), l_min, l_max


def largest_remainder_allocation(bin_counts: dict[str, int], M: int) -> dict[str, int]:
    """Integer targets per bin summing to M, proportional to bin_counts."""
    tot = sum(bin_counts.values())
    if tot <= 0:
        raise ValueError("bin_counts sum to zero")
    keys = list(bin_counts.keys())
    raw = {k: M * bin_counts[k] / tot for k in keys}
    base = {k: int(math.floor(raw[k])) for k in keys}
    rem = M - sum(base.values())
    if rem < 0:
        for k in sorted(keys, key=lambda x: base[x] - raw[x], reverse=True):
            while rem < 0 and base[k] > 0:
                base[k] -= 1
                rem += 1
    elif rem > 0:
        order = sorted(keys, key=lambda k: raw[k] - base[k], reverse=True)
        i = 0
        while rem > 0:
            base[order[i % len(order)]] += 1
            rem -= 1
            i += 1
    return base


def sample_proportional_whole_proteins(
    lengths_counter: Counter[int],
    coding_fa: Path,
    rng: random.Random,
    *,
    bin_width: int,
    target_n: int,
    max_proteins_scan: int,
    min_protein_len: int,
    top_up: bool = True,
) -> tuple[list[tuple[str, str, str]], pd.DataFrame]:
    """
    Build length bins from the reference histogram; pick ``target_n`` whole proteome
    sequences whose lengths fall into bins with counts matching reference bin fractions.

    Returns (coding_mps same schema as other samplers, summary dataframe).
    """
    if target_n < 1:
        raise ValueError("target_n must be >= 1")
    bin_counts, l_min, l_max = ref_lengths_to_bins(lengths_counter, bin_width)
    targets = largest_remainder_allocation(bin_counts, target_n)

    def bin_for_len(L: int) -> str:
        if L < l_min or L > l_max:
            return ""
        lo, hi = length_to_bin_edges(L, l_min, bin_width)
        return _bin_key(lo, hi)

    pool: dict[str, list[tuple[str, str]]] = {k: [] for k in bin_counts}
    scanned = 0
    idx = SeqIO.index(str(coding_fa), "fasta")
    try:
        ids = list(idx.keys())
        rng.shuffle(ids)
        for raw_id in ids:
            scanned += 1
            if max_proteins_scan > 0 and scanned > max_proteins_scan:
                break
            rec = idx[raw_id]
            seq = clean_aa(str(rec.seq))
            L = len(seq)
            if L < min_protein_len:
                continue
            b = bin_for_len(L)
            if not b or b not in pool:
                continue
            pid = raw_id.replace(" ", "_")[:120]
            pool[b].append((pid, seq))
    finally:
        idx.close()

    for b in pool:
        rng.shuffle(pool[b])

    remaining: dict[str, list[tuple[str, str]]] = {b: list(v) for b, v in pool.items()}
    chosen: list[tuple[str, str, str]] = []
    rows_summary: list[dict[str, object]] = []

    for b in sorted(bin_counts.keys(), key=lambda x: int(x.split("-")[0])):
        want = int(targets.get(b, 0))
        have = remaining[b]
        n_avail = len(have)
        take = min(want, n_avail)
        shortfall = want - take
        pick: list[tuple[str, str]] = []
        for _ in range(take):
            pick.append(have.pop())
        for i, (pid, seq) in enumerate(pick):
            oid = f"coding_pwhole|{b}|{pid}|{i}"
            chosen.append((oid, "CODING", seq))
        ref_n = int(bin_counts[b])
        ref_frac = ref_n / sum(bin_counts.values())
        rows_summary.append(
            {
                "bin": b,
                "bin_width_aa": bin_width,
                "ref_parent_count": ref_n,
                "ref_fraction": ref_frac,
                "target_sample_count": want,
                "proteome_candidates_in_bin": n_avail,
                "sampled_count_phase1": take,
                "shortfall_vs_target_phase1": shortfall,
            }
        )

    ref_frac_map = {b: bin_counts[b] / sum(bin_counts.values()) for b in bin_counts}
    bin_order = [b for b in sorted(bin_counts.keys(), key=lambda x: int(x.split("-")[0]))]
    if top_up:
        deficit = target_n - len(chosen)
        while deficit > 0:
            pools_left = [b for b in bin_order if remaining[b]]
            if not pools_left:
                break
            w = [ref_frac_map[b] for b in pools_left]
            s = float(sum(w))
            if s <= 0:
                pick_bin = rng.choice(pools_left)
            else:
                r = rng.random() * s
                acc = 0.0
                pick_bin = pools_left[-1]
                for b, wb in zip(pools_left, w):
                    acc += wb
                    if r <= acc:
                        pick_bin = b
                        break
            pid, seq = remaining[pick_bin].pop()
            oid = f"coding_pwhole|{pick_bin}|{pid}|top{target_n - deficit}"
            chosen.append((oid, "CODING", seq))
            deficit -= 1

    tot_bin: Counter[str] = Counter()
    for oid, _, _ in chosen:
        parts = oid.split("|")
        if len(parts) >= 2 and parts[0] == "coding_pwhole":
            tot_bin[parts[1]] += 1
    for row in rows_summary:
        b = str(row["bin"])
        total = int(tot_bin.get(b, 0))
        p1 = int(row["sampled_count_phase1"])
        row["sampled_count_total"] = total
        row["top_up_extra_vs_phase1"] = max(0, total - p1)

    if len(chosen) < target_n:
        print(
            f"Warning: proportional_whole filled {len(chosen)}/{target_n} parents "
            f"(proteome scan capped or sparse bins). See coding_proportional_bin_summary.csv.",
            file=sys.stderr,
        )

    summary = pd.DataFrame(rows_summary)
    tot_s = int(summary["sampled_count_total"].sum())
    if tot_s > 0:
        summary["sample_fraction_total"] = summary["sampled_count_total"].astype(float) / float(tot_s)
        summary["delta_sample_minus_ref_fraction"] = summary["sample_fraction_total"] - summary[
            "ref_fraction"
        ].astype(float)
    return chosen, summary


def write_proportional_vs_reference_report(
    bin_summary: pd.DataFrame,
    *,
    ref_parent_n: int,
    sampled_parent_n: int,
    peptide_label: str,
    out_md: Path,
) -> None:
    """Overall + per-bin QC: reference vs sampled bin fractions (Markdown)."""
    p = [float(x) for x in bin_summary["ref_fraction"]]
    q = [float(x) for x in bin_summary["sample_fraction_total"]]
    tvd = 0.5 * sum(abs(pi - qi) for pi, qi in zip(p, q))
    l2 = math.sqrt(sum((pi - qi) ** 2 for pi, qi in zip(p, q)))
    rc = [float(x) for x in bin_summary["ref_parent_count"]]
    st = [float(x) for x in bin_summary["sampled_count_total"]]
    nbin = len(rc)
    mean_r = sum(rc) / max(nbin, 1)
    mean_s = sum(st) / max(nbin, 1)
    num = sum((rc[i] - mean_r) * (st[i] - mean_s) for i in range(nbin))
    den_r = sum((x - mean_r) ** 2 for x in rc)
    den_s = sum((x - mean_s) ** 2 for x in st)
    if den_r > 0 and den_s > 0:
        r_bins = num / math.sqrt(den_r * den_s)
    else:
        r_bins = float("nan")

    lines = [
        "# Proportional whole-proteome sample vs reference length bins\n",
        f"- **Reference peptide set:** {peptide_label}\n",
        f"- **Reference parent MPs (with sequence):** {ref_parent_n}\n",
        f"- **Sampled whole proteome parents:** {sampled_parent_n}\n",
        "\n## Overall match (bin fractions)\n\n",
        "| Metric | Value |\n",
        "|--------|-------|\n",
        f"| Total variation distance (½Σ|p−q|) | {tvd:.4f} |\n",
        f"| L2 distance (‖p−q‖₂) | {l2:.4f} |\n",
        f"| Pearson r (ref bin count vs sampled bin count) | {r_bins:.4f} |\n",
        "\nPer-bin detail: `coding_proportional_bin_summary.csv` "
        "(`ref_fraction`, `sample_fraction_total`, `delta_sample_minus_ref_fraction`).\n",
    ]
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")


def write_fasta(path: Path, id_seq: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as w:
        for uid, pep in id_seq:
            w.write(f">{uid}\n{pep}\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="9-mer FASTA for NetMHCpan: significant lncRNA peptides vs length-matched coding."
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--peptide-faa",
        type=Path,
        default=SIG_FAA,
        help=(
            "Reference micropeptide FASTA whose length distribution drives controls (default: "
            "data/significant_lnc_peptides.faa ≈ 501 analyzed significant MPs). "
            "Do NOT pass data/smprot_tcga_filtered_peptides.faa unless you intend the full "
            "TCGA-matrix filtered set (~2.6k parents), not significant-only."
        ),
    )
    ap.add_argument(
        "--peptide-tsv",
        type=Path,
        default=SIG_TSV,
        help=(
            "TSV of smPEP_ID rows that define the reference cohort (default: data/significant_lnc_peptides.tsv). "
            "Must match the cohort implied by --peptide-faa."
        ),
    )
    ap.add_argument(
        "--coding-fa",
        type=Path,
        default=CODING_FA_DEFAULT,
        help="Human coding proteome FASTA (default: data/known_proteins.fasta, UniProt-style headers).",
    )
    ap.add_argument(
        "--max-proteins",
        type=int,
        default=0,
        help="Max proteome records to scan (0 = no limit; shuffle order varies with --seed).",
    )
    ap.add_argument(
        "--strict-fasta-match",
        action="store_true",
        help="Require every TSV smPEP_ID to appear in the FASTA (default: allow gaps, write a skip table).",
    )
    ap.add_argument(
        "--coding-control-mode",
        choices=("fragments", "whole_protein", "proportional_whole"),
        default="fragments",
        help=(
            "How to build coding controls from --coding-fa: "
            "'fragments' = random contiguous substring of each length from longer proteins (default); "
            "'whole_protein' = entire FASTA entry only if its length exactly matches (often infeasible for short MPs); "
            "'proportional_whole' = whole entries binned by length to match reference bin proportions (see --bin-width-aa, --proportional-target-n)."
        ),
    )
    ap.add_argument(
        "--proportional-target-n",
        type=int,
        default=None,
        help="For proportional_whole: number of proteome parents to draw (default: number of reference parents in FASTA).",
    )
    ap.add_argument(
        "--bin-width-aa",
        type=int,
        default=5,
        help="For proportional_whole: equal-width inclusive length bins aligned from reference min length.",
    )
    ap.add_argument(
        "--min-protein-len",
        type=int,
        default=9,
        help="For proportional_whole: ignore proteome entries shorter than this (no sliding 9-mers).",
    )
    ap.add_argument(
        "--no-proportional-top-up",
        action="store_true",
        help="For proportional_whole: skip the refill pass; sample at most min(target,available) per bin only.",
    )
    ap.add_argument(
        "--coding-ninemers-out",
        type=Path,
        default=None,
        help="For proportional_whole only: path for control 9-mer FASTA (default: "
        "data/netmhc/ninemers_coding_proportional_whole.fasta). Ignored for fragments/whole_protein.",
    )
    args = ap.parse_args()
    rng = random.Random(args.seed)

    if not args.peptide_faa.exists():
        print(
            "Missing peptide FASTA. Build it from the full significant peptide table, for example:\n"
            "  python export_tcga_filtered_peptides_fasta.py "
            "--peptides-tsv data/significant_lnc_peptides_full.tsv "
            "--transcripts-fa /path/to/gencode.*.transcripts.fa.gz\n"
            "Then: python sync_significant_lnc_peptides_with_fasta.py\n",
            file=sys.stderr,
        )
        raise SystemExit(1)

    if not args.peptide_tsv.exists():
        raise SystemExit(f"Missing {args.peptide_tsv}")

    if not args.coding_fa.exists():
        raise SystemExit(f"Missing coding proteome FASTA: {args.coding_fa}")

    OUT.mkdir(parents=True, exist_ok=True)

    expected_ids = load_expected_smpep_ids(args.peptide_tsv)
    tsv_df = pd.read_csv(args.peptide_tsv, sep="\t", dtype=str)
    sig_mps = parse_sig_peptide_fasta(args.peptide_faa)
    by_id = {sid: (sym, seq) for sid, sym, seq in sig_mps}

    missing = sorted(expected_ids - set(by_id))
    if missing:
        miss_df = tsv_df[tsv_df["smPEP_ID"].astype(str).str.strip().isin(missing)]
        miss_path = OUT / "sig_peptides_missing_from_fasta.tsv"
        miss_df.to_csv(miss_path, sep="\t", index=False)
        msg = (
            f"{args.peptide_faa} is missing {len(missing)} smPEP_ID row(s) from "
            f"{args.peptide_tsv.name} (listed in {miss_path})."
        )
        if args.strict_fasta_match:
            raise SystemExit(msg + " Re-run export or drop --strict-fasta-match.")
        print(f"Warning: {msg} Proceeding with {len(expected_ids) - len(missing)} FASTA records.", file=sys.stderr)

    extra = set(by_id) - expected_ids
    if extra:
        print(
            f"Note: ignoring {len(extra)} FASTA record(s) not listed in {args.peptide_tsv.name}.",
            file=sys.stderr,
        )

    tsv_ids = tsv_df["smPEP_ID"].astype(str).str.strip()
    ordered_ids: list[str] = []
    seen: set[str] = set()
    for sid in tsv_ids:
        if sid in seen:
            continue
        seen.add(sid)
        if sid in expected_ids:
            ordered_ids.append(sid)

    if len(ordered_ids) != len(expected_ids):
        raise SystemExit(
            f"TSV unique smPEP_ID count ({len(ordered_ids)}) != unique set size ({len(expected_ids)})."
        )

    ordered_ids_faa = [sid for sid in ordered_ids if sid in by_id]
    if not ordered_ids_faa:
        raise SystemExit("No peptide sequences overlap TSV and FASTA.")

    ordered: list[tuple[str, str, str]] = [
        (sid, by_id[sid][0], by_id[sid][1]) for sid in ordered_ids_faa
    ]

    lengths = Counter(len(s) for _, _, s in ordered)
    pd.DataFrame(
        [{"length_aa": L, "n_parent_mps": c} for L, c in sorted(lengths.items())]
    ).to_csv(OUT / "sig_mp_length_distribution.csv", index=False)

    pd.DataFrame(
        [{"smPEP_ID": a, "GeneSymbol": b, "length_aa": len(c), "sequence": c} for a, b, c in ordered]
    ).to_csv(OUT / "sig_parent_micropeptides.csv", index=False)

    if args.coding_control_mode == "whole_protein":
        coding_mps = sample_whole_protein_matched_coding_mps(
            lengths, args.coding_fa, rng, args.max_proteins
        )
    elif args.coding_control_mode == "proportional_whole":
        M = args.proportional_target_n if args.proportional_target_n is not None else len(ordered)
        coding_mps, bin_summary = sample_proportional_whole_proteins(
            lengths,
            args.coding_fa,
            rng,
            bin_width=max(1, int(args.bin_width_aa)),
            target_n=M,
            max_proteins_scan=args.max_proteins,
            min_protein_len=max(1, int(args.min_protein_len)),
            top_up=not args.no_proportional_top_up,
        )
        bin_summary.to_csv(OUT / "coding_proportional_bin_summary.csv", index=False)
    else:
        coding_mps = sample_length_matched_coding_mps(
            lengths, args.coding_fa, rng, args.max_proteins
        )
    (OUT / "coding_control_sampling_mode.txt").write_text(
        f"{args.coding_control_mode}"
        + (
            f"\nbin_width_aa={args.bin_width_aa}\nproportional_target_n={args.proportional_target_n or len(ordered)}\n"
            if args.coding_control_mode == "proportional_whole"
            else "\n"
        ),
        encoding="utf-8",
    )
    if args.coding_control_mode == "proportional_whole":
        parent_csv_path = OUT / "coding_proportional_whole_parent_micropeptides.csv"
    else:
        parent_csv_path = OUT / "coding_control_parent_micropeptides.csv"
    pd.DataFrame(
        [{"control_id": a, "label": b, "length_aa": len(c), "sequence": c} for a, b, c in coding_mps]
    ).to_csv(parent_csv_path, index=False)

    sig_9: list[tuple[str, str]] = []
    for sid, sym, seq in ordered:
        pid = f"{sid}|{sym}"
        sig_9.extend(sliding_ninemers(seq, pid, "SIG"))

    cod_9: list[tuple[str, str]] = []
    for oid, _, seq in coding_mps:
        cod_9.extend(sliding_ninemers(seq, oid, "COD"))

    write_fasta(OUT / "ninemers_sig_lnc.fasta", sig_9)
    if args.coding_control_mode == "proportional_whole":
        cod_fa_path = args.coding_ninemers_out or (OUT / "ninemers_coding_proportional_whole.fasta")
        write_fasta(cod_fa_path, cod_9)
        write_proportional_vs_reference_report(
            bin_summary,
            ref_parent_n=len(ordered),
            sampled_parent_n=len(coding_mps),
            peptide_label=str(args.peptide_faa),
            out_md=OUT / "coding_proportional_vs_reference_report.md",
        )
    else:
        cod_fa_path = OUT / "ninemers_coding_control.fasta"
        write_fasta(cod_fa_path, cod_9)

    meta = pd.DataFrame(
        [
            {
                "set": "significant_lncRNA_MP",
                "n_parent_mps": len(ordered),
                "n_ninemers": len(sig_9),
                "deduplicated_unique_ninemers": len({p for _, p in sig_9}),
                "ninemers_fasta": str(OUT / "ninemers_sig_lnc.fasta"),
            },
            {
                "set": (
                    "coding_proportional_whole_MP"
                    if args.coding_control_mode == "proportional_whole"
                    else "coding_control_MP"
                ),
                "n_parent_mps": len(coding_mps),
                "n_ninemers": len(cod_9),
                "deduplicated_unique_ninemers": len({p for _, p in cod_9}),
                "ninemers_fasta": str(cod_fa_path.resolve()),
            },
        ]
    )
    meta.to_csv(OUT / "ninemers_summary.csv", index=False)

    print(meta.to_string(index=False))
    print(f"Wrote {OUT / 'ninemers_sig_lnc.fasta'} ({len(sig_9)} records)")
    print(f"Wrote {cod_fa_path} ({len(cod_9)} records)")
    if args.coding_control_mode == "proportional_whole":
        print(f"Wrote {OUT / 'coding_proportional_vs_reference_report.md'}")
        print(f"Wrote {parent_csv_path.name} (proportional whole-protein parents)")
    print("Alleles: data/netmhc/hla_european27_class1.txt (edit to match your paper).")
    print("Next: run NetMHCpan-4.1 on Linux/WSL; see data/netmhc/README_netmhc.md")


if __name__ == "__main__":
    main()
