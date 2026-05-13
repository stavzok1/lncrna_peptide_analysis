"""
Histograms: significant lnc micropeptide lengths vs **proteome sampling** controls.

The “sampling from proteome” distribution is the empirical distribution of **parent
micropeptide lengths** in ``data/netmhc/coding_control_parent_micropeptides.csv`` —
the length-matched contiguous fragments drawn from the coding FASTA by
``prepare_netmhc_tr_vs_coding_epitopes.py`` (same design as NetMHC coding controls).

Significant lengths come from ``data/netmhc/sig_parent_micropeptides.csv``.

Optional: ``--full-proteome-protein-lengths`` adds a scan of the whole proteome FASTA
(entry lengths), which is **not** the sampling distribution used for controls.

Writes PNG + CSV under ``data/netmhc/figures/`` by default.
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

import figure_palettes as pal

DATA = ROOT / "data"
DEFAULT_PROTEOME = DATA / "known_proteins.fasta"
DEFAULT_SIG_CSV = DATA / "netmhc" / "sig_parent_micropeptides.csv"
DEFAULT_CONTROL_CSV = DATA / "netmhc" / "coding_control_parent_micropeptides.csv"
DEFAULT_OUT_DIR = DATA / "netmhc" / "figures"

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
AA20_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def clean_len(seq: str) -> int:
    s = "".join(c for c in str(seq).upper() if c in AA20)
    return len(s)


def load_lengths_csv(csv_path: Path) -> np.ndarray:
    df = pd.read_csv(csv_path)
    if "length_aa" in df.columns:
        return df["length_aa"].astype(int).to_numpy()
    if "sequence" not in df.columns:
        raise SystemExit(f"{csv_path} needs length_aa or sequence column")
    return np.array([clean_len(s) for s in df["sequence"].astype(str)], dtype=np.int64)


def load_sequences_csv(csv_path: Path) -> list[str]:
    df = pd.read_csv(csv_path)
    if "sequence" not in df.columns:
        raise SystemExit(f"{csv_path} needs a sequence column for AA composition")
    return [str(s) for s in df["sequence"].astype(str)]


def aa_frequencies(seqs: list[str]) -> np.ndarray:
    """Length-20 vector of AA frequencies (normalized to sum 1) over pooled sequences."""
    counts = np.zeros(20, dtype=np.float64)
    for s in seqs:
        for c in str(s).upper():
            if c in AA20:
                counts[AA20_ORDER.index(c)] += 1
    tot = counts.sum()
    if tot <= 0:
        return counts
    return counts / tot


def load_proteome_lengths(fa_path: Path, max_records: int) -> np.ndarray:
    lengths: list[int] = []
    for i, rec in enumerate(SeqIO.parse(fa_path, "fasta"), start=1):
        L = clean_len(str(rec.seq))
        if L > 0:
            lengths.append(L)
        if max_records > 0 and i >= max_records:
            break
    return np.asarray(lengths, dtype=np.int64)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Histogram: significant MP lengths vs proteome-sampled control lengths "
            "(coding_control_parent_micropeptides.csv)."
        )
    )
    ap.add_argument(
        "--sig-csv",
        type=Path,
        default=DEFAULT_SIG_CSV,
        help="sig_parent_micropeptides.csv",
    )
    ap.add_argument(
        "--control-csv",
        type=Path,
        default=DEFAULT_CONTROL_CSV,
        help="coding_control_parent_micropeptides.csv (sampled from --coding-fa in prep).",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Directory for PNG and CSV outputs.",
    )
    ap.add_argument(
        "--full-proteome-protein-lengths",
        action="store_true",
        help="Add a third panel: full proteome FASTA entry lengths (not the control sampler).",
    )
    ap.add_argument(
        "--proteome-fa",
        type=Path,
        default=DEFAULT_PROTEOME,
        help="Used only with --full-proteome-protein-lengths.",
    )
    ap.add_argument(
        "--max-proteome-records",
        type=int,
        default=0,
        help="Cap FASTA records when using --full-proteome-protein-lengths (0 = all).",
    )
    ap.add_argument(
        "--proteome-xmax",
        type=int,
        default=0,
        help="Clip full-proteome panel x-axis (0 = 99th percentile).",
    )
    args = ap.parse_args()

    if not args.sig_csv.exists():
        raise SystemExit(f"Missing {args.sig_csv}. Run: python prepare_netmhc_tr_vs_coding_epitopes.py")
    if not args.control_csv.exists():
        raise SystemExit(f"Missing {args.control_csv}. Run: python prepare_netmhc_tr_vs_coding_epitopes.py")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    sig_L = load_lengths_csv(args.sig_csv)
    samp_L = load_lengths_csv(args.control_csv)

    if sig_L.size == 0 or samp_L.size == 0:
        raise SystemExit("Empty length arrays from CSV.")

    lo = int(min(sig_L.min(), samp_L.min()))
    hi = int(max(sig_L.max(), samp_L.max()))
    bins = np.arange(lo, hi + 2) - 0.5

    sig_counts = (
        pd.Series(sig_L).value_counts().sort_index().rename_axis("length_aa").reset_index(name="n_sig")
    )
    ctrl_counts = (
        pd.Series(samp_L).value_counts().sort_index().rename_axis("length_aa").reset_index(name="n_proteome_sample")
    )
    sig_counts.merge(ctrl_counts, on="length_aa", how="outer").sort_values("length_aa").to_csv(
        args.out_dir / "sig_vs_proteome_sample_length_counts.csv", index=False
    )

    sig_seqs = load_sequences_csv(args.sig_csv)
    ctrl_seqs = load_sequences_csv(args.control_csv)
    f_sig = aa_frequencies(sig_seqs)
    f_ctrl = aa_frequencies(ctrl_seqs)
    pd.DataFrame(
        {"aa": AA20_ORDER, "freq_sig": f_sig, "freq_proteome_sample": f_ctrl, "delta_sig_minus_sample": f_sig - f_ctrl}
    ).to_csv(args.out_dir / "sig_vs_proteome_sample_aa_frequencies.csv", index=False)

    n_panels = 3 if args.full_proteome_protein_lengths else 2
    fig, axes = plt.subplots(n_panels, 1, figsize=(9, 3.2 * n_panels), constrained_layout=True)
    if n_panels == 2:
        ax_len, ax_aa = axes[0], axes[1]
    else:
        ax_len, ax_aa, ax_prot = axes[0], axes[1], axes[2]

    ax_len.hist(
        [sig_L, samp_L],
        bins=bins,
        label=["Significant lnc MPs", "Proteome-sampled controls"],
        color=[pal.SIG_LNC, pal.CODING_CONTROL],
        edgecolor="white",
        linewidth=0.5,
        alpha=0.85,
    )
    ax_len.set_xlabel("Parent micropeptide length (aa)")
    ax_len.set_ylabel("Count")
    ax_len.set_title(
        "Length: significant vs length-matched fragments sampled from coding proteome "
        f"(n_sig={len(sig_L):,}, n_sample={len(samp_L):,})"
    )
    ax_len.legend(loc="upper right", framealpha=0.9)
    ax_len.grid(True, alpha=0.25)

    x = np.arange(20)
    w = 0.35
    ax_aa.bar(x - w / 2, f_sig, width=w, label="Significant", color=pal.SIG_LNC, edgecolor="white", linewidth=0.3)
    ax_aa.bar(x + w / 2, f_ctrl, width=w, label="Proteome sample", color=pal.CODING_CONTROL, edgecolor="white", linewidth=0.3)
    ax_aa.set_xticks(x)
    ax_aa.set_xticklabels(AA20_ORDER, fontsize=9)
    ax_aa.set_ylabel("Pooled AA frequency")
    ax_aa.set_title("Amino-acid composition (pooled sequences): significant vs proteome-sampled controls")
    ax_aa.legend(loc="upper right", framealpha=0.9)
    ax_aa.grid(True, axis="y", alpha=0.25)

    if args.full_proteome_protein_lengths:
        if not args.proteome_fa.exists():
            raise SystemExit(f"Missing proteome FASTA: {args.proteome_fa}")
        prot_L = load_proteome_lengths(args.proteome_fa, args.max_proteome_records)
        if prot_L.size == 0:
            raise SystemExit("No proteome lengths read.")
        xmax_prot = args.proteome_xmax
        if xmax_prot <= 0:
            xmax_prot = int(np.percentile(prot_L, 99))
            xmax_prot = max(xmax_prot, hi + 1)
        prot_clip = prot_L[prot_L <= xmax_prot]
        bins1 = min(80, max(20, int(np.sqrt(max(len(prot_clip), 1)))))
        ax_prot.hist(prot_clip, bins=bins1, color=pal.LENGTH_HIST_ALT, edgecolor="white", linewidth=0.4, alpha=0.9)
        ax_prot.set_xlim(0, xmax_prot)
        ax_prot.set_xlabel("Full FASTA entry length (aa)")
        ax_prot.set_ylabel("Count (proteins)")
        ax_prot.set_title(
            f"Reference: whole proteome entries in {args.proteome_fa.name} "
            f"(n={len(prot_L):,}; plotted ≤{xmax_prot} aa, n={len(prot_clip):,})"
        )
        ax_prot.grid(True, alpha=0.25)
        edges = np.linspace(0, xmax_prot, min(101, xmax_prot + 1))
        if len(edges) < 2:
            edges = np.array([0, xmax_prot, xmax_prot + 1], dtype=float)
        hist_p, edges = np.histogram(prot_clip, bins=edges)
        pd.DataFrame(
            {"bin_left_aa": edges[:-1], "bin_right_aa": edges[1:], "n_proteins": hist_p}
        ).to_csv(args.out_dir / "full_proteome_length_histogram_bins.csv", index=False)

    fig.suptitle(
        "Significant lnc micropeptides vs proteome sampling (NetMHC prep controls)",
        fontsize=12,
        y=1.01,
    )
    out_png = args.out_dir / "proteome_sampling_vs_sig_mp_histograms.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    report = args.out_dir / "proteome_sampling_vs_sig_mp_report.txt"
    report.write_text(
        "\n".join(
            [
                f"sig_csv: {args.sig_csv}",
                f"control_csv (proteome sample): {args.control_csv}",
                f"n_sig: {len(sig_L)}",
                f"n_proteome_sample: {len(samp_L)}",
                f"sig_length_median: {float(np.median(sig_L)):.2f}",
                f"sample_length_median: {float(np.median(samp_L)):.2f}",
                f"figure: {out_png}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote {out_png}")
    print(f"Wrote {args.out_dir / 'sig_vs_proteome_sample_length_counts.csv'}")
    print(f"Wrote {args.out_dir / 'sig_vs_proteome_sample_aa_frequencies.csv'}")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
