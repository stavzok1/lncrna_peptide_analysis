"""
Length distribution of **whole** proteins used as the Fig. 5 coding control sample.

**Default:** reads ``length_aa`` from ``data/netmhc/coding_proportional_whole_parent_micropeptides.csv``
(the ~501 **whole** UniProt parents drawn with ``prepare_netmhc_tr_vs_coding_epitopes.py
--coding-control-mode proportional_whole``, matching ``netmhcpan_coding_proportional_whole.xls``).

Optional: ``--full-reference-proteome`` scans ``data/known_proteins.fasta`` (all records).

Writes PNG + binned CSV + summary TXT under ``figures/`` by default.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent
import figure_palettes as pal

DATA = ROOT / "data"
NETMHC = DATA / "netmhc"
DEFAULT_FIG5_PARENT_CSV = NETMHC / "coding_proportional_whole_parent_micropeptides.csv"
DEFAULT_FULL_PROTEOME_FA = DATA / "known_proteins.fasta"
FIGURES = ROOT / "figures"


def lengths_from_fig5_parent_csv(path: Path) -> tuple[np.ndarray, list[str]]:
    df = pd.read_csv(path)
    if "length_aa" not in df.columns:
        raise SystemExit(f"{path}: expected column length_aa")
    lengths = df["length_aa"].astype(int).to_numpy()
    if "control_id" in df.columns:
        ids = df["control_id"].astype(str).tolist()
    else:
        ids = [str(i) for i in range(len(df))]
    return lengths, ids


def lengths_from_fasta(path: Path, collect_ids: bool) -> tuple[np.ndarray, list[str] | None]:
    lengths: list[int] = []
    ids: list[str] | None = [] if collect_ids else None
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).replace(" ", "").replace("\n", "").upper()
        lengths.append(len(seq))
        if ids is not None:
            ids.append(rec.id)
    if not lengths:
        raise SystemExit(f"No FASTA records read from {path}")
    return np.asarray(lengths, dtype=np.int64), ids


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Protein length histogram: Fig 5 coding proportional-whole sample (default) or full proteome FASTA."
    )
    ap.add_argument(
        "--fig5-coding-parent-csv",
        type=Path,
        default=DEFAULT_FIG5_PARENT_CSV,
        help="Fig 5 coding control: parent micropeptide table with length_aa (default: proportional_whole cohort).",
    )
    ap.add_argument(
        "--full-reference-proteome",
        action="store_true",
        help="Instead of the Fig 5 CSV, scan the full reference proteome FASTA (--proteome-fa).",
    )
    ap.add_argument(
        "--proteome-fa",
        type=Path,
        default=DEFAULT_FULL_PROTEOME_FA,
        help="Used only with --full-reference-proteome (default: data/known_proteins.fasta).",
    )
    ap.add_argument("--out-dir", type=Path, default=FIGURES, help="Directory for PNG + CSV + summary TXT.")
    ap.add_argument(
        "--out-stem",
        type=str,
        default="fig5_coding_proportional_whole_parent_length_distribution",
        help="Output filename stem.",
    )
    ap.add_argument("--bin-width", type=int, default=5, help="Histogram bin width in amino acids (default 5).")
    ap.add_argument(
        "--max-length-cap",
        type=int,
        default=None,
        metavar="L",
        help="Optional: extend histogram grid; default covers min..max of the data.",
    )
    ap.add_argument(
        "--per-protein-csv",
        type=Path,
        default=None,
        help="Optional path to write every parent id + length.",
    )
    args = ap.parse_args()

    if args.full_reference_proteome:
        if not args.proteome_fa.is_file():
            raise SystemExit(f"Missing proteome FASTA: {args.proteome_fa}")
        arr, ids = lengths_from_fasta(args.proteome_fa, collect_ids=args.per_protein_csv is not None)
        title_suffix = args.proteome_fa.name
    else:
        if not args.fig5_coding_parent_csv.is_file():
            raise SystemExit(
                f"Missing {args.fig5_coding_parent_csv}\n"
                "Build proportional-whole coding inputs with:\n"
                "  python prepare_netmhc_tr_vs_coding_epitopes.py --coding-control-mode proportional_whole --max-proteins 0\n"
                "Or pass --full-reference-proteome to use the full proteome FASTA instead."
            )
        arr, ids = lengths_from_fig5_parent_csv(args.fig5_coding_parent_csv)
        source_line = f"fig5_coding_parent_csv: {args.fig5_coding_parent_csv.resolve()}"
        title_suffix = "Fig 5 coding (proportional whole proteins)"

    n = int(arr.size)
    mean = float(arr.mean())
    sd = float(arr.std(ddof=0))
    qs = np.percentile(arr, [5, 25, 50, 75, 95, 99]).astype(float)
    mn, mx = int(arr.min()), int(arr.max())

    bw_user = max(int(args.bin_width), 1)
    max_bins = 200
    bw = bw_user
    if mx > 0 and int(np.ceil(mx / bw)) > max_bins:
        bw = int(np.ceil(mx / max_bins))

    cap = args.max_length_cap
    if cap is None:
        cap = int(np.ceil(mx / bw) * bw)
    else:
        cap = int(cap)
    edges = np.arange(0, cap + bw, bw, dtype=np.int64)
    if edges[-1] < mx:
        edges = np.append(edges, int(np.ceil(mx / bw) * bw) + bw)
    hist, bin_edges = np.histogram(arr, bins=edges)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    df_hist = pd.DataFrame(
        {
            "bin_low_aa": bin_edges[:-1].astype(int),
            "bin_high_exclusive_aa": bin_edges[1:].astype(int),
            "n_proteins": hist.astype(int),
        }
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.out_stem
    csv_path = args.out_dir / f"{stem}_histogram.csv"
    txt_path = args.out_dir / f"{stem}_summary.txt"
    png_path = args.out_dir / f"{stem}.png"

    df_hist.to_csv(csv_path, index=False)

    lines = [
        source_line,
        f"n_proteins: {n:,}",
        f"length_aa: min={mn:,} max={mx:,} mean={mean:.2f} sd={sd:.2f}",
        f"quantiles (aa): p5={qs[0]:.0f} p25={qs[1]:.0f} p50={qs[2]:.0f} p75={qs[3]:.0f} p95={qs[4]:.0f} p99={qs[5]:.0f}",
        f"histogram_bin_width_aa: {bw}",
        f"histogram_csv: {csv_path}",
        f"figure_png: {png_path}",
    ]
    if args.per_protein_csv is not None:
        pd.DataFrame({"parent_id": ids, "length_aa": arr}).to_csv(args.per_protein_csv, index=False)
        lines.append(f"per_protein_csv: {args.per_protein_csv.resolve()}")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    ax.bar(centers, hist, width=float(bw) * 0.92, align="center", color=pal.LENGTH_HIST_SINGLE, edgecolor="white", linewidth=0.35)
    ax.set_xlabel("Protein length (amino acids)")
    ax.set_ylabel("Number of proteins")
    ax.set_title(f"Whole-protein length distribution ({title_suffix}, n={n:,})")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.grid(True, axis="y", linestyle=":", alpha=0.45)
    fig.tight_layout()
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)

    print("\n".join(lines))


if __name__ == "__main__":
    main()
