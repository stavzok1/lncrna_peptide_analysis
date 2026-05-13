"""
CD-HIT clustering (WSL) of TCGA-filtered peptide sequences per gene, at 60% and
90% global sequence identity (protein), then plots per gene.

CD-HIT 4.8.1 is appropriate for UniRef-like greedy clustering; identity is controlled
with -c (1.0 = identical). Word length -n and related flags follow the official CD-HIT
protein guide; see docs/analysis_params.md.

Inputs: data/smprot_tcga_filtered_peptides.faa (headers >smPEP|GeneSymbol|...)

Outputs under data/cdhit_clustering/<GENE>/:
  - gene.fa, cdhit_0.6.clstr, cdhit_0.9.clstr, summary.csv
  - <GENE>_cdhit_clustering.png
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
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure_palettes as pal

DATA = ROOT / "data"
FAA = DATA / "smprot_tcga_filtered_peptides.faa"
OUT_ROOT = DATA / "cdhit_clustering"


def win_path_to_wsl(p: Path) -> str:
    s = str(p.resolve())
    if len(s) >= 2 and s[1] == ":":
        drive = s[0].lower()
        rest = s[2:].replace("\\", "/")
        return f"/mnt/{drive}{rest}"
    raise ValueError(f"Need Windows path with drive letter, got {s}")


def parse_faa_by_gene(faa_path: Path, gene: str) -> list[tuple[str, str]]:
    """Return list of (seq_id, aa_sequence) for records where GeneSymbol == gene."""
    out: list[tuple[str, str]] = []
    cur_id: str | None = None
    with faa_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                parts = line[1:].split("|")
                cur_id = parts[0] if parts else line[1:]
                sym = parts[1] if len(parts) > 1 else ""
                if sym != gene:
                    cur_id = None
            elif cur_id is not None:
                out.append((cur_id, line.upper()))
                cur_id = None
    return out


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as w:
        for sid, seq in records:
            w.write(f">{sid}\n{seq}\n")


def cdhit_params(identity: float) -> tuple[int, int]:
    """
    Return ``(-n`` word length, ``-l`` min length) for protein ``cd-hit``.

    ``-n`` uses the official CD-HIT protein bands (weizhongli/cdhit wiki,
    "Choice of word size" for ``cd-hit``): 5 for 0.7–1.0, 4 for 0.6–0.7,
    3 for 0.5–0.6, 2 for 0.4–0.5. ``-l`` is set to cd-hit's default (10);
    input peptides are already >=10 aa upstream.
    """
    if identity >= 0.7:
        return 5, 10
    if identity >= 0.6:
        return 4, 10
    if identity >= 0.5:
        return 3, 10
    return 2, 10


def run_cdhit_wsl(in_fa: Path, out_prefix: Path, identity: float) -> Path:
    """Run cd-hit via WSL; return path to .clstr file."""
    wsl_in = win_path_to_wsl(in_fa)
    wsl_out = win_path_to_wsl(out_prefix)
    n, min_len = cdhit_params(identity)
    clstr = Path(str(out_prefix) + ".clstr")
    if clstr.exists():
        clstr.unlink()
    cmd = [
        "wsl",
        "cd-hit",
        "-i",
        wsl_in,
        "-o",
        wsl_out,
        "-c",
        str(identity),
        "-n",
        str(n),
        "-d",
        "0",
        "-T",
        "0",
        "-l",
        str(min_len),
        "-M",
        "16000",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"cd-hit failed ({r.returncode}): {r.stderr or r.stdout}\nCMD: {' '.join(cmd)}"
        )
    return clstr


def parse_clstr(clstr_path: Path) -> list[int]:
    """Return list of cluster sizes (number of member sequences per cluster)."""
    sizes: list[int] = []
    cur = 0
    with clstr_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                if cur > 0:
                    sizes.append(cur)
                cur = 0
            elif re.match(r"^\d+", line):
                cur += 1
        if cur > 0:
            sizes.append(cur)
    return sizes


def plot_gene_figure(
    gene: str,
    sizes_60: list[int],
    sizes_90: list[int],
    out_png: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    def one_panel(ax, sizes: list[int], title: str) -> None:
        if not sizes:
            ax.text(0.5, 0.5, "No clusters", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            return
        s = np.array(sorted(sizes, reverse=True), dtype=float)
        q = np.quantile(s, [0.25, 0.5, 0.75])
        x = np.arange(1, len(s) + 1)
        ax.bar(x, s, width=0.85, color=pal.OI_SKY_BLUE, edgecolor="none")
        for qi, val in zip(("Q1", "Median", "Q3"), q):
            ax.axhline(val, linestyle="--", linewidth=1.0, alpha=0.85, label=f"{qi}={val:.1f}")
        cum = np.cumsum(s)
        tot = cum[-1]
        ax2 = ax.twinx()
        ax2.plot(x, cum / tot, color=pal.OI_ORANGE, linewidth=2.0, alpha=0.9, label="Cumulative peptide fraction")
        ax2.set_ylabel("Cumulative fraction of peptides", color=pal.OI_ORANGE)
        ax2.set_ylim(0, 1.05)
        ax2.tick_params(axis="y", labelcolor=pal.OI_ORANGE)
        for thr in (0.25, 0.5, 0.75):
            idx = int(np.searchsorted(cum / tot, thr, side="left"))
            if idx < len(x):
                ax2.axvline(x[idx], color="gray", linestyle=":", alpha=0.45)
        ax.set_xlabel("Cluster rank (by size, largest first)")
        ax.set_ylabel("Peptides in cluster")
        ax.set_title(f"{title}\n(n_clusters={len(s)}, n_peptides={int(tot)})")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    one_panel(axes[0], sizes_60, f"{gene}  CD-HIT  ~60% identity (-c 0.6)")
    one_panel(axes[1], sizes_90, f"{gene}  CD-HIT  ~90% identity (-c 0.9)")
    fig.suptitle(f"{gene}: TCGA-filtered peptide clustering (CD-HIT 4.8.1)", fontsize=12)
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="+", default=["TTN-AS1", "PVT1"])
    ap.add_argument("--faa", type=Path, default=FAA)
    args = ap.parse_args()

    if not args.faa.exists():
        raise SystemExit(f"Missing {args.faa}; run export_tcga_filtered_peptides_fasta.py first.")

    for gene in args.genes:
        recs = parse_faa_by_gene(args.faa, gene)
        if not recs:
            print(f"No peptides for {gene} in {args.faa.name}; skip.")
            continue
        gdir = OUT_ROOT / gene.replace(" ", "_")
        gdir.mkdir(parents=True, exist_ok=True)
        in_fa = gdir / f"{gene.replace('-', '_')}_peptides.fa"
        write_fasta(recs, in_fa)
        print(f"{gene}: {len(recs)} sequences -> {in_fa}")

        rows = []
        fig_sizes: dict[str, list[int]] = {}
        for ident, tag in ((0.6, "0.6"), (0.9, "0.9")):
            prefix = gdir / f"cdhit_{tag}"
            clstr = run_cdhit_wsl(in_fa, prefix, ident)
            sizes = parse_clstr(clstr)
            fig_sizes[tag] = sizes
            if sum(sizes) != len(recs):
                print(
                    f"  note: -c {tag} cluster member sum={sum(sizes)} vs input={len(recs)} "
                    "(CD-HIT can omit duplicate sequences or very short outliers)."
                )
            rows.append(
                {
                    "gene": gene,
                    "identity_cutoff": tag,
                    "n_peptides_input": len(recs),
                    "n_clusters": len(sizes),
                    "largest_cluster": max(sizes) if sizes else 0,
                    "singleton_clusters": sum(1 for x in sizes if x == 1),
                }
            )
            print(f"  -c {tag}: clusters={len(sizes)}, sizes (first 10 desc)={sorted(sizes, reverse=True)[:10]}")

        pd.DataFrame(rows).to_csv(gdir / "cdhit_summary.csv", index=False)
        plot_gene_figure(gene, fig_sizes["0.6"], fig_sizes["0.9"], gdir / f"{gene.replace('-', '_')}_cdhit_clustering.png")
        print(f"Wrote figure -> {gdir / (gene.replace('-', '_') + '_cdhit_clustering.png')}")


if __name__ == "__main__":
    main()
