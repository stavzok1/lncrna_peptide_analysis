"""
Scatter: allele population frequency (x) vs **SB epitopes** or unique 9-mers (y) for
**strong binders** by predicted IC50 < threshold (default 150 nM). Default **y** is **SB epitopes per allele**
(prediction rows per allele with IC50 < cutoff; each row is one 9-mer × allele); use ``--y-metric unique`` for distinct 9-mers per allele.

NetMHCpan-4.2 wide XLS (``-xls 1 -BA 1``) does **not** print a literal ``IC50_nM`` column; it
prints ``BA_score`` in (0, 1) and ``BA_rank``. Those come from the **same binding-affinity
(BA) head** as the IC50-ranked predictions; the BA score is the usual **log-scale affinity
compressed to [0, 1] with a 50 000 nM ceiling** used in the NetMHCpan / MHCflurry family
(NetMHCpan-4.0/4.x binding models, Jurtz *et al.*, *J Immunol* 2017). The **inverse** used
everywhere in tooling that reconstructs nM from that score is::

    IC50_nM = 50000 ** (1 - BA_score)

equivalently ``10 ** ((1 - BA_score) * log10(50000))``. This is **not new physics**—it is the
same information as ``BA_score``, rewritten into nM so you can apply an IC50 cutoff (e.g.
150 nM). You **do not** need to rerun NetMHCpan solely to obtain nM from an existing wide
XLS unless you want a different export (e.g. per-row text output where the binary prints
``Affinity`` explicitly). **NetMHCpan’s built-in SB/WB labels** in EL mode use **%Rank**
defaults, not IC50 < 150 nM—that threshold is your analysis rule on the reconstructed IC50.

Inputs
~~~~~~
- ``--netmhc-xls``: tab-separated NetMHCpan wide output (first line ``#...netMHCpan...``).
- ``--freq-file``: text file with allele + frequency (tab or comma). Recognized columns:
  ``allele`` / ``HLA`` / first column, and ``freq`` / ``frequency`` / second column.
  Alleles are matched to XLS alleles after normalization (``HLA-A*02:01`` vs ``HLA-A02:01``).

Writes PNG + CSV summary under ``data/netmhc/figures/`` by default (default filenames use a
``fig5a_`` prefix). The same files are copied to repo-root ``figures/``.
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy import stats

ROOT = Path(__file__).resolve().parent
import figure_palettes as pal

NET = ROOT / "data" / "netmhc"
FIGS = NET / "figures"
REPO_FIGURES = ROOT / "figures"
_DEFAULT_ALLELE_FREQ_CSV = FIGS / "fig5a_epitopes_vs_allele_frequency_ic50_sb.csv"
_ALT_ALLELE_FREQ_CSV = FIGS / "epitopes_vs_allele_frequency_ic50_sb.csv"


def _mirror_to_repo_figures(*paths: Path) -> None:
    """Copy PNG/CSV outputs into repo-root ``figures/`` alongside other manuscript figures."""
    REPO_FIGURES.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.is_file():
            shutil.copy2(p, REPO_FIGURES / p.name)


def ba_score_to_ic50_nm(ba: float) -> float:
    """Map NetMHCpan wide-format BA_score to predicted IC50 (nM)."""
    x = float(ba)
    if not math.isfinite(x):
        return float("nan")
    x = min(1.0 - 1e-12, max(1e-12, x))
    return float(50000.0 ** (1.0 - x))


def normalize_allele(s: str) -> str:
    t = str(s).strip().upper()
    t = t.replace("*", "")
    t = re.sub(r"\s+", "", t)
    if not t.startswith("HLA-"):
        t = "HLA-" + t.lstrip("-")
    return t


def alleles_from_xls_header_line1(line: str) -> list[str]:
    m = re.search(r"-a\s+([^\s]+(?:,[^\s]+)*)", line)
    if not m:
        raise ValueError("Could not parse -a allele list from first line of XLS")
    raw = m.group(1).strip()
    return [normalize_allele(a) for a in raw.split(",") if a.strip()]


def parse_wide_netmhc_xls(path: Path) -> tuple[list[str], np.ndarray, list[str]]:
    """
    Return (allele_names, ba_matrix, peptides) where ba_matrix has shape (n_rows, n_alleles)
    with BA_score per allele (NaN if missing), and peptides are the 9-mer strings per row
    (standard AA only; same filter as downstream IC50 logic).
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError(f"Too few lines in {path}")
    alleles = alleles_from_xls_header_line1(lines[0])
    n_per = 10  # core,icore,EL_score,EL_rank,BA_score,BA_rank,Pathogen_*,Neo_*
    ba_col_idx = [3 + i * n_per + 4 for i in range(len(alleles))]

    rows: list[list[float]] = []
    peps: list[str] = []
    for ln in lines[3:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 3 + n_per:
            continue
        pep = parts[1].strip().upper()
        if len(pep) != 9 or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in pep):
            continue
        peps.append(pep)
        vals: list[float] = []
        for j in ba_col_idx:
            if j >= len(parts):
                vals.append(float("nan"))
                continue
            try:
                vals.append(float(parts[j]))
            except ValueError:
                vals.append(float("nan"))
        rows.append(vals)
    return alleles, np.asarray(rows, dtype=np.float64), peps


def display_allele(norm: str) -> str:
    """HLA-A02:01 -> HLA-A*02:01 for plot labels."""
    m = re.match(r"HLA-([AB])(\d{2}:\d{2})$", norm)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return norm


def load_freq_table(path: Path) -> dict[str, float]:
    df = pd.read_csv(path, sep=None, engine="python", header=0, comment="#")
    cols = [str(c).strip().lower() for c in df.columns]
    df.columns = cols
    if "allele" in cols and ("frequency" in cols or "freq" in cols):
        a_col, f_col = "allele", "frequency" if "frequency" in cols else "freq"
    elif len(cols) >= 2:
        a_col, f_col = cols[0], cols[1]
    else:
        raise SystemExit(f"Need at least two columns in {path}; got {cols}")
    out: dict[str, float] = {}
    for _, r in df.iterrows():
        key = normalize_allele(str(r[a_col]))
        out[key] = float(str(r[f_col]).replace(",", "."))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Epitope SB count vs allele frequency (IC50 from BA_score).")
    ap.add_argument("--netmhc-xls", type=Path, default=NET / "netmhcpan_sig_lnc.xls")
    ap.add_argument(
        "--freq-file",
        type=Path,
        default=None,
        help=(
            "TSV/CSV: allele, frequency. Default: data/netmhc/figures/fig5a_epitopes_vs_allele_frequency_ic50_sb.csv "
            "(or epitopes_vs_allele_frequency_ic50_sb.csv) when present."
        ),
    )
    ap.add_argument("--ic50-nm", type=float, default=150.0, help="Strong binder if predicted IC50 < this (nM).")
    ap.add_argument(
        "--y-metric",
        choices=("instances", "unique"),
        default="instances",
        help="Per-allele y-axis: SB prediction rows vs unique 9-mers (default instances).",
    )
    ap.add_argument(
        "--out-png",
        type=Path,
        default=FIGS / "fig5a_epitopes_vs_allele_frequency_ic50_sb.png",
    )
    ap.add_argument(
        "--out-csv",
        type=Path,
        default=FIGS / "fig5a_epitopes_vs_allele_frequency_ic50_sb.csv",
    )
    ap.add_argument(
        "--no-repo-mirror",
        action="store_true",
        help="Do not copy PNG/CSV into repo-root figures/ (use when writing under figures/... subfolders).",
    )
    args = ap.parse_args()

    freq_path = args.freq_file
    if freq_path is None:
        for candidate in (_DEFAULT_ALLELE_FREQ_CSV, _ALT_ALLELE_FREQ_CSV):
            if candidate.is_file():
                freq_path = candidate
                break
    if freq_path is None or not freq_path.is_file():
        raise SystemExit(
            "Missing allele frequency table. Pass --freq-file, or place one of:\n"
            f"  {_DEFAULT_ALLELE_FREQ_CSV}\n"
            f"  {_ALT_ALLELE_FREQ_CSV}"
        )

    alleles, ba_mat, peps = parse_wide_netmhc_xls(args.netmhc_xls)
    ic50 = np.vectorize(ba_score_to_ic50_nm)(ba_mat)
    strong = ic50 < args.ic50_nm

    freqs = load_freq_table(freq_path)
    xs: list[float] = []
    y_uniq: list[int] = []
    y_inst: list[int] = []
    labels: list[str] = []
    union_sb: set[str] = set()

    for i, al in enumerate(alleles):
        f = freqs.get(al)
        if f is None:
            alt = al.replace("HLA-", "")
            f = freqs.get(normalize_allele(alt))
        if f is None:
            print(f"Warning: no frequency for {al}, skipping", flush=True)
            continue
        mask = strong[:, i]
        uniq = {peps[j] for j in range(len(peps)) if mask[j]}
        n_u = len(uniq)
        n_i = int(mask.sum())
        union_sb |= uniq
        xs.append(float(f))
        y_uniq.append(n_u)
        y_inst.append(n_i)
        labels.append(display_allele(al))

    xs_arr = np.asarray(xs, dtype=np.float64)
    y_u_arr = np.asarray(y_uniq, dtype=np.float64)
    y_i_arr = np.asarray(y_inst, dtype=np.float64)
    ys_arr = y_i_arr if args.y_metric == "instances" else y_u_arr
    rho, pval = stats.spearmanr(xs_arr, ys_arr)
    sum_per_allele_unique = int(y_u_arr.sum())
    sum_per_allele_instances = int(y_i_arr.sum())

    FIGS.mkdir(parents=True, exist_ok=True)
    df_out = pd.DataFrame(
        {
            "allele": labels,
            "allele_frequency": xs_arr,
            f"n_sb_row_instances_ic50_lt_{int(args.ic50_nm)}nm": y_i_arr.astype(int),
            f"n_unique_epitopes_ic50_lt_{int(args.ic50_nm)}nm": y_u_arr.astype(int),
        }
    )
    df_out.to_csv(args.out_csv, index=False)

    fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=150)
    ax.scatter(xs_arr, ys_arr, s=36, c=pal.SCATTER_DEFAULT, alpha=0.85, edgecolors="none")
    for x, y, lab in zip(xs_arr, ys_arr, labels):
        ax.annotate(lab.replace("HLA-", ""), (x, y), fontsize=6, xytext=(3, 2), textcoords="offset points")

    ax.set_xlabel("Allele frequency")
    ax.set_ylabel(
        "SB epitopes per allele" if args.y_metric == "instances" else "Unique 9-mers per allele"
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    xls_name = args.netmhc_xls.name.lower()
    if "coding" in xls_name:
        sample_line = "Proteome control sample"
    else:
        sample_line = "Tr-lncRNA MPs"
    y_note = (
        f"Σ SB epitopes × alleles={sum_per_allele_instances}"
        if args.y_metric == "instances"
        else f"Σ per-allele unique 9-mers={sum_per_allele_unique}"
    )
    ax.set_title(
        f"Epitope production vs allele frequency ({sample_line})\n"
        f"Spearman ρ={rho:.3f}, p={pval:.2g}, {y_note}, "
        f"n_union={len(union_sb)}, n_alleles={len(xs_arr)}"
    )
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)
    if not args.no_repo_mirror:
        _mirror_to_repo_figures(args.out_png, args.out_csv)
    print(f"Wrote {args.out_png}\n{args.out_csv}")
    print(f"Spearman rho={rho:.4f} p={pval:.4g} per-allele points={len(xs_arr)}")


if __name__ == "__main__":
    main()
