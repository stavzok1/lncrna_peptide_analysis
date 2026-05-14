"""
**Manuscript / catalog Figure 6** — TTN-AS1 (smPEP_ID 108065): allele coverage along the
parent peptide, histogram, along-sequence epitope vs allele tracks, and sequence logos.

This is the **single-peptide** NetMHCpan story for TTN-AS1. For **cohort-level** allele vs
epitope interplay (significant lncRNA MPs vs coding control, optionally with IEDB filters),
use the **Figure 5** NetMHC scripts documented in ``docs/figure_catalog.md``.

**Sensitivity / grids (supplement) for the same XLS and sequence:**

- NetMHC-only SB sweeps: ``supplement/plot_figure6_ttn_as1_sb_sensitivity.py`` (default under ``data/netmhc/figures/fig6_ttn_as1_sensitivity/``; orchestrator mirrors under ``figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_wide_netmhc_sb_sweeps/``).
- IEDB+NetMHC **1D + LOO** on the merged TTN table: ``supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py`` (orchestrator: ``fig6_ttn_merged_iedb_1d_sensitivity_loo/``).
- IEDB+NetMHC **Cartesian** SB threshold grid on the merged TTN table: ``supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py``.

**Gating modes (``--gating``):**

- ``iedb_sb`` (**default**): NetMHCpan wide rows joined to an **IEDB peptide_table CSV** on
  ``stable_key`` (see ``--iedb-csv`` + ``--iedb-parent-input-seq-id``). SB uses **one**
  threshold bundle from the CLI (``--sb-mode`` + imm/proc/EL/IC50); **IC50 uses the merged IEDB
  IC50 column when available** (else BA-derived fallback in ``sb_mask_spec``). Default IC50 cap
  matches ``FIG5_IEDB_IC50_MAX_NM_DEFAULT`` (**150 nM**). This is **not** a multi-threshold figure by
  itself—use the **supplement scripts above** for NetMHC-only sweeps, merged IEDB 1D+LOO, and merged IEDB Cartesian grids.
- ``netmhc``: SB from **wide XLS only** — **BA_rank** and/or **IC50 from BA_score**
  (plus optional **EL_rank**). No IEDB immunogenicity / processing.

**Coverage outputs (``--coverage-output``):** default **``instances``** — **panels A–B**
(sequence heatmap, histogram) use **distinct SB alleles per site**; **panel C** line plots
use **SB epitope instance hits** (track 1), **distinct SB alleles** (track 2), and their
**overlay** (track 3). **Panels D/E logos** remain **instance-weighted**. Optional **``unique``**
or **``both``**; ``both`` writes **separate files** ``*_instances*`` and ``*_unique*`` (see
``-o`` stem). Use ``--also-write-unique`` when ``--coverage-output`` is ``instances`` to add
the unique-only companion files (same A–C metrics as ``instances`` for this figure). **Repo
orchestrators** default to **instances only** under ``figures/``; Fig 6 **unique** is opt-in
(``--write-fig6-unique-supplement`` on ``generate_canonical_manuscript_figures.py``,
``--include-fig6-unique-split`` on ``export_publication_figures.py``, ``--include-fig6-unique`` on
``regenerate_all_figures.py``, or ``--also-write-unique`` on ``generate_netmhc_figure_bundle.py``).

**Panel B (histogram + stats):** with **``instances``**, the histogram counts positions by
**distinct SB alleles** per site (union across overlapping 9-mers). With **``unique``**, same
allele-based histogram and stats.

**Panels C (``instances`` / ``unique``):** line tracks are **(A)** SB epitope **instance** load
(peptide×allele hit sum per site), **(B)** **distinct** SB alleles per site, **(C)** overlay of
those two curves. (Heatmap / histogram still use distinct alleles per site.)

**Panels D/E (logos):** **instance-weighted** — each SB peptide×allele hit adds one copy of the
9-mer to the frequency matrix (not deduplicated by sequence). Titles report ``N`` = total
weighted rows (and the subtitle still reports how many **distinct** 9-mer sequences those hits
come from).

Strong binders (configurable, default **BA_rank**):

- ``ba_rank``: NetMHCpan wide-XLS **BA_rank** ≤ cutoff (default **0.5**), i.e. top 0.5 %
  predicted binders for that allele (rank is a 0–100 %-tile style score; lower = stronger).
- ``--require-el-rank``: additionally require **EL_rank** ≤ the same cutoff by default
  (``--el-rank-pct`` overrides EL only). EL = elution / ligand model head in the same export.
- ``ic50``: ``IC50_nM = 50000 ** (1 - BA_score)`` and SB if IC50 < threshold nM (same as
  ``plot_netmhc_epitopes_vs_hla_frequency.py``); ``--require-el-rank`` can be combined.

``--split-panels``: write **five** PNGs next to ``-o`` (must end with ``.png``):
``{stem}_a.png`` … ``{stem}_e.png`` (sequence, histogram, three-track coverage, two logos).
Default ``-o`` is under repo-root ``figures/`` (see ``--out``).

Requires: matplotlib, numpy, pandas. Panels **D/E** use **sequence logos** when
``logomaker`` is importable and ``logomaker.Logo`` succeeds; otherwise they fall back to
**stacked per-position bars** (same frequency matrix). Install: ``pip install logomaker``.
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
_MS = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_REPO / "scripts"), str(_MS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import REPO_ROOT, DATA, FIGURES, NETMHC_DATA, NETMHC_FIGURES
from figure_export import add_publication_args, save_figure_bundle

ROOT = REPO_ROOT


import argparse
import math
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, to_rgb
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator

import figure_palettes as pal

NET = ROOT / "data" / "netmhc"
REPO_FIGURES = ROOT / "figures"
TTN_IEDB_SYNTHETIC_CSV = NET / "ttn_as1_iedb_companion_synthetic.csv"
TTN_IEDB_PARENT_SEQ_ID_DEFAULT = "108065|TTN-AS1|synthetic"
_SCRIPTS = ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from netmhc_sb_core import (  # noqa: E402
    FIG5_IEDB_EL_RANK_MAX_DEFAULT,
    FIG5_IEDB_IC50_MAX_NM_DEFAULT,
    FIG5_IEDB_IMM_MIN_DEFAULT,
    FIG5_IEDB_PROC_MIN_DEFAULT,
    ba_score_min_for_ic50_lt,
    pick_iedb_ic50_column,
    sb_mask_spec,
    sb_spec_from_mode,
)

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AAS_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def ba_score_to_ic50_nm(ba: float) -> float:
    x = float(ba)
    if not math.isfinite(x):
        return float("nan")
    x = min(1.0 - 1e-12, max(1e-12, x))
    return float(50000.0 ** (1.0 - x))


def normalize_allele(s: str) -> str:
    t = str(s).strip().upper().replace("*", "")
    t = re.sub(r"\s+", "", t)
    if not t.startswith("HLA-"):
        t = "HLA-" + t.lstrip("-")
    return t


def display_allele(norm: str) -> str:
    m = re.match(r"HLA-([AB])(\d{2}:\d{2})$", norm)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return norm


def normalize_hla_netmhc_to_iedb(name: str) -> str:
    """HLA-A01:01 -> HLA-A*01:01 (IEDB stable_key allele token)."""
    s = name.strip()
    if not s:
        return s
    up = s.upper()
    m = re.fullmatch(r"HLA-([ABCEFG])(\d{2}:\d{2})", up)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return s


def alleles_from_xls_header_line1(line: str) -> list[str]:
    m = re.search(r"-a\s+([^\s]+(?:,[^\s]+)*)", line)
    if not m:
        raise ValueError("Could not parse -a allele list from first line of XLS")
    raw = m.group(1).strip()
    return [normalize_allele(a) for a in raw.split(",") if a.strip()]


def parse_wide_netmhc_xls_rows(
    path: Path,
) -> tuple[list[str], list[int], list[str], np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (alleles, starts, peptides, ba_matrix, ba_rank_matrix, el_rank_matrix).

    Per-allele block (10 columns): core, icore, EL_score, EL_rank, BA_score, BA_rank, ...
    Rows use only standard-AA peptides (drops NetMHC X-padding rows).
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError(f"Too few lines in {path}")
    alleles = alleles_from_xls_header_line1(lines[0])
    n_per = 10
    el_rank_col_idx = [3 + i * n_per + 3 for i in range(len(alleles))]
    ba_col_idx = [3 + i * n_per + 4 for i in range(len(alleles))]
    ba_rank_col_idx = [3 + i * n_per + 5 for i in range(len(alleles))]

    rows_ba: list[list[float]] = []
    rows_br: list[list[float]] = []
    rows_elr: list[list[float]] = []
    peps: list[str] = []
    starts: list[int] = []
    for ln in lines[3:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 3 + n_per:
            continue
        pep = parts[1].strip().upper()
        if len(pep) != 9 or any(c not in STANDARD_AA for c in pep):
            continue
        pid = parts[2].strip()
        if not pid.isdigit():
            raise ValueError(f"Unexpected peptide ID (expected digits): {pid!r}")
        start = int(pid)
        vals_elr: list[float] = []
        vals_ba: list[float] = []
        vals_br: list[float] = []
        for j in el_rank_col_idx:
            if j >= len(parts):
                vals_elr.append(float("nan"))
            else:
                try:
                    vals_elr.append(float(parts[j]))
                except ValueError:
                    vals_elr.append(float("nan"))
        for j in ba_col_idx:
            if j >= len(parts):
                vals_ba.append(float("nan"))
            else:
                try:
                    vals_ba.append(float(parts[j]))
                except ValueError:
                    vals_ba.append(float("nan"))
        for j in ba_rank_col_idx:
            if j >= len(parts):
                vals_br.append(float("nan"))
            else:
                try:
                    vals_br.append(float(parts[j]))
                except ValueError:
                    vals_br.append(float("nan"))
        rows_ba.append(vals_ba)
        rows_br.append(vals_br)
        rows_elr.append(vals_elr)
        peps.append(pep)
        starts.append(start)
    return (
        alleles,
        starts,
        peps,
        np.asarray(rows_ba, dtype=np.float64),
        np.asarray(rows_br, dtype=np.float64),
        np.asarray(rows_elr, dtype=np.float64),
    )


def strong_binder_mask(
    ba: np.ndarray,
    ba_rank: np.ndarray,
    criterion: str,
    ic50_nm: float,
    ba_rank_pct: float,
    *,
    el_rank: np.ndarray | None = None,
    require_el_rank: bool = False,
    el_rank_pct: float | None = None,
) -> np.ndarray:
    """Shape (n_rows, n_alleles) boolean: predicted strong binder for that peptide–allele."""
    if require_el_rank:
        if el_rank is None:
            raise ValueError("require_el_rank is True but el_rank matrix is missing")
        el_cut = float(ba_rank_pct if el_rank_pct is None else el_rank_pct)
        el_ok = np.isfinite(el_rank) & (el_rank <= el_cut)
    else:
        el_ok = np.ones_like(ba_rank, dtype=bool)

    if criterion == "ic50":
        ic = np.vectorize(ba_score_to_ic50_nm)(ba)
        base = np.isfinite(ic) & (ic < ic50_nm)
        return base & el_ok
    if criterion == "ba_rank":
        base = np.isfinite(ba_rank) & (ba_rank <= ba_rank_pct)
        return base & el_ok
    raise ValueError(f"Unknown --sb-criterion: {criterion!r} (use ic50 or ba_rank)")


def per_position_metrics(
    full: str,
    starts: list[int],
    peps: list[str],
    sb: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    L = len(full)
    n_alleles = sb.shape[1]
    allele_cov = np.zeros(L, dtype=np.int32)
    epitope_cov = np.zeros(L, dtype=np.int32)
    for pos in range(L):
        ep_set: set[str] = set()
        al_set: set[int] = set()
        for si, s in enumerate(starts):
            if not (s <= pos < s + 9):
                continue
            row_sb = sb[si]
            if not row_sb.any():
                continue
            ep_set.add(peps[si])
            for ai in range(n_alleles):
                if row_sb[ai]:
                    al_set.add(ai)
        epitope_cov[pos] = len(ep_set)
        allele_cov[pos] = len(al_set)
    return epitope_cov, allele_cov


def per_position_sb_instance_load(full: str, starts: list[int], peps: list[str], sb: np.ndarray) -> np.ndarray:
    """
    Per residue: sum of SB peptide×allele hits over overlapping 9-mers (each allele prediction
    counts as one epitope instance at that site).
    """
    L = len(full)
    out = np.zeros(L, dtype=np.int32)
    for pos in range(L):
        t = 0
        for si, s in enumerate(starts):
            if not (s <= pos < s + 9):
                continue
            rw = sb[si]
            if rw.any():
                t += int(rw.sum())
        out[pos] = t
    return out


def build_ttn_long_for_iedb_merge(
    starts: list[int],
    peps: list[str],
    alleles: list[str],
    ba: np.ndarray,
    ba_rank: np.ndarray,
    el_rank: np.ndarray,
    parent_input_seq_id: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for si in range(len(starts)):
        pos0 = int(starts[si])
        pep = str(peps[si]).strip().upper()
        s1 = pos0 + 1
        e1 = pos0 + len(pep)
        for ai, al_xls in enumerate(alleles):
            al_i = normalize_hla_netmhc_to_iedb(str(al_xls))
            stable = f"{parent_input_seq_id}:{s1}:{e1}:{pep}:{al_i}"
            rows.append(
                {
                    "stable_key": stable,
                    "Peptide": pep,
                    "EL_rank": float(el_rank[si, ai]),
                    "BA_score": float(ba[si, ai]),
                    "_si": si,
                    "_ai": ai,
                }
            )
    return pd.DataFrame(rows)


def merge_iedb_csv(long_df: pd.DataFrame, iedb_csv: Path) -> pd.DataFrame:
    iedb = pd.read_csv(iedb_csv, low_memory=False)
    ren = {c: f"iedb_{c}" for c in iedb.columns}
    iedb2 = iedb.rename(columns=ren)
    out = long_df.merge(iedb2, left_on="stable_key", right_on="iedb_stable_key", how="left")
    return out.drop(columns=["iedb_stable_key"], errors="ignore")


def sb_matrix_from_iedb_merged(
    merged: pd.DataFrame,
    *,
    spec: SBSpec,
    iedb_ic50_col: str | None,
    n_si: int,
    n_ai: int,
) -> np.ndarray:
    bm = ba_score_min_for_ic50_lt(spec.ic50_max_nm)
    ok = sb_mask_spec(merged, spec, ba_min=bm, iedb_ic50_col=iedb_ic50_col)
    sb = np.zeros((n_si, n_ai), dtype=bool)
    for i in range(len(merged)):
        r = merged.iloc[i]
        sb[int(r["_si"]), int(r["_ai"])] = bool(ok.iloc[i])
    return sb


def collect_sb_epitopes(peps: list[str], sb: np.ndarray, allele_index: int | None = None) -> list[str]:
    """One copy of each 9-mer row with ≥1 SB (sensitivity / unique counts; logos use instance-weighted)."""
    out: list[str] = []
    for si, pep in enumerate(peps):
        if allele_index is None:
            hit = bool(sb[si].any())
        else:
            hit = bool(sb[si, allele_index])
        if hit:
            out.append(pep)
    return out


def collect_sb_epitopes_instance_weighted(peps: list[str], sb: np.ndarray, allele_index: int | None = None) -> list[str]:
    """
    For logos / composition: append each 9-mer once per **SB peptide×allele hit** (all alleles),
    or once per row if ``allele_index`` is set and that allele is SB (same as unweighted for a
    single boolean column).
    """
    out: list[str] = []
    for si, pep in enumerate(peps):
        if allele_index is None:
            k = int(np.sum(sb[si]))
        else:
            k = 1 if bool(sb[si, allele_index]) else 0
        if k:
            out.extend([pep] * k)
    return out


def freq_matrix(seqs: list[str]) -> pd.DataFrame:
    """Rows = epitope positions 1..9, columns = AA frequencies."""
    if not seqs:
        return pd.DataFrame(np.zeros((9, 20)), columns=AAS_ORDER)
    counts = np.zeros((9, 20), dtype=np.float64)
    for s in seqs:
        for i, c in enumerate(s):
            j = AAS_ORDER.index(c)
            counts[i, j] += 1.0
    row_sum = counts.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    prob = counts / row_sum
    return pd.DataFrame(prob, columns=AAS_ORDER)


_LOGO_COLOR_SCHEME_FIG6 = {
    "RKDENQ": [0, 0, 1],
    "SGHTAP": [0, 0.5, 0],
    "YVMCLFIW": list(to_rgb(pal.OI_ORANGE)),
}


def _logo_aa_color_for_scheme(aa: str) -> tuple[float, float, float]:
    for group, col in _LOGO_COLOR_SCHEME_FIG6.items():
        if aa in group:
            return (float(col[0]), float(col[1]), float(col[2]))
    return to_rgb(pal.OI_ORANGE)


_LOGO_FALLBACK_AA_COLORS = [_logo_aa_color_for_scheme(aa) for aa in AAS_ORDER]


def draw_logo(ax, seqs: list[str], title: str) -> None:
    df = freq_matrix(seqs)
    try:
        import logomaker

        logo = logomaker.Logo(
            df,
            ax=ax,
            color_scheme=_LOGO_COLOR_SCHEME_FIG6,
            vpad=0.05,
            width=0.85,
        )
        logo.style_spines(visible=False)
        ax.set_ylabel("Frequency")
    except Exception as exc:
        # Fallback when logomaker is missing or raises (e.g. bad matrix / version quirks).
        print(
            f"[plot_figure6_ttn_as1_allele_coverage] logomaker logo failed ({type(exc).__name__}: {exc!r}); "
            "using per-position stacked-bar fallback. Install with: pip install logomaker",
            flush=True,
        )
        # Fallback: stacked bars per position
        x = np.arange(9)
        bottom = np.zeros(9)
        for j, aa in enumerate(AAS_ORDER):
            h = df[aa].values
            ax.bar(
                x,
                h,
                bottom=bottom,
                color=_LOGO_FALLBACK_AA_COLORS[j],
                width=0.9,
                label=aa if h.max() > 0.35 else None,
            )
            bottom += h
        ax.set_ylabel("Frequency")
        ax.set_xticks(x)
        ax.set_xticklabels([str(i + 1) for i in x])
    ax.set_xlabel("Epitope position (1–9)")
    ax.set_title(title)
    ax.set_ylim(0, 1.0)


def stats_block_text(full: str, distinct_sb_alleles_per_site: np.ndarray, n_panel_alleles: int) -> str:
    """
    Summary for panel B: **distinct SB alleles per residue**, i.e. the union of allele indices
    with SB on any overlapping 9-mer (from ``per_position_metrics`` → ``allele_cov``). This is
    never greater than ``n_panel_alleles``.
    """
    cov = distinct_sb_alleles_per_site
    covered = int((cov > 0).sum())
    L = len(full)
    mx = int(cov.max())
    return (
        f"Total positions: {L}\n"
        f"Covered positions: {covered} ({100.0 * covered / L:.1f}%)\n"
        f"Uncovered positions: {L - covered} ({100.0 * (L - covered) / L:.1f}%)\n"
        f"Max distinct SB alleles: {mx}\n"
        f"Mean distinct SB alleles: {cov.mean():.2f}\n"
        f"Median distinct SB alleles: {float(np.median(cov)):.1f}\n"
        f"HLA panel size: {n_panel_alleles} alleles"
    )


def panel_sequence_grid(
    ax,
    full: str,
    allele_cov: np.ndarray,
    title: str,
    subtitle: str,
    *,
    colorbar_label: str = "Allele coverage (unique)",
) -> None:
    L = len(full)
    vmax = max(int(allele_cov.max()), 1)
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = plt.get_cmap(pal.sequential_heatmap())
    n_per_row = 41
    rows = int(np.ceil(L / n_per_row))
    ax.set_xlim(0, n_per_row)
    ax.set_ylim(0, rows)
    ax.invert_yaxis()
    for i, aa in enumerate(full):
        r, c = divmod(i, n_per_row)
        col = cmap(norm(int(allele_cov[i])))
        ax.add_patch(Rectangle((c, r), 0.95, 0.95, facecolor=col, edgecolor="0.5", linewidth=0.3))
        ax.text(c + 0.47, r + 0.5, aa, ha="center", va="center", fontsize=7, color="0.1")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title + "\n" + subtitle, fontsize=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(colorbar_label)
    cb.ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{int(round(v))}"))
    cb.ax.tick_params(axis="x", bottom=False, labelbottom=False)


def metric_output_path(out: Path, metric_tag: str) -> Path:
    """Insert ``_instances`` / ``_unique`` before extension unless already present."""
    parent, stem, suf = out.parent, out.stem, out.suffix
    if stem.endswith(f"_{metric_tag}"):
        return out
    return parent / f"{stem}_{metric_tag}{suf}"


def _render_figure6_single_combined(
    out: Path,
    full: str,
    alleles: list[str],
    subtitle: str,
    sb_all_logo: list[str],
    sb_ex_logo: list[str],
    ex_norm: str,
    *,
    heatmap_cov: np.ndarray,
    hist_cov: np.ndarray,
    c1_y: np.ndarray,
    c2_y: np.ndarray,
    c3_y1: np.ndarray,
    c3_y2: np.ndarray,
    c3_leg: tuple[str, str],
    c1_title: str,
    c2_title: str,
    c1_ylabel: str,
    c2_ylabel: str,
    hist_xlabel: str,
    heatmap_cb_label: str,
    stats_cov: np.ndarray,
    metric_suffix: str,
    stats_txt: str | None = None,
    publication_dir: Path | None = None,
    publication_tiff_kind: str = "color",
    figures_root: Path = REPO_FIGURES,
) -> None:
    stats_use = stats_txt if stats_txt is not None else stats_block_text(full, stats_cov, len(alleles))
    fig = plt.figure(figsize=(12, 14), dpi=150)
    outer = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.15, 1.0], hspace=0.32)
    gs_top = outer[0].subgridspec(1, 2, width_ratios=[2.1, 1.0], wspace=0.15)
    ax_a = fig.add_subplot(gs_top[0, 0])
    ax_b = fig.add_subplot(gs_top[0, 1])
    panel_sequence_grid(
        ax_a,
        full,
        heatmap_cov,
        f"Amino acid sequence (TTN-AS1, {metric_suffix})",
        subtitle,
        colorbar_label=heatmap_cb_label,
    )
    bins = np.arange(-0.5, float(np.max(hist_cov)) + 1.5, 1.0)
    ax_b.hist(hist_cov, bins=bins, color=pal.HIST_FILL, edgecolor="white")
    ax_b.set_xlabel(hist_xlabel)
    ax_b.set_ylabel("Number of positions")
    ax_b.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax_b.xaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax_b.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{int(round(v))}"))
    ax_b.set_title("Distribution of coverage across positions")
    ax_b.text(0.98, 0.97, stats_use, transform=ax_b.transAxes, va="top", ha="right", fontsize=8, family="monospace")

    gsc = outer[1].subgridspec(3, 1, hspace=0.08)
    x = np.arange(1, len(full) + 1)
    ax_c1 = fig.add_subplot(gsc[0, 0])
    ax_c2 = fig.add_subplot(gsc[1, 0], sharex=ax_c1)
    ax_c3 = fig.add_subplot(gsc[2, 0], sharex=ax_c1)
    ax_c1.plot(x, c1_y, color=pal.LINE_PRIMARY, lw=1.2)
    ax_c1.set_ylabel(c1_ylabel)
    ax_c1.set_title(c1_title)
    ax_c1.grid(True, alpha=0.25)
    ax_c2.plot(x, c2_y, color=pal.LINE_SECONDARY, lw=1.2)
    ax_c2.set_ylabel(c2_ylabel)
    ax_c2.set_title(c2_title)
    ax_c2.grid(True, alpha=0.25)
    ax_c3.plot(x, c3_y1, color=pal.LINE_PRIMARY, lw=1.0, label=c3_leg[0])
    ax_c3.plot(x, c3_y2, color=pal.LINE_SECONDARY, lw=1.0, label=c3_leg[1])
    ax_c3.set_ylabel("Overlay")
    ax_c3.set_xlabel("Position")
    ax_c3.legend(loc="upper right", fontsize=8)
    ax_c3.grid(True, alpha=0.25)
    for _ax in (ax_c1, ax_c2, ax_c3):
        _ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    plt.setp(ax_c1.get_xticklabels(), visible=False)
    plt.setp(ax_c2.get_xticklabels(), visible=False)

    gs_bot = outer[2].subgridspec(1, 2, wspace=0.25)
    ax_d = fig.add_subplot(gs_bot[0, 0])
    ax_e = fig.add_subplot(gs_bot[0, 1])
    draw_logo(
        ax_d,
        sb_all_logo,
        f"9-mer epitope logo (instance-weighted, N={len(sb_all_logo)})",
    )
    draw_logo(
        ax_e,
        sb_ex_logo,
        f"9-mer epitope logo — {display_allele(ex_norm)} (instance-weighted, N={len(sb_ex_logo)})",
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    save_figure_bundle(
        fig,
        out,
        png_dpi=150,
        publication_dir=publication_dir,
        publication_tiff_kind=publication_tiff_kind,
        figures_root=figures_root,
        bbox_inches="tight",
    )
    plt.close(fig)
    print(f"Wrote {out}")


def _render_figure6_split_panels(
    out: Path,
    full: str,
    alleles: list[str],
    subtitle: str,
    sb_all_logo: list[str],
    sb_ex_logo: list[str],
    ex_norm: str,
    *,
    heatmap_cov: np.ndarray,
    hist_cov: np.ndarray,
    c1_y: np.ndarray,
    c2_y: np.ndarray,
    c3_y1: np.ndarray,
    c3_y2: np.ndarray,
    c3_leg: tuple[str, str],
    c1_title: str,
    c2_title: str,
    c1_ylabel: str,
    c2_ylabel: str,
    hist_xlabel: str,
    heatmap_cb_label: str,
    stats_cov: np.ndarray,
    metric_suffix: str,
    stats_txt: str | None = None,
    publication_dir: Path | None = None,
    publication_tiff_kind: str = "color",
    figures_root: Path = REPO_FIGURES,
) -> None:
    if out.suffix.lower() != ".png":
        raise SystemExit("--split-panels requires -o to end with .png")
    stem, suf, parent = out.stem, out.suffix, out.parent

    def _p(letter: str) -> Path:
        return parent / f"{stem}_{letter}{suf}"

    def _save_pub(fig, path_png: Path) -> None:
        save_figure_bundle(
            fig,
            path_png,
            png_dpi=150,
            publication_dir=publication_dir,
            publication_tiff_kind=publication_tiff_kind,
            figures_root=figures_root,
            bbox_inches="tight",
        )

    stats_use = stats_txt if stats_txt is not None else stats_block_text(full, stats_cov, len(alleles))

    fig_a, ax_a = plt.subplots(figsize=(14, 4.8), dpi=150)
    panel_sequence_grid(
        ax_a,
        full,
        heatmap_cov,
        f"TTN-AS1 sequence ({metric_suffix})",
        subtitle,
        colorbar_label=heatmap_cb_label,
    )
    _save_pub(fig_a, _p("a"))
    plt.close(fig_a)

    fig_b, ax_b = plt.subplots(figsize=(7, 5.5), dpi=150)
    bins = np.arange(-0.5, float(np.max(hist_cov)) + 1.5, 1.0)
    ax_b.hist(hist_cov, bins=bins, color=pal.HIST_FILL, edgecolor="white")
    ax_b.set_xlabel(hist_xlabel)
    ax_b.set_ylabel("Number of positions")
    ax_b.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax_b.xaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    ax_b.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{int(round(v))}"))
    ax_b.set_title("Distribution across positions")
    ax_b.text(0.98, 0.97, stats_use, transform=ax_b.transAxes, va="top", ha="right", fontsize=8, family="monospace")
    _save_pub(fig_b, _p("b"))
    plt.close(fig_b)

    fig_c = plt.figure(figsize=(12, 7.5), dpi=150)
    gsc = fig_c.add_gridspec(3, 1, hspace=0.12)
    x = np.arange(1, len(full) + 1)
    ax_c1 = fig_c.add_subplot(gsc[0, 0])
    ax_c2 = fig_c.add_subplot(gsc[1, 0], sharex=ax_c1)
    ax_c3 = fig_c.add_subplot(gsc[2, 0], sharex=ax_c1)
    ax_c1.plot(x, c1_y, color=pal.LINE_PRIMARY, lw=1.2)
    ax_c1.set_ylabel(c1_ylabel)
    ax_c1.set_title(c1_title)
    ax_c1.grid(True, alpha=0.25)
    ax_c2.plot(x, c2_y, color=pal.LINE_SECONDARY, lw=1.2)
    ax_c2.set_ylabel(c2_ylabel)
    ax_c2.set_title(c2_title)
    ax_c2.grid(True, alpha=0.25)
    ax_c3.plot(x, c3_y1, color=pal.LINE_PRIMARY, lw=1.0, label=c3_leg[0])
    ax_c3.plot(x, c3_y2, color=pal.LINE_SECONDARY, lw=1.0, label=c3_leg[1])
    ax_c3.set_ylabel("Overlay")
    ax_c3.set_xlabel("Position")
    ax_c3.legend(loc="upper right", fontsize=8)
    ax_c3.grid(True, alpha=0.25)
    for _ax in (ax_c1, ax_c2, ax_c3):
        _ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=3))
    plt.setp(ax_c1.get_xticklabels(), visible=False)
    plt.setp(ax_c2.get_xticklabels(), visible=False)
    _save_pub(fig_c, _p("c"))
    plt.close(fig_c)

    fig_d, ax_d = plt.subplots(figsize=(8, 4.2), dpi=150)
    draw_logo(
        ax_d,
        sb_all_logo,
        f"9-mer epitope logo (instance-weighted, N={len(sb_all_logo)})",
    )
    _save_pub(fig_d, _p("d"))
    plt.close(fig_d)

    fig_e, ax_e = plt.subplots(figsize=(8, 4.2), dpi=150)
    draw_logo(
        ax_e,
        sb_ex_logo,
        f"9-mer epitope logo — {display_allele(ex_norm)} (instance-weighted, N={len(sb_ex_logo)})",
    )
    _save_pub(fig_e, _p("e"))
    plt.close(fig_e)
    for ch in "abcde":
        print(f"Wrote {_p(ch)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--gating",
        choices=("netmhc", "iedb_sb"),
        default="iedb_sb",
        help="netmhc: SB from wide XLS only. iedb_sb: merge XLS rows to IEDB CSV on stable_key "
        "and apply cohort-style SB (immuno/processing/EL/IC50 via netmhc_sb_core). Default: iedb_sb.",
    )
    ap.add_argument(
        "--netmhc-xls",
        type=Path,
        default=NET / "netmhcpan_ttn_as1_108065.xls",
        help="NetMHCpan wide XLS (-xls 1 -BA 1)",
    )
    ap.add_argument(
        "--parent-fasta",
        type=Path,
        default=None,
        help="Optional single-record FASTA for the 79 aa parent (else use SmProt export)",
    )
    ap.add_argument(
        "--iedb-csv",
        type=Path,
        default=TTN_IEDB_SYNTHETIC_CSV,
        help="IEDB peptide_table CSV (with stable_key) when --gating iedb_sb (default: synthetic companion).",
    )
    ap.add_argument(
        "--iedb-parent-input-seq-id",
        type=str,
        default=TTN_IEDB_PARENT_SEQ_ID_DEFAULT,
        help="input_seq_id token matching IEDB stable_key prefix (default for synthetic TTN-AS1 companion).",
    )
    ap.add_argument(
        "--imm-min",
        type=float,
        default=FIG5_IEDB_IMM_MIN_DEFAULT,
        help="IEDB-gated SB: immunogenicity score > this (default matches Fig 5 netmhc_sb_core)",
    )
    ap.add_argument(
        "--proc-min",
        type=float,
        default=FIG5_IEDB_PROC_MIN_DEFAULT,
        help="IEDB-gated SB: processing score > this (default matches Fig 5 netmhc_sb_core)",
    )
    ap.add_argument(
        "--el-rank-max",
        type=float,
        default=FIG5_IEDB_EL_RANK_MAX_DEFAULT,
        help="IEDB-gated SB: EL %%rank cutoff vs NetMHC EL_rank column (default matches Fig 5)",
    )
    ap.add_argument(
        "--el-rank-lte",
        action="store_true",
        help="IEDB-gated SB: use EL_rank ≤ --el-rank-max (default: strict <)",
    )
    ap.add_argument(
        "--ic50-max-nm",
        type=float,
        default=FIG5_IEDB_IC50_MAX_NM_DEFAULT,
        help="IEDB-gated SB: IC50 upper bound (nM) when IEDB IC50 column is available (default matches Fig 5)",
    )
    ap.add_argument(
        "--sb-mode",
        choices=("full", "no_ic50", "ic50_only"),
        default="full",
        help="IEDB SB: full stack; no_ic50 (IEDB+EL, no binding gate); ic50_only (binding/IC50 only).",
    )
    ap.add_argument(
        "--sb-criterion",
        choices=("ba_rank", "ic50"),
        default="ba_rank",
        help="How to define strong binders (default: BA_rank in top %%ile band)",
    )
    ap.add_argument(
        "--ba-rank-pct",
        type=float,
        default=0.5,
        help="SB if BA_rank <= this value (NetMHCpan %% rank; default 0.5 = top 0.5%%)",
    )
    ap.add_argument(
        "--require-el-rank",
        action="store_true",
        help="Also require EL_rank <= cutoff (default same as --ba-rank-pct unless --el-rank-pct)",
    )
    ap.add_argument(
        "--el-rank-pct",
        type=float,
        default=None,
        help="EL_rank cutoff when --require-el-rank (default: same as --ba-rank-pct)",
    )
    ap.add_argument(
        "--ic50-nm",
        type=float,
        default=150.0,
        help="With --sb-criterion ic50: SB if predicted IC50 < this (nM)",
    )
    ap.add_argument(
        "--example-allele",
        type=str,
        default="HLA-A*32:01",
        help="Allele for panel E (logo of SB 9-mers for that allele only)",
    )
    ap.add_argument(
        "--coverage-output",
        choices=("instances", "unique", "both"),
        default="instances",
        help="instances: A–C use distinct SB alleles/epitopes per site; D/E logos stay instance-weighted. "
        "unique: same A–C curves as instances; both: write *_instances* and *_unique* stems.",
    )
    ap.add_argument(
        "--also-write-unique",
        action="store_true",
        help="When --coverage-output is instances, also write the *_unique companion files.",
    )
    ap.add_argument(
        "-o",
        "--out",
        type=Path,
        default=REPO_FIGURES / "fig6_ttn_as1_allele_coverage.png",
        help="Output PNG (or stem for --split-panels). Default: figures/fig6_ttn_as1_allele_coverage.png",
    )
    ap.add_argument(
        "--split-panels",
        action="store_true",
        help="Write five PNGs: {stem}_a.png … {stem}_e.png (uses -o path stem/suffix/parent)",
    )
    add_publication_args(ap)
    args = ap.parse_args()

    if args.gating == "iedb_sb":
        if not args.iedb_csv or not args.iedb_csv.is_file():
            raise SystemExit("--gating iedb_sb requires a readable --iedb-csv")
        if not args.iedb_parent_input_seq_id:
            raise SystemExit("--gating iedb_sb requires --iedb-parent-input-seq-id")

    if args.parent_fasta:
        from Bio import SeqIO

        rec = next(SeqIO.parse(args.parent_fasta, "fasta"))
        full = str(rec.seq).upper().strip()
    else:
        full = (
            "ISATDRICENTSMSRLGIILRHHLASPASHFKMIANDSTSSITDWLIPLYFHAVPGGQCDNWSARRTRNFEWILGYSRL"
        )

    alleles, starts, peps, ba, ba_rank, el_rank = parse_wide_netmhc_xls_rows(args.netmhc_xls)

    if args.gating == "netmhc":
        sb = strong_binder_mask(
            ba,
            ba_rank,
            args.sb_criterion,
            args.ic50_nm,
            args.ba_rank_pct,
            el_rank=el_rank,
            require_el_rank=args.require_el_rank,
            el_rank_pct=args.el_rank_pct,
        )
    else:
        long_df = build_ttn_long_for_iedb_merge(
            starts,
            peps,
            alleles,
            ba,
            ba_rank,
            el_rank,
            str(args.iedb_parent_input_seq_id).strip(),
        )
        merged = merge_iedb_csv(long_df, args.iedb_csv)
        n_hit = int(merged["iedb_score"].notna().sum()) if "iedb_score" in merged.columns else 0
        if n_hit == 0:
            print(
                "[plot_figure6_ttn_as1_allele_coverage] Warning: no IEDB rows matched stable_key; "
                "SB matrix will be all false unless --iedb-parent-input-seq-id matches IEDB input_seq_id.",
                flush=True,
            )
        spec = sb_spec_from_mode(
            args.sb_mode,
            imm_min=float(args.imm_min),
            proc_min=float(args.proc_min),
            el_max=float(args.el_rank_max),
            el_lte=bool(args.el_rank_lte),
            ic50_max_nm=float(args.ic50_max_nm),
        )
        iedb_ic50_col = pick_iedb_ic50_column(merged.columns)
        ba_min = ba_score_min_for_ic50_lt(spec.ic50_max_nm)
        sb = sb_matrix_from_iedb_merged(
            merged,
            spec=spec,
            iedb_ic50_col=iedb_ic50_col,
            n_si=len(starts),
            n_ai=len(alleles),
        )

    _, allele_cov = per_position_metrics(full, starts, peps, sb)
    inst_load = per_position_sb_instance_load(full, starts, peps, sb)

    ex_norm = normalize_allele(args.example_allele)
    if ex_norm not in alleles:
        raise SystemExit(f"Example allele {args.example_allele!r} not in XLS (after normalize).")
    ex_i = alleles.index(ex_norm)

    sb_all_logo = collect_sb_epitopes_instance_weighted(peps, sb, None)
    sb_ex_logo = collect_sb_epitopes_instance_weighted(peps, sb, ex_i)
    n_inst = len(sb_all_logo)
    n_seq_u = len(set(sb_all_logo))

    subtitle = (
        f"Length: {len(full)} AA | Total unique alleles in panel: {len(alleles)} | "
        f"SB epitope–allele instances (logo weights): {n_inst} "
        f"({n_seq_u} distinct 9-mer sequences)"
    )

    want_instances = args.coverage_output in ("instances", "both")
    want_unique = args.coverage_output in ("unique", "both") or (
        args.coverage_output == "instances" and args.also_write_unique
    )

    pack_instances = dict(
        metric_suffix="instances",
        heatmap_cov=allele_cov,
        hist_cov=allele_cov,
        c1_y=inst_load,
        c2_y=allele_cov,
        c3_y1=inst_load,
        c3_y2=allele_cov,
        c3_leg=("SB epitope instances (hits)", "Distinct SB alleles"),
        c1_title="SB epitope instance hits per site (sum over overlapping 9-mer×allele SB)",
        c2_title="Distinct SB alleles per site",
        c1_ylabel="Instance hits",
        c2_ylabel="Unique alleles",
        hist_xlabel="Distinct SB alleles per position",
        heatmap_cb_label="Distinct SB alleles (per site)",
        stats_cov=allele_cov,
        stats_txt=stats_block_text(full, allele_cov, len(alleles)),
    )
    pack_unique = dict(
        metric_suffix="unique",
        heatmap_cov=allele_cov,
        hist_cov=allele_cov,
        c1_y=inst_load,
        c2_y=allele_cov,
        c3_y1=inst_load,
        c3_y2=allele_cov,
        c3_leg=("SB epitope instances (hits)", "Distinct SB alleles"),
        c1_title="SB epitope instance hits per site (sum over overlapping 9-mer×allele SB)",
        c2_title="Distinct SB alleles per site",
        c1_ylabel="Instance hits",
        c2_ylabel="Unique alleles",
        hist_xlabel="Number of distinct alleles (SB at position)",
        heatmap_cb_label="Unique alleles (SB)",
        stats_cov=allele_cov,
    )

    jobs: list[tuple[str, dict[str, object]]] = []
    if want_instances:
        jobs.append(("instances", pack_instances))
    if want_unique:
        jobs.append(("unique", pack_unique))
    if not jobs:
        raise SystemExit("No coverage outputs selected (internal error).")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    for tag, pack in jobs:
        outp = metric_output_path(args.out, tag)
        kw = dict(pack)
        stats_txt_kw = kw.pop("stats_txt", None)
        if args.split_panels:
            _render_figure6_split_panels(
                outp,
                full,
                alleles,
                subtitle,
                sb_all_logo,
                sb_ex_logo,
                ex_norm,
                stats_txt=stats_txt_kw,
                publication_dir=args.publication_dir,
                publication_tiff_kind=args.publication_tiff_kind,
                figures_root=REPO_FIGURES,
                **kw,
            )
        else:
            _render_figure6_single_combined(
                outp,
                full,
                alleles,
                subtitle,
                sb_all_logo,
                sb_ex_logo,
                ex_norm,
                stats_txt=stats_txt_kw,
                publication_dir=args.publication_dir,
                publication_tiff_kind=args.publication_tiff_kind,
                figures_root=REPO_FIGURES,
                **kw,
            )


if __name__ == "__main__":
    main()
