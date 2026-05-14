"""
Shared SB (stringent binding-style) filter spec for NetMHC + IEDB merged TSVs.

Used by ``plot_fig5abc_netmhc_sb_triple.py``, ``netmhc_sb_sensitivity_robustness.py``,
``plot_fig5_netmhc_sb_combination_grid.py``, ``netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py``,
``plot_fig6_ttn_merged_iedb_sb_combination_grid.py``, and ``plot_figure6_ttn_as1_allele_coverage.py`` (IEDB SB).
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

import pandas as pd

LOG10_50000 = math.log10(50000.0)

# Manuscript defaults for merged IEDB + NetMHC SB (Fig 5ABC, sensitivity / combo baselines, TTN iedb_sb).
# Immuno / processing align with the common literature-style IEDB cuts; EL and IC50 are still
# tightened vs an unfiltered screen (EL %%rank top band + strict IC50 < threshold).
FIG5_IEDB_IMM_MIN_DEFAULT = 0.1
FIG5_IEDB_PROC_MIN_DEFAULT = 1.5
FIG5_IEDB_EL_RANK_MAX_DEFAULT = 1.0
FIG5_IEDB_IC50_MAX_NM_DEFAULT = 150.0


def ba_score_min_for_ic50_lt(ic50_max_nm: float) -> float:
    """Minimum BA_score so that IC50 < ic50_max_nm (strict upper bound on IC50)."""
    return 1.0 - math.log10(float(ic50_max_nm)) / LOG10_50000


def pick_iedb_ic50_column(columns: Iterable[str]) -> str | None:
    cols = list(columns)
    preferred = "iedb_netmhcpan_ba_ic50"
    if preferred in cols:
        return preferred
    cands = [c for c in cols if c.startswith("iedb_") and "ic50" in c.lower()]
    if not cands:
        return None
    cands.sort(key=lambda s: (0 if "netmhcpan" in s.lower() and "ba" in s.lower() else 1, s))
    return cands[0]


@dataclass(frozen=True)
class SBSpec:
    imm_min: float
    proc_min: float
    el_max: float
    el_lte: bool
    ic50_max_nm: float
    use_imm: bool = True
    use_proc: bool = True
    use_el: bool = True
    use_ic50: bool = True


def sb_spec_from_mode(
    sb_mode: str,
    *,
    imm_min: float,
    proc_min: float,
    el_max: float,
    el_lte: bool,
    ic50_max_nm: float,
) -> SBSpec:
    """
    Build ``SBSpec`` for merged cohort / TTN IEDB joins.

    - ``full``: IEDB immuno + processing + NetMHC EL_rank + IC50/BA binding (default Fig 5 stack).
    - ``no_ic50``: same as ``full`` but **no** IC50 / BA binding gate (presentation-only IEDB + EL).
    - ``ic50_only``: **only** predicted binding (IEDB ``netmhcpan_ba_ic50`` if present, else local
      ``BA_score`` vs ``ic50_max_nm``); IEDB immuno / processing / EL gates are **off**.
    """
    m = str(sb_mode).strip().lower()
    if m == "full":
        return SBSpec(
            imm_min=imm_min,
            proc_min=proc_min,
            el_max=el_max,
            el_lte=el_lte,
            ic50_max_nm=ic50_max_nm,
        )
    if m == "no_ic50":
        return SBSpec(
            imm_min=imm_min,
            proc_min=proc_min,
            el_max=el_max,
            el_lte=el_lte,
            ic50_max_nm=ic50_max_nm,
            use_ic50=False,
        )
    if m == "ic50_only":
        return SBSpec(
            imm_min=imm_min,
            proc_min=proc_min,
            el_max=el_max,
            el_lte=el_lte,
            ic50_max_nm=ic50_max_nm,
            use_imm=False,
            use_proc=False,
            use_el=False,
            use_ic50=True,
        )
    raise ValueError(f"Unknown sb_mode {sb_mode!r} (use full, no_ic50, ic50_only)")


def sb_mask_spec(
    df: pd.DataFrame,
    spec: SBSpec,
    *,
    ba_min: float,
    iedb_ic50_col: str | None,
) -> pd.Series:
    el = pd.to_numeric(df["EL_rank"], errors="coerce")
    imm = pd.to_numeric(df["iedb_score"], errors="coerce")
    proc = pd.to_numeric(df["iedb_processing_score"], errors="coerce")
    pep = df["Peptide"].astype(str)

    ok = ~pep.str.contains("X", case=False, na=False)
    if spec.use_imm or spec.use_proc:
        ok &= imm.notna() & proc.notna()
    if spec.use_el:
        ok &= el.notna()

    if spec.use_imm:
        ok &= imm > spec.imm_min
    if spec.use_proc:
        ok &= proc > spec.proc_min
    if spec.use_el:
        el_cmp = (el <= spec.el_max) if spec.el_lte else (el < spec.el_max)
        ok &= el_cmp

    if spec.use_ic50:
        if iedb_ic50_col is not None:
            ic50 = pd.to_numeric(df[iedb_ic50_col], errors="coerce")
            ok &= ic50.notna() & (ic50 < float(spec.ic50_max_nm))
        else:
            ba = pd.to_numeric(df["BA_score"], errors="coerce")
            ok &= ba.notna() & (ba > ba_min)
    return ok


def sb_mask_fig5_defaults(
    df: pd.DataFrame,
    *,
    el_max: float,
    el_lte: bool,
    ba_min: float,
    ic50_max_nm: float,
    iedb_ic50_col: str | None,
    imm_min: float = FIG5_IEDB_IMM_MIN_DEFAULT,
    proc_min: float = FIG5_IEDB_PROC_MIN_DEFAULT,
    sb_mode: str = "full",
) -> pd.Series:
    """SB mask for merged Fig 5 tables. ``sb_mode``: see ``sb_spec_from_mode``."""
    spec = sb_spec_from_mode(
        sb_mode,
        imm_min=imm_min,
        proc_min=proc_min,
        el_max=el_max,
        el_lte=el_lte,
        ic50_max_nm=ic50_max_nm,
    )
    return sb_mask_spec(df, spec, ba_min=ba_min, iedb_ic50_col=iedb_ic50_col)


def spec_label(spec: SBSpec) -> str:
    elop = "<=" if spec.el_lte else "<"
    parts = []
    if spec.use_imm:
        parts.append(f"imm>{spec.imm_min:g}")
    else:
        parts.append("imm=OFF")
    if spec.use_proc:
        parts.append(f"proc>{spec.proc_min:g}")
    else:
        parts.append("proc=OFF")
    if spec.use_el:
        parts.append(f"EL{elop}{spec.el_max:g}%")
    else:
        parts.append("EL=OFF")
    if spec.use_ic50:
        parts.append(f"IC50<{spec.ic50_max_nm:g}nM")
    else:
        parts.append("IC50=OFF")
    return ", ".join(parts)


def spec_profile_id(spec: SBSpec) -> str:
    """Stable filesystem token (no spaces)."""
    elop = "le" if spec.el_lte else "lt"
    return (
        f"imm{spec.imm_min:g}_proc{spec.proc_min:g}_el{elop}{spec.el_max:g}_ic50lt{spec.ic50_max_nm:g}"
        .replace(".", "p")
    )
