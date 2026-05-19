"""
TCGA primary-tumor cohort filters shared with tr_lncrna_de_analysis / tr_limma_de.R.

Used by sample embedding scripts (Fig 1B t-SNE, supplement PCA) so plots use the same
cancer-type sample-size rule as DE / peptide-fraction bar charts (> N samples per type).
"""

from __future__ import annotations

import pandas as pd

# Match pipeline/tr_lncrna_de_analysis.py and tr_limma_de.R (strictly > N).
MIN_SAMPLES_CANCER = 100


def cancer_types_over_n(meta: pd.DataFrame, n: int = MIN_SAMPLES_CANCER) -> list[str]:
    """Cancer types with more than ``n`` samples in ``meta``."""
    if "cancer_type" not in meta.columns:
        raise ValueError("meta must include a cancer_type column")
    counts = meta.groupby("cancer_type", observed=True).size()
    return counts[counts > n].index.astype(str).tolist()


def filter_samples_min_cancer_count(
    df: pd.DataFrame,
    *,
    min_samples_per_cancer: int = MIN_SAMPLES_CANCER,
    cancer_col: str = "cancer_type",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Restrict to samples whose ``cancer_col`` value has > ``min_samples_per_cancer`` rows in ``df``.

    Returns (filtered dataframe, sorted list of retained cancer types).
    Set ``min_samples_per_cancer=0`` to keep all samples.
    """
    if min_samples_per_cancer <= 0:
        types = sorted(df[cancer_col].astype(str).unique())
        return df.copy(), types
    keep_types = cancer_types_over_n(df[[cancer_col]], min_samples_per_cancer)
    if not keep_types:
        mx = int(df.groupby(cancer_col, observed=True).size().max())
        raise ValueError(
            f"No cancer types with >{min_samples_per_cancer} samples (largest type has {mx})"
        )
    out = df.loc[df[cancer_col].astype(str).isin(keep_types)].copy()
    return out, sorted(keep_types)
