from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "scripts"))

from tcga_cohort_filters import MIN_SAMPLES_CANCER, cancer_types_over_n, filter_samples_min_cancer_count


def test_cancer_types_over_n_strict_gt() -> None:
    meta = pd.DataFrame({"cancer_type": ["A"] * 100 + ["B"] * 101})
    assert cancer_types_over_n(meta, 100) == ["B"]


def test_filter_zero_keeps_all() -> None:
    meta = pd.DataFrame({"cancer_type": ["A", "B"], "x": [1, 2]})
    out, types = filter_samples_min_cancer_count(meta, min_samples_per_cancer=0)
    assert len(out) == 2
    assert set(types) == {"A", "B"}


def test_default_matches_tr_constant() -> None:
    assert MIN_SAMPLES_CANCER == 100
