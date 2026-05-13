"""Guard against drift between z-score Python filters and limma R filters (shared literals)."""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _load_tr_py():
    path = REPO / "pipeline" / "tr_lncrna_de_analysis.py"
    spec = importlib.util.spec_from_file_location("tr_lncrna_de_analysis", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _r_float(name: str, text: str) -> float:
    m = re.search(rf"^{re.escape(name)}\s*<-\s*([0-9.]+)\s*$", text, re.MULTILINE)
    assert m, f"missing R assignment {name}"
    return float(m.group(1))


def _r_int(name: str, text: str) -> int:
    m = re.search(rf"^{re.escape(name)}\s*<-\s*(\d+)L?\s*$", text, re.MULTILINE)
    assert m, f"missing R assignment {name}"
    return int(m.group(1))


def _r_transition_names_in_block(text: str, start_marker: str, end_marker: str) -> list[str]:
    i = text.index(start_marker)
    j = text.index(end_marker, i)
    return re.findall(r'name\s*=\s*"([^"]+)"', text[i:j])


def test_shared_filter_constants_match_r() -> None:
    py = _load_tr_py()
    r_text = (REPO / "tr_limma_de.R").read_text(encoding="utf-8")
    assert py.LOG2_THRESH == _r_float("LOG2_THRESH", r_text)
    assert py.EXPR_FRAC_MIN == _r_float("EXPR_FRAC_MIN", r_text)
    assert py.MIN_SAMPLES_CANCER == _r_int("MIN_SAMPLES_CANCER", r_text)


def test_stage_and_m_transition_names_match_r() -> None:
    py = _load_tr_py()
    r_text = (REPO / "tr_limma_de.R").read_text(encoding="utf-8")
    py_stage = [t[0] for t in py.STAGE_TRANSITIONS]
    py_m = [t[0] for t in py.M_TRANSITIONS]
    r_stage = _r_transition_names_in_block(
        r_text, "STAGE_TRANSITIONS <- list(", "M_TRANSITIONS <- list("
    )
    r_m = _r_transition_names_in_block(
        r_text, "M_TRANSITIONS <- list(", "cancer_types_over_n <- function"
    )
    assert py_stage == r_stage
    assert py_m == r_m
