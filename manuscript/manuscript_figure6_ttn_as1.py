#!/usr/bin/env python3
"""
Catalog **Figure 6** (TTN-AS1, smPEP 108065): forwards CLI to the implementation module.

Implementation lives in ``manuscript/plot_figure6_ttn_as1_allele_coverage.py``. Run from repository root::

    python manuscript/manuscript_figure6_ttn_as1.py --split-panels -o figures/fig6_ttn_as1_split.png
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


import importlib.util

def main() -> None:
    path = ROOT / "manuscript" / "plot_figure6_ttn_as1_allele_coverage.py"
    spec = importlib.util.spec_from_file_location("_fig6_ttn_as1_impl", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    mod.main()


if __name__ == "__main__":
    main()
