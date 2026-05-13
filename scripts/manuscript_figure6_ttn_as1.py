#!/usr/bin/env python3
"""
Catalog **Figure 6** (TTN-AS1, smPEP 108065): forwards CLI to the implementation module.

Implementation lives at repo root: ``plot_figure6_ttn_as1_allele_coverage.py`` (historical path
imports such as ``build_ttn_iedb_companion_csv`` depend on). Run from repository root::

    python scripts/manuscript_figure6_ttn_as1.py --split-panels -o figures/fig6_ttn_as1_split.png
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    path = ROOT / "plot_figure6_ttn_as1_allele_coverage.py"
    spec = importlib.util.spec_from_file_location("_fig6_ttn_as1_impl", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    mod.main()


if __name__ == "__main__":
    main()
