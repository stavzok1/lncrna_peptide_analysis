#!/usr/bin/env python3
"""
Regenerate the **full** recommended figure pipeline from the repository root.

This is a thin orchestrator over the steps documented in ``docs/figure_generation_overview.md``
(§ “Full regeneration”). It does **not** rebuild upstream ``data/`` tables (expression matrices,
NetMHC merges, FASTAs); ensure those inputs exist first.

Default order:

1. ``generate_catalog_figures.py`` — Fig 1B + catalog Fig 2–3 for both peptide modes + Fig 3C–3D + 4A
2. ``generate_netmhc_figure_bundle.py`` — merged NetMHC Fig 5–6 (+ random-fragment mirrors unless skipped)
3. ``generate_netmhc_fig5_fig6_supplement.py`` — Fig 5–6 supplement tree (1D LOO, Cartesian grids, TTN sweeps)
4. ``generate_canonical_manuscript_figures.py`` — OpenTSNE Fig 1B supplement, main-text NetMHC subset (pass this orchestrator ``--include-fig6-unique`` to also emit Fig 6 unique under ``figures/supplementary/figure6_ttn_as1/``)
5. ``export_publication_figures.py`` — PDF/TIFF under ``figures/publication/`` (same ``--include-fig6-unique`` adds Fig 6 unique publication exports)

Pass ``--include-fig6-unique`` on **this** script to forward ``--write-fig6-unique-supplement`` / ``--include-fig6-unique-split`` to steps 4–5. **Default:** no Fig 6 unique files are created.

Usage::

    python regenerate_all_figures.py
    python regenerate_all_figures.py --strict
    python regenerate_all_figures.py --skip-catalog --strict
"""
from __future__ import annotations

import argparse
import sys

from orchestrate_subprocess import call_echo
from repo_paths import REPO_ROOT


def run_root(script: str, args: list[str]) -> int:
    cmd = [sys.executable, str(REPO_ROOT / script), *args]
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="Forward --strict to each step and exit non-zero on first failure.")
    ap.add_argument("--skip-catalog", action="store_true", help="Skip generate_catalog_figures.py.")
    ap.add_argument("--skip-netmhc-bundle", action="store_true", help="Skip generate_netmhc_figure_bundle.py.")
    ap.add_argument("--skip-netmhc-fig5-fig6-supplement", action="store_true", help="Skip generate_netmhc_fig5_fig6_supplement.py.")
    ap.add_argument("--skip-canonical", action="store_true", help="Skip generate_canonical_manuscript_figures.py.")
    ap.add_argument("--skip-export", action="store_true", help="Skip export_publication_figures.py.")
    ap.add_argument(
        "--include-fig6-unique",
        action="store_true",
        help="Forward --write-fig6-unique-supplement to generate_canonical_manuscript_figures.py and "
        "--include-fig6-unique-split to export_publication_figures.py (default: no Fig 6 unique anywhere).",
    )
    args = ap.parse_args()

    strict = ["--strict"] if args.strict else []
    fig6_unique: list[str] = []
    if args.include_fig6_unique:
        fig6_unique = ["--write-fig6-unique-supplement"]
    fig6_export: list[str] = []
    if args.include_fig6_unique:
        fig6_export = ["--include-fig6-unique-split"]
    steps: list[tuple[str, list[str]]] = []
    if not args.skip_catalog:
        steps.append(("generate_catalog_figures.py", strict))
    if not args.skip_netmhc_bundle:
        steps.append(("generate_netmhc_figure_bundle.py", strict))
    if not args.skip_netmhc_fig5_fig6_supplement:
        steps.append(("generate_netmhc_fig5_fig6_supplement.py", strict))
    if not args.skip_canonical:
        steps.append(("generate_canonical_manuscript_figures.py", [*strict, *fig6_unique]))
    if not args.skip_export:
        steps.append(("export_publication_figures.py", [*strict, *fig6_export]))

    if not steps:
        print("No steps selected.", file=sys.stderr)
        return 1

    for script, extra in steps:
        code = run_root(script, extra)
        if code != 0:
            print(f"Failed: {script} (exit {code})", file=sys.stderr)
            return int(code)
    print("regenerate_all_figures: all selected steps finished OK.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
