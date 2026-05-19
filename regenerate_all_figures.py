#!/usr/bin/env python3
"""
Regenerate **all** manuscript and supplementary figures from the repository root.

Does **not** rebuild upstream ``data/`` tables (expression matrices, NetMHC merges, FASTAs).

Default order:

1. ``generate_canonical_manuscript_figures.py`` — main-text panels under ``figures/``
2. ``generate_supplementary_figures.py`` — all panels under ``figures/supplementary/``
3. ``export_publication_figures.py`` — optional PDF/TIFF under ``figures/publication/`` (gitignored)

Pass ``--include-fig6-unique`` to also write Fig 6 unique supplement panels (step 2) and
publication exports for them (step 3). **Default:** no Fig 6 unique files.

Usage::

    python regenerate_all_figures.py
    python regenerate_all_figures.py --strict
    python regenerate_all_figures.py --skip-export --strict
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
    ap.add_argument("--skip-canonical", action="store_true", help="Skip generate_canonical_manuscript_figures.py.")
    ap.add_argument("--skip-supplementary", action="store_true", help="Skip generate_supplementary_figures.py.")
    ap.add_argument("--skip-export", action="store_true", help="Skip export_publication_figures.py.")
    ap.add_argument(
        "--include-fig6-unique",
        action="store_true",
        help="Forward --include-fig6-unique to generate_supplementary_figures.py and "
        "--include-fig6-unique-split to export_publication_figures.py.",
    )
    args = ap.parse_args()

    strict = ["--strict"] if args.strict else []
    fig6_unique: list[str] = ["--include-fig6-unique"] if args.include_fig6_unique else []
    fig6_export: list[str] = ["--include-fig6-unique-split"] if args.include_fig6_unique else []

    steps: list[tuple[str, list[str]]] = []
    if not args.skip_canonical:
        steps.append(("generate_canonical_manuscript_figures.py", strict))
    if not args.skip_supplementary:
        steps.append(("generate_supplementary_figures.py", [*strict, *fig6_unique]))
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
