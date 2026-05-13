#!/usr/bin/env python3
"""
Render Fig 6 split panels (A–E) from the **NetMHCpan-4.1** TTN-AS1 XLS in::

    data/netmhc/predictions/ttn_as1_smpep108065_netmhc41/netmhcpan_ttn_as1_108065_netmhc41.xls

Writes under::

    figures/manuscript_netmhc/fig6/ttn_netmhc41_ba_el/

Override paths with ``--netmhc-xls`` / ``--out-dir``. This script only orchestrates
``plot_figure6_ttn_as1_allele_coverage.py`` (NetMHC ``--gating netmhc`` by default).
"""
from __future__ import annotations

from pathlib import Path
import sys

_REPO = Path(__file__).resolve().parent.parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from repo_paths import DATA, FIGURES, REPO_ROOT, NETMHC_DATA, NETMHC_FIGURES, MANUSCRIPT_DIR

ROOT = REPO_ROOT


import argparse
import subprocess

DEFAULT_XLS = NETMHC_DATA / "predictions/ttn_as1_smpep108065_netmhc41/netmhcpan_ttn_as1_108065_netmhc41.xls"
DEFAULT_OUT_DIR = FIGURES / "manuscript_netmhc/fig6/ttn_netmhc41_ba_el"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--netmhc-xls",
        type=Path,
        default=DEFAULT_XLS,
        help="Wide XLS from NetMHCpan-4.1 (-xls 1 -BA 1).",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Directory for fig6_ttn_split_*.png (created if missing).",
    )
    ap.add_argument(
        "--no-also-write-unique",
        action="store_true",
        help="Do not emit *_unique* split panels (default: emit them).",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the command instead of running.",
    )
    args = ap.parse_args()

    xls: Path = args.netmhc_xls
    if not xls.is_file():
        raise SystemExit(f"Missing XLS (run the 4.1 shell script first): {xls}")

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "fig6_ttn_split.png"

    also_unique = not bool(args.no_also_write_unique)

    cmd = [
        sys.executable,
        str(MANUSCRIPT_DIR / "plot_figure6_ttn_as1_allele_coverage.py"),
        "--gating",
        "netmhc",
        "--netmhc-xls",
        str(xls),
        "--split-panels",
        "-o",
        str(out_png),
    ]
    if also_unique:
        cmd.append("--also-write-unique")

    print("+", " ".join(cmd))
    if args.dry_run:
        return
    r = subprocess.run(cmd, cwd=str(ROOT))
    if r.returncode != 0:
        raise SystemExit(r.returncode)
    print(f"Wrote split panels under {out_dir}")


if __name__ == "__main__":
    main()
