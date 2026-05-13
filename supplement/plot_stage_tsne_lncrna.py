"""
Deprecated: use ``manuscript/plot_figure1b_tsne_stage_lncrna.py`` (Figure 1B).

This shim forwards argv unchanged so old commands keep working.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
_MS = _REPO / "manuscript" / "plot_figure1b_tsne_stage_lncrna.py"


def main() -> None:
    cmd = [sys.executable, str(_MS), *sys.argv[1:]]
    print(f"[deprecated] {Path(__file__).name} → {_MS.name}", flush=True)
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
