#!/usr/bin/env python3
"""
Build a minimal IEDB-style **peptide_table** CSV for TTN-AS1 (smPEP 108065) so
``plot_figure6_ttn_as1_allele_coverage.py --gating iedb_sb`` can merge on ``stable_key``.

The cohort ``iedb_sig_lnc_peptides_el_tap_immuno.csv`` does not contain TTN ninemers; this
companion carries **synthetic** immunogenicity / processing scores that satisfy default Fig 5
IEDB gates, and **netmhcpan_ba_ic50** from the same BA scores as the wide XLS (so IC50 gating
tracks NetMHCpan BA).

Outputs columns compatible with ``plot_figure6_ttn_as1_allele_coverage.merge_iedb_csv`` /
``netmhc_sb_core.sb_mask_spec`` (expects ``score`` → ``iedb_score``, ``processing_score`` →
``iedb_processing_score`` after merge rename).

Example::

    python scripts/build_ttn_iedb_companion_csv.py \\
        --out data/netmhc/ttn_as1_iedb_companion_synthetic.csv \\
        --parent-input-seq-id '108065|TTN-AS1|synthetic'
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


import argparse
import math
import re

import pandas as pd

_SCRIPTS = ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from plot_figure6_ttn_as1_allele_coverage import (  # noqa: E402
    build_ttn_long_for_iedb_merge,
    normalize_hla_netmhc_to_iedb,
    parse_wide_netmhc_xls_rows,
)


def ba_score_to_ic50_nm(ba: float) -> float:
    x = float(ba)
    if not math.isfinite(x):
        return float("nan")
    x = min(1.0 - 1e-12, max(1e-12, x))
    return float(50000.0 ** (1.0 - x))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--netmhc-xls",
        type=Path,
        default=ROOT / "data" / "netmhc" / "netmhcpan_ttn_as1_108065.xls",
    )
    ap.add_argument(
        "--parent-input-seq-id",
        type=str,
        default="108065|TTN-AS1|synthetic",
        help="Must match ``--iedb-parent-input-seq-id`` for Fig 6 merge.",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "data" / "netmhc" / "ttn_as1_iedb_companion_synthetic.csv",
    )
    ap.add_argument(
        "--imm-score",
        type=float,
        default=0.25,
        help="Synthetic IEDB ``score`` (immunogenicity proxy; > default imm cut).",
    )
    ap.add_argument(
        "--proc-score",
        type=float,
        default=2.0,
        help="Synthetic ``processing_score`` (> default proc cut).",
    )
    args = ap.parse_args()

    alleles, starts, peps, ba, ba_rank, el_rank = parse_wide_netmhc_xls_rows(args.netmhc_xls)
    long_df = build_ttn_long_for_iedb_merge(
        starts, peps, alleles, ba, ba_rank, el_rank, str(args.parent_input_seq_id).strip()
    )

    rows: list[dict[str, object]] = []
    for i in range(len(long_df)):
        r = long_df.iloc[i]
        si, ai = int(r["_si"]), int(r["_ai"])
        pep = str(r["Peptide"]).upper()
        elp = float(el_rank[si, ai])
        ic50 = ba_score_to_ic50_nm(float(ba[si, ai]))
        sk = str(r["stable_key"])
        al_disp = normalize_hla_netmhc_to_iedb(str(alleles[ai]))
        pos0 = int(starts[si])
        rows.append(
            {
                "sequence_number": i + 1,
                "peptide": pep,
                "start": pos0 + 1,
                "end": pos0 + 9,
                "length": 9,
                "allele": al_disp,
                "peptide_index": si + 1,
                "netmhcpan_el_percentile": elp,
                "netmhcpan_ba_ic50": ic50,
                "score": float(args.imm_score),
                "processing_score": float(args.proc_score),
                "iedb_chunk_index": 1,
                "input_seq_id": str(args.parent_input_seq_id).strip(),
                "stable_key": sk,
            }
        )

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
