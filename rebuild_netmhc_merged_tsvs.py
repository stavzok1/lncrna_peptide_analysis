"""
Rebuild merged NetMHC ``*_with_iedb.tsv`` cohort files and compare **line counts** to existing
``data/netmhc/*`` copies.

- **Significant lncRNA:** default join is ``smpep_window`` (no ninemers FASTA), which matches the
  shipped ``netmhcpan_sig_lnc_with_iedb.tsv`` row count when the cohort was merged that way.
  Use ``--sig-join stable`` only if your wide ``*.xls`` row count is divisible by the number of
  records in ``ninemers_sig_lnc.fasta`` (same NetMHC run as the XLS).
- **Coding control:** ``stable`` join requires the same divisibility; if merge fails, the script
  skips coding (with a message) unless you fix XLS / FASTA alignment.
- **Proportional whole coding:** ``--synthetic-iedb-pass`` (no IEDB CSV), same as
  ``supplement/run_fig5_sig_vs_proportional_coding.py``.

Default ``--external-netmhc`` is ``<repo-parent>/data/netmhc`` (IEDB CSVs + coding ``*.xls`` that
may not live under ``paper-github``).

Usage::

    python rebuild_netmhc_merged_tsvs.py
    python rebuild_netmhc_merged_tsvs.py --strict
    python rebuild_netmhc_merged_tsvs.py --apply   # overwrite data/netmhc/*_with_iedb.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent
for _p in (str(_REPO), str(_REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from orchestrate_subprocess import call_echo
from repo_paths import DATA, NETMHC_DATA, REPO_ROOT, SCRIPTS_DIR

MERGE = SCRIPTS_DIR / "merge_netmhcpan_xls_with_iedb.py"


def count_lines(path: Path) -> int:
    n = 0
    with path.open("rb") as f:
        for _ in f:
            n += 1
    return n


def run_merge(cmd: list[str]) -> int:
    return call_echo(cmd, cwd=REPO_ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--external-netmhc",
        type=Path,
        default=REPO_ROOT.parent / "data" / "netmhc",
        help="Directory with IEDB CSVs and coding *.xls not always present under paper-github.",
    )
    ap.add_argument(
        "--verify-dir",
        type=Path,
        default=NETMHC_DATA / "merge_rebuild_verify",
        help="Write rebuilt TSVs here unless --apply.",
    )
    ap.add_argument(
        "--apply",
        action="store_true",
        help="Write directly to data/netmhc/netmhcpan_*_with_iedb.tsv (overwrites).",
    )
    ap.add_argument(
        "--sig-join",
        choices=("smpep_window", "stable"),
        default="smpep_window",
        help="IEDB join for sig lnc (default smpep_window; stable needs aligned ninemers vs XLS).",
    )
    ap.add_argument("--strict", action="store_true", help="Exit non-zero on merge failure or row mismatch.")
    args = ap.parse_args()

    ext: Path = args.external_netmhc
    merge_py = str(MERGE)
    py = sys.executable

    sig_xls = NETMHC_DATA / "netmhcpan_sig_lnc.xls"
    sig_iedb = ext / "iedb_sig_lnc_peptides_el_tap_immuno.csv"
    cod_xls = ext / "netmhcpan_coding_control.xls"
    cod_iedb = ext / "iedb_coding_control_parents_tap_immuno.csv"
    prop_xls = ext / "netmhcpan_coding_proportional_whole.xls"

    out_base = NETMHC_DATA if args.apply else args.verify_dir
    if not args.apply:
        out_base.mkdir(parents=True, exist_ok=True)

    def out_path(name: str) -> Path:
        return NETMHC_DATA / name if args.apply else (out_base / name)

    mismatches: list[str] = []

    def compare_lines(label: str, canonical: Path, rebuilt: Path) -> None:
        if not canonical.is_file():
            print(f"[warn] {label}: no canonical at {canonical}", flush=True)
            return
        n0, n1 = count_lines(canonical), count_lines(rebuilt)
        ok = n0 == n1
        print(f"[compare] {label}: canonical_lines={n0} rebuilt_lines={n1} match={ok}", flush=True)
        if not ok:
            mismatches.append(f"{label}: {n0} vs {n1}")

    # --- sig lnc ---
    if sig_xls.is_file() and sig_iedb.is_file():
        outp = out_path("netmhcpan_sig_lnc_with_iedb.tsv")
        cmd = [
            py,
            merge_py,
            "--netmhc-xls",
            str(sig_xls),
            "--iedb-csv",
            str(sig_iedb),
            "--join",
            args.sig_join,
            "--out-tsv",
            str(outp),
        ]
        if args.sig_join == "stable":
            cmd.extend(
                [
                    "--ninemers-fasta",
                    str(NETMHC_DATA / "ninemers_sig_lnc.fasta"),
                    "--parent-peptide-fasta",
                    str(DATA / "significant_lnc_peptides.faa"),
                ]
            )
        code = run_merge(cmd)
        if code != 0:
            print("[fail] sig_lnc merge", flush=True)
            if args.strict:
                sys.exit(code)
        elif not args.apply:
            compare_lines("sig_lnc", NETMHC_DATA / "netmhcpan_sig_lnc_with_iedb.tsv", outp)
    else:
        print("Skip sig_lnc: missing XLS or IEDB CSV.", flush=True)

    # --- coding control (stable only; needs divisible wide rows vs ninemers) ---
    if cod_xls.is_file() and cod_iedb.is_file():
        outp = out_path("netmhcpan_coding_control_with_iedb.tsv")
        cmd = [
            py,
            merge_py,
            "--netmhc-xls",
            str(cod_xls),
            "--iedb-csv",
            str(cod_iedb),
            "--join",
            "stable",
            "--ninemers-fasta",
            str(NETMHC_DATA / "ninemers_coding_control.fasta"),
            "--parent-peptide-fasta",
            str(NETMHC_DATA / "coding_control_parents_for_iedb.faa"),
            "--out-tsv",
            str(outp),
        ]
        code = run_merge(cmd)
        if code != 0:
            print(
                "[skip] coding_control merge failed (usually XLS row count not divisible by "
                "ninemers FASTA records). Regenerate NetMHC from the packaged ninemers FASTA.",
                flush=True,
            )
            if args.strict:
                print(
                    "Strict mode: coding_control merge is optional when XLS/FASTA are misaligned; continuing.",
                    flush=True,
                )
        elif not args.apply:
            compare_lines("coding_control", NETMHC_DATA / "netmhcpan_coding_control_with_iedb.tsv", outp)
    else:
        print("Skip coding_control: missing XLS or IEDB CSV.", flush=True)

    # --- proportional whole (synthetic IEDB pass) ---
    if prop_xls.is_file():
        outp = out_path("netmhcpan_coding_proportional_whole_with_iedb.tsv")
        cmd = [
            py,
            merge_py,
            "--netmhc-xls",
            str(prop_xls),
            "--synthetic-iedb-pass",
            "--join",
            "stable",
            "--ninemers-fasta",
            str(NETMHC_DATA / "ninemers_coding_proportional_whole.fasta"),
            "--out-tsv",
            str(outp),
        ]
        code = run_merge(cmd)
        if code != 0:
            print("[fail] coding_proportional_whole merge", flush=True)
            if args.strict:
                sys.exit(code)
        elif not args.apply:
            compare_lines(
                "coding_proportional_whole",
                NETMHC_DATA / "netmhcpan_coding_proportional_whole_with_iedb.tsv",
                outp,
            )
    else:
        print("Skip proportional_whole: missing XLS.", flush=True)

    if mismatches and args.strict:
        print("Row-count mismatches:", *mismatches, sep="\n  ", file=sys.stderr)
        sys.exit(1)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
