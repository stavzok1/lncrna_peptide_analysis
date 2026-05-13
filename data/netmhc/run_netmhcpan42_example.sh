#!/usr/bin/env bash
# Example NetMHCpan-4.2 run from repo root (WSL paths).
#
# NetMHCpan requires the install root in env var NETMHCpan (note spelling), e.g.
#   /home/you/netMHCpan-4.2
# so it can open data/version, data/MHC_pseudo.dat, etc. Without it you get:
#   Unable to open(r) file $NETMHCpan/data/version
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

export NETMHCpan="${NETMHCpan:-$HOME/netMHCpan-4.2}"
NETMHCPAN="${NETMHCPAN:-$NETMHCpan/Linux_x86_64/bin/netMHCpan-4.2}"
ALLELES="${ALLELES:-data/netmhc/hla_european27_class1.txt}"

# The main binary lives under Linux_x86_64/bin, but the program shells out to
# $NETMHCpan/bin/nnalign_* and $NETMHCpan/bin/estimate_PCC. If bin/ is missing,
# you get: sh: 1: .../netMHCpan-4.2/bin/nnalign_...: not found
LINUX_BIN="$NETMHCpan/Linux_x86_64/bin"
if [[ -d "$LINUX_BIN" && ! -d "$NETMHCpan/bin" ]]; then
  ln -sfn "$LINUX_BIN" "$NETMHCpan/bin"
fi

# NetMHCpan builds work dirs under $TMPDIR (see -h: -tdir). On WSL, an empty or
# Windows-backed TMPDIR while cwd is /mnt/c/... often triggers: "Cannot make tmpdir. Exit"
export TMPDIR="${TMPDIR:-$HOME/tmp}"
mkdir -p "$TMPDIR"

# NetMHCpan-4.2 -a expects comma-separated allele names (see -h: "multiple alleles with ,").
# Passing a filesystem path (with or without a "file:" prefix) is treated as a single allele
# ID, which triggers: "... cannot be found in hla_pseudo list".
ALLELE_ARG="$(tr -d '\r' < "$ALLELES" | grep -v '^[[:space:]]*$' | paste -sd, -)"

# Full optional outputs in 4.2: BA + pathogen (IEDB) + neo (CEDAR). Add -context 1 if you want context encoding.
FLAGS=( -inptype 0 -a "$ALLELE_ARG" -l 9 -BA 1 -pathogen 1 -neo 1 -context 0 -xls 1 )

"$NETMHCPAN" -f data/netmhc/ninemers_sig_lnc.fasta "${FLAGS[@]}" -xlsfile data/netmhc/netmhcpan_sig_lnc.xls
"$NETMHCPAN" -f data/netmhc/ninemers_coding_control.fasta "${FLAGS[@]}" -xlsfile data/netmhc/netmhcpan_coding_control.xls

echo "Done. Outputs: data/netmhc/netmhcpan_sig_lnc.xls data/netmhc/netmhcpan_coding_control.xls"
