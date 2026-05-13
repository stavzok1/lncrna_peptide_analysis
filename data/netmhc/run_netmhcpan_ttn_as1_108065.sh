#!/usr/bin/env bash
# NetMHCpan-4.2 for TTN-AS1 sliding 9-mers (smPEP_ID 108065 / SPROHSA108065).
# For NetMHCpan-4.1 + separate output folder, use:
#   data/netmhc/predictions/ttn_as1_smpep108065_netmhc41/run_netmhcpan41_ttn_as1_108065.sh
# Run from repo root in WSL: bash data/netmhc/run_netmhcpan_ttn_as1_108065.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

export NETMHCpan="${NETMHCpan:-$HOME/netMHCpan-4.2}"
NETMHCPAN="${NETMHCPAN:-$NETMHCpan/Linux_x86_64/bin/netMHCpan-4.2}"
ALLELES="${ALLELES:-data/netmhc/hla_european27_class1.txt}"

LINUX_BIN="$NETMHCpan/Linux_x86_64/bin"
if [[ -d "$LINUX_BIN" && ! -d "$NETMHCpan/bin" ]]; then
  ln -sfn "$LINUX_BIN" "$NETMHCpan/bin"
fi

export TMPDIR="${TMPDIR:-$HOME/tmp}"
mkdir -p "$TMPDIR"

ALLELE_ARG="$(tr -d '\r' < "$ALLELES" | grep -v '^[[:space:]]*$' | paste -sd, -)"
FLAGS=( -inptype 0 -a "$ALLELE_ARG" -l 9 -BA 1 -pathogen 1 -neo 1 -context 0 -xls 1 )

"$NETMHCPAN" -f data/netmhc/ttn_as1_108065_ninemers.fasta "${FLAGS[@]}" \
  -xlsfile data/netmhc/netmhcpan_ttn_as1_108065.xls

echo "Done. Output: data/netmhc/netmhcpan_ttn_as1_108065.xls"
