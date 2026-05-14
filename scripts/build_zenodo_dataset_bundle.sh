#!/usr/bin/env bash
# Run the Zenodo bundler from Git Bash or WSL (delegates to Windows PowerShell).
# Usage:
#   ./scripts/build_zenodo_dataset_bundle.sh "/c/Users/you/Desktop/masters/UNDEFINED" --zip
#   ./scripts/build_zenodo_dataset_bundle.sh /mnt/c/Users/you/Desktop/masters/UNDEFINED --zip
set -euo pipefail

usage() {
  echo "Usage: $0 <DataParent> [--zip]" >&2
  echo "  DataParent = directory that contains a 'data' subfolder (e.g. UNDEFINED)." >&2
  echo "  Git Bash:  $0 /c/Users/you/.../UNDEFINED --zip" >&2
  echo "  WSL:       $0 /mnt/c/Users/you/.../UNDEFINED --zip" >&2
  exit 1
}

[[ $# -ge 1 ]] || usage
DATA_PARENT="$1"
ZIP_ARGS=()
[[ "${2:-}" == "--zip" ]] && ZIP_ARGS=(-Zip)
[[ -n "${2:-}" && "${2:-}" != "--zip" ]] && usage

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PS_SCRIPT="$REPO_ROOT/scripts/build_zenodo_dataset_bundle.ps1"

find_powershell() {
  if command -v powershell.exe >/dev/null 2>&1; then
    command -v powershell.exe
    return
  fi
  local wsl_ps="/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
  if [[ -x "$wsl_ps" ]]; then
    echo "$wsl_ps"
    return
  fi
  return 1
}

to_windows_path() {
  local p="$1"
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -aw "$p"
  elif command -v wslpath >/dev/null 2>&1; then
    wslpath -w "$p"
  else
    # Last resort: hope PowerShell accepts this path
    printf '%s' "$p"
  fi
}

PS_BIN="$(find_powershell)" || {
  echo "error: powershell.exe not found. Use Windows PowerShell or cmd instead:" >&2
  echo '  cd paper-github' >&2
  echo '  .\scripts\build_zenodo_dataset_bundle.ps1 -DataParent "C:\path\to\UNDEFINED" -Zip' >&2
  exit 1
}

DATA_WIN="$(to_windows_path "$DATA_PARENT")"
SCRIPT_WIN="$(to_windows_path "$PS_SCRIPT")"

exec "$PS_BIN" -NoProfile -ExecutionPolicy Bypass -File "$SCRIPT_WIN" -DataParent "$DATA_WIN" "${ZIP_ARGS[@]}"
