#!/usr/bin/env bash
set -euo pipefail

rscript=/d/R/R-4.5.0/bin/Rscript.exe
if [[ ! -f "$rscript" ]]; then
  printf 'Required Rscript missing: %s\n' "$rscript" >&2
  exit 2
fi
if [[ $# -eq 0 ]]; then
  printf 'Usage: bash qc/r450_rscript.sh SCRIPT.R [args...]\n' >&2
  exit 2
fi
for arg in "$@"; do
  if [[ "$arg" == --vanilla || "$arg" == --no-init-file ]]; then
    printf '%s disables the required R exit profile\n' "$arg" >&2
    exit 2
  fi
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export R_PROFILE_USER="$(cygpath -m "$script_dir/r450_exit_profile.R")"
export CLI_NO_THREAD=1
exec "$rscript" "$@"
