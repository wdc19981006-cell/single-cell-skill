#!/usr/bin/env bash
set -euo pipefail

# DEPRECATED COMPATIBILITY ENTRYPOINT. New code invokes run_r45.py directly.
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
exec /c/Python312/python.exe "$(cygpath -m "$root/runtime/r45/run_r45.py")" "$@"
