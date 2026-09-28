#!/usr/bin/env bash
set -euo pipefail

bash qc/r450_rscript.sh -e 'stopifnot(isTRUE(getOption("r45.runtime.active")), identical(.libPaths(), Sys.getenv("R45_EXPECTED_LIBRARY"))); library(cli); cat("PASS\n")' >/dev/null

if bash qc/r450_rscript.sh -e 'library(cli); stop("intentional failure")' >/dev/null 2>&1; then
  printf 'R failure was hidden by the launcher\n' >&2
  exit 1
else
  status=$?
  if [[ "$status" -ne 1 ]]; then
    printf 'Expected R failure exit 1, got %s\n' "$status" >&2
    exit 1
  fi
fi

if bash qc/r450_rscript.sh --vanilla -e '1+1' >/dev/null 2>&1; then
  printf 'Launcher unexpectedly accepted --vanilla\n' >&2
  exit 1
else
  status=$?
  if [[ "$status" -ne 2 ]]; then
    printf 'Expected launcher rejection exit 2, got %s\n' "$status" >&2
    exit 1
  fi
fi

printf 'PASS: R 4.5.0 launcher preserves success and failure exit codes.\n'
