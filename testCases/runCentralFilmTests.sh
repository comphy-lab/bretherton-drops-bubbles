#!/bin/bash
# Unit tests for the solver-independent central-film observation contract.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/bretherton-central-film-XXXXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
CC_BIN="${CC:-cc}"

"${CC_BIN}" \
  -std=c99 -Wall -Wextra -Werror -pedantic \
  "${SCRIPT_DIR}/centralFilmObserver.c" -lm \
  -o "${WORK_DIR}/centralFilmObserver"

"${WORK_DIR}/centralFilmObserver"
