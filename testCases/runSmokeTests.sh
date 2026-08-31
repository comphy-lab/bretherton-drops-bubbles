#!/bin/bash
# runSmokeTests.sh
#
# Smoke test for the production solver.
#
# A smoke test proves only that simulationCases/bretherton.c compiles and
# integrates a few dozen time steps at coarse resolution without blowing
# up. It compares against nothing and therefore supports no claim about
# the discretisation or the physics. The exact-solution comparisons live
# in verificationCases/ and are driven by
# verificationCases/runVerification.sh.
#
# Usage:
#   bash testCases/runSmokeTests.sh
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

for arg in "$@"; do
  case "$arg" in
    -h|--help)
      sed -n '2,15p' "${BASH_SOURCE[0]}"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $arg" >&2
      exit 1
      ;;
  esac
done

if [[ -f "${REPO_ROOT}/.project_config" ]]; then
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/.project_config"
fi

if ! command -v qcc >/dev/null 2>&1; then
  echo "ERROR: qcc not found in PATH." >&2
  exit 1
fi

echo "-----------------------------------------"
echo "Smoke test: simulationCases/bretherton.c"
echo "-----------------------------------------"

SMOKE_CASE_DIR="${REPO_ROOT}/simulationCases/9999"
rm -rf "$SMOKE_CASE_DIR"

FAILED=0
if bash "${REPO_ROOT}/runSimulation.sh" "testCases/smoke.params"; then
  SMOKE_LOG="${SMOKE_CASE_DIR}/c9999-log"
  if [[ -f "$SMOKE_LOG" ]] && \
     [[ "$(grep -cv '^#' "$SMOKE_LOG")" -ge 10 ]] && \
     ! grep -qi "blew up" "$SMOKE_LOG"; then
    echo "bretherton.c smoke test: PASS ($(grep -cv '^#' "$SMOKE_LOG") logged steps)"
  else
    echo "bretherton.c smoke test: FAIL (missing or truncated log)" >&2
    FAILED=1
  fi
else
  echo "bretherton.c smoke test: FAIL (runner exited non-zero)" >&2
  FAILED=1
fi

exit $FAILED
