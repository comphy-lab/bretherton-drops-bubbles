#!/bin/bash
# runTests.sh
#
# Single entry point for the repository's evidence suite. Runs the
# pure-C unit tests, verification cases (exact-solution comparisons), and then
# the smoke test (compiles and integrates a few dozen steps, compares against
# nothing).
#
# Usage:
#   bash runTests.sh                    # unit, verification, then smoke tests
#   bash runTests.sh --unit             # central-film observer unit tests only
#   bash runTests.sh --verification     # verification cases only
#   bash runTests.sh --smoke            # smoke test only
#
# Environment:
#   VERIFICATION_THREADS=N  Build the verification cases with -fopenmp and
#                           run them on N threads. Serial by default, because
#                           a recorded verification result should be bitwise
#                           reproducible. The full serial suite is dominated
#                           by staticFilmTube (about 49 min; 10 min on 16
#                           threads).
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUN_UNIT=1
RUN_VERIFICATION=1
RUN_SMOKE=1
SELECTED_MODE=""
for arg in "$@"; do
  case "$arg" in
    --unit|--verification|--smoke)
      if [[ -n "${SELECTED_MODE}" ]]; then
        echo "ERROR: suite selectors are mutually exclusive." >&2
        exit 1
      fi
      SELECTED_MODE="${arg}"
      ;;
    -h|--help)
      sed -n '2,12p' "${BASH_SOURCE[0]}"
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $arg" >&2
      exit 1
      ;;
  esac
done

case "${SELECTED_MODE}" in
  --unit)         RUN_VERIFICATION=0; RUN_SMOKE=0 ;;
  --verification) RUN_UNIT=0; RUN_SMOKE=0 ;;
  --smoke)        RUN_UNIT=0; RUN_VERIFICATION=0 ;;
esac

FAILED=0

if [[ $RUN_UNIT -eq 1 ]]; then
  bash "${REPO_ROOT}/testCases/runCentralFilmTests.sh" || FAILED=1
fi

if [[ $RUN_VERIFICATION -eq 1 ]]; then
  bash "${REPO_ROOT}/verificationCases/runVerification.sh" || FAILED=1
fi

if [[ $RUN_SMOKE -eq 1 ]]; then
  bash "${REPO_ROOT}/testCases/runSmokeTests.sh" || FAILED=1
  echo ""
fi

echo "========================================="
if [[ $FAILED -eq 0 ]]; then
  echo "All evidence cases passed."
else
  echo "Some evidence cases FAILED." >&2
fi
echo "========================================="
exit $FAILED
