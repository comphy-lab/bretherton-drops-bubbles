#!/bin/bash
# runSmokeTests.sh
#
# Compile-and-run checks for the Basilisk cases in this repository.
# Basilisk "tests" are not unit tests: a smoke test proves that a case
# compiles and integrates a few time steps; the two verification cases
# additionally compare against exact solutions and print PASS/FAIL.
#
# Usage:
#   bash testCases/runSmokeTests.sh [--skip-main]
#
#   --skip-main   Skip the (slower) bretherton.c smoke run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SKIP_MAIN=0
for arg in "$@"; do
  case "$arg" in
    --skip-main) SKIP_MAIN=1 ;;
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

if [[ -f "${REPO_ROOT}/.project_config" ]]; then
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/.project_config"
fi

if ! command -v qcc >/dev/null 2>&1; then
  echo "ERROR: qcc not found in PATH." >&2
  exit 1
fi

FAILED=0

run_verification() {
  local src="$1"
  local name
  name="$(basename "$src" .c)"
  local build_dir="${SCRIPT_DIR}/build-${name}"

  echo "-----------------------------------------"
  echo "Verification: ${name}"
  echo "-----------------------------------------"
  rm -rf "$build_dir"
  mkdir -p "$build_dir"
  cp "$src" "$build_dir/"
  (
    cd "$build_dir"
    qcc -I"${REPO_ROOT}/src-local" -O2 -Wall -disable-dimensions \
        "$(basename "$src")" -o "$name" -lm
    ./"$name" > out.log 2> err.log
  )
  if grep -q "^PASS$" "${build_dir}/out.log"; then
    echo "${name}: PASS"
    grep -v "^PASS$" "${build_dir}/out.log" | sed 's/^/  /'
  else
    echo "${name}: FAIL" >&2
    sed 's/^/  /' "${build_dir}/out.log" >&2 || true
    tail -5 "${build_dir}/err.log" | sed 's/^/  /' >&2 || true
    FAILED=1
  fi
  echo ""
}

run_verification "${SCRIPT_DIR}/embedAxiVofAdvection.c"
run_verification "${SCRIPT_DIR}/laplaceEmbedTube.c"

if [[ $SKIP_MAIN -eq 0 ]]; then
  echo "-----------------------------------------"
  echo "Smoke test: simulationCases/bretherton.c"
  echo "-----------------------------------------"
  SMOKE_CASE_DIR="${REPO_ROOT}/simulationCases/9999"
  rm -rf "$SMOKE_CASE_DIR"
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
  echo ""
fi

echo "========================================="
if [[ $FAILED -eq 0 ]]; then
  echo "All smoke/verification tests passed."
else
  echo "Some tests FAILED." >&2
fi
echo "========================================="
exit $FAILED
