#!/bin/bash
# runVerification.sh
#
# Build-and-run driver for the verification cases in this directory.
#
# Every case in verificationCases/ compares against an exact solution of
# the equations the solver implements and prints a bare PASS or FAIL line
# on stdout, exiting non-zero on failure. This script builds each case in
# its own build-<name>/ directory, runs it, and reports the measured
# errors the case printed.
#
# A case whose source mentions SKIP_EMBED_GUARDS also carries a negative
# control: the same source built with -DSKIP_EMBED_GUARDS must FAIL. The
# driver builds that variant too and treats a passing guard-less build as
# a failure of the suite, because it would mean the guard under test does
# nothing.
#
# Usage:
#   bash verificationCases/runVerification.sh [case.c ...]
#
#   With no arguments, runs every .c in this directory.
#
# Environment:
#   VERIFICATION_THREADS  If set to N > 1, build with -fopenmp and run on N
#                         threads. Unset (the default) builds serial, which
#                         is bitwise reproducible; OpenMP reductions sum in
#                         a nondeterministic order and shift results at the
#                         1e-5 level. Use threads to iterate, serial to
#                         record a result. staticFilmTube is the expensive
#                         case: about 49 min serial, 10 min on 16 threads.
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ -f "${REPO_ROOT}/.project_config" ]]; then
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/.project_config"
fi

if ! command -v qcc >/dev/null 2>&1; then
  echo "ERROR: qcc not found in PATH." >&2
  exit 1
fi

QCC_FLAGS=(-I"${REPO_ROOT}/src-local" -O2 -Wall -disable-dimensions)

# Optional OpenMP. Serial is the default because a verification result
# should be reproducible; threaded reductions are not bitwise stable.
THREADS="${VERIFICATION_THREADS:-1}"
if [[ ! "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "ERROR: VERIFICATION_THREADS must be a positive integer." >&2
  exit 1
fi
if [[ "$THREADS" -gt 1 ]]; then
  QCC_FLAGS+=(-fopenmp)
  export OMP_NUM_THREADS="$THREADS"
  echo "OpenMP build: ${THREADS} threads (results are not bitwise reproducible)"
  echo ""
fi

FAILED=0

# Build one case into its own directory and run it.
# Prints the case's own stdout, indented, on success.
build_and_run() {
  local src="$1" name="$2" build_dir="$3"
  shift 3
  rm -rf "$build_dir"
  mkdir -p "$build_dir"
  cp "$src" "$build_dir/"
  (
    cd "$build_dir"
    qcc "${QCC_FLAGS[@]}" "$@" "$(basename "$src")" -o "$name" -lm
    ./"$name" > out.log 2> err.log
  )
}

run_case() {
  local src="$1"
  local name build_dir
  name="$(basename "$src" .c)"
  build_dir="${SCRIPT_DIR}/build-${name}"

  echo "-----------------------------------------"
  echo "Verification: ${name}"
  echo "-----------------------------------------"

  if build_and_run "$src" "$name" "$build_dir" && \
     grep -q "^PASS$" "${build_dir}/out.log"; then
    echo "${name}: PASS"
    # A case whose only stdout is the PASS line makes grep exit 1, which
    # would abort the driver under `set -e`.
    grep -v "^PASS$" "${build_dir}/out.log" | sed 's/^/  /' || true
  else
    echo "${name}: FAIL" >&2
    sed 's/^/  /' "${build_dir}/out.log" >&2 || true
    tail -5 "${build_dir}/err.log" | sed 's/^/  /' >&2 || true
    FAILED=1
  fi

  # Machine-readable artefacts the case may have written.
  local csv
  while IFS= read -r csv; do
    echo "  raw errors: ${csv#"${REPO_ROOT}/"}"
  done < <(find "$build_dir" -maxdepth 1 -name '*.csv' -print | sort)

  # Negative control, when the case declares one.
  if grep -q "SKIP_EMBED_GUARDS" "$src"; then
    local nc_dir="${SCRIPT_DIR}/build-${name}-noguards"
    echo "  negative control: rebuilding with -DSKIP_EMBED_GUARDS (must FAIL)"
    # Build status, exit status and output are three different things. A
    # qcc failure or a crash is not evidence that the guard assertion
    # fired, so accept the control only when the binary built, exited
    # non-zero, and said FAIL for itself.
    rm -rf "$nc_dir"
    mkdir -p "$nc_dir"
    cp "$src" "$nc_dir/"
    local nc_built=0 nc_status=0
    if ( cd "$nc_dir" && qcc "${QCC_FLAGS[@]}" -DSKIP_EMBED_GUARDS \
           "$(basename "$src")" -o "$name" -lm ) 2>"${nc_dir}/build.log"; then
      nc_built=1
      ( cd "$nc_dir" && ./"$name" > out.log 2> err.log ) || nc_status=$?
    fi

    if [[ $nc_built -eq 0 ]]; then
      echo "${name} negative control: BUILD FAILED, so the control proves nothing" >&2
      tail -5 "${nc_dir}/build.log" | sed 's/^/    /' >&2 || true
      FAILED=1
    elif [[ $nc_status -eq 0 ]]; then
      echo "${name} negative control: UNEXPECTED PASS without the embed guards" >&2
      echo "  the guard under test is not detecting the coupling it claims to" >&2
      FAILED=1
    elif ! grep -q "^FAIL$" "${nc_dir}/out.log"; then
      echo "${name} negative control: exited ${nc_status} without printing FAIL" >&2
      echo "  it crashed rather than failing its own assertion" >&2
      tail -5 "${nc_dir}/err.log" | sed 's/^/    /' >&2 || true
      FAILED=1
    else
      echo "  negative control: FAILED as required (exit ${nc_status})"
      grep -iE "^(P[0-9]|probe|violation)" "${nc_dir}/out.log" | sed 's/^/    /' || true
    fi
  fi

  echo ""
}

if [[ $# -gt 0 ]]; then
  for src in "$@"; do
    run_case "$src"
  done
else
  shopt -s nullglob
  cases=("${SCRIPT_DIR}"/*.c)
  shopt -u nullglob
  if [[ ${#cases[@]} -eq 0 ]]; then
    echo "ERROR: no verification cases found in ${SCRIPT_DIR}." >&2
    exit 1
  fi
  for src in "${cases[@]}"; do
    run_case "$src"
  done
fi

echo "========================================="
if [[ $FAILED -eq 0 ]]; then
  echo "All verification cases passed."
else
  echo "Some verification cases FAILED." >&2
fi
echo "========================================="
exit $FAILED
