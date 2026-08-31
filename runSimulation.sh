#!/bin/bash
# runSimulation.sh
#
# Run a single bretherton-drops-bubbles simulation from the repository root.
# The script creates simulationCases/<CaseNo>/, copies the parameter file and
# source file, compiles the selected case, and runs it.
#
# Usage:
#   bash runSimulation.sh [params_file] [--exec exec_code] [--threads N]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage: bash runSimulation.sh [params_file] [OPTIONS]

Arguments:
  params_file    Parameter file path (default: default.params)

Options:
  --exec FILE    C source in simulationCases/ (default: bretherton.c)
  --threads N    Thread count; N=1 runs serial (default: 1)
  -h, --help     Show this help message
EOF
}

get_param_value() {
  local key="$1"
  local file="$2"
  awk -F '=' -v key="$key" '
    /^[[:space:]]*#/ { next }
    {
      k = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", k)
      if (k == key) {
        v = $2
        sub(/[[:space:]]*#.*/, "", v)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
        print v
        exit
      }
    }
  ' "$file"
}

# Defaults
EXEC_CODE="bretherton.c"
PARAM_FILE="default.params"
PARAM_FILE_SET=0
OMP_THREADS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --exec)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: --exec requires a file name." >&2
        usage
        exit 1
      fi
      EXEC_CODE="$2"
      shift 2
      ;;
    --exec=*)
      EXEC_CODE="${1#*=}"
      shift
      ;;
    --threads)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: $1 requires a positive integer value." >&2
        usage
        exit 1
      fi
      OMP_THREADS="$2"
      shift 2
      ;;
    --threads=*)
      OMP_THREADS="${1#*=}"
      shift
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      if [[ $PARAM_FILE_SET -eq 0 ]]; then
        PARAM_FILE="$1"
        PARAM_FILE_SET=1
        shift
      else
        echo "ERROR: Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      ;;
  esac
done

if [[ $# -gt 0 ]]; then
  echo "ERROR: Unexpected trailing arguments: $*" >&2
  usage
  exit 1
fi

if [[ ! "$OMP_THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --threads must be a positive integer, got: $OMP_THREADS" >&2
  exit 1
fi

USE_OPENMP=0
if [[ "$OMP_THREADS" -gt 1 ]]; then
  USE_OPENMP=1
fi

if [[ "$EXEC_CODE" != *.c ]]; then
  EXEC_CODE="${EXEC_CODE}.c"
fi

if [[ ! "$PARAM_FILE" = /* ]]; then
  PARAM_FILE="${SCRIPT_DIR}/${PARAM_FILE}"
fi

if [[ -f "${SCRIPT_DIR}/.project_config" ]]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/.project_config"
fi

if ! command -v qcc >/dev/null 2>&1; then
  echo "ERROR: qcc not found in PATH." >&2
  echo "Hint: source your Basilisk environment or provide .project_config." >&2
  exit 1
fi

if [[ ! -f "$PARAM_FILE" ]]; then
  echo "ERROR: Parameter file not found: $PARAM_FILE" >&2
  exit 1
fi

SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${EXEC_CODE}"
if [[ ! -f "$SRC_FILE_ORIG" ]]; then
  echo "ERROR: Source file not found: $SRC_FILE_ORIG" >&2
  exit 1
fi

CASE_NO="$(get_param_value "CaseNo" "$PARAM_FILE")"
if [[ -z "$CASE_NO" ]]; then
  echo "ERROR: CaseNo not found in parameter file: $PARAM_FILE" >&2
  exit 1
fi

if [[ ! "$CASE_NO" =~ ^[0-9]+$ ]]; then
  echo "ERROR: CaseNo must be numeric, got: $CASE_NO" >&2
  exit 1
fi

if [[ "$CASE_NO" -lt 1000 ]]; then
  echo "ERROR: CaseNo must be >= 1000 for consistent sorting, got: $CASE_NO" >&2
  exit 1
fi

CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_NO}"
SRC_FILE_LOCAL="${EXEC_CODE}"
EXECUTABLE_NAME="${EXEC_CODE%.c}"
CASE_LOG_FILE="c${CASE_NO}-log"

echo "========================================="
echo "bretherton-drops-bubbles - Case Runner"
echo "========================================="
echo "Source file: ${EXEC_CODE}"
echo "Parameter file: ${PARAM_FILE}"
echo "CaseNo: ${CASE_NO}"
echo "Case directory: ${CASE_DIR}"
if [[ $USE_OPENMP -eq 1 ]]; then
  echo "Run mode: OpenMP (threads=${OMP_THREADS})"
else
  echo "Run mode: Serial"
fi
echo "Expected log file: ${CASE_LOG_FILE}"
echo "========================================="
echo ""

mkdir -p "$CASE_DIR"
cp "$PARAM_FILE" "$CASE_DIR/case.params"
cp "$SRC_FILE_ORIG" "$CASE_DIR/$SRC_FILE_LOCAL"

cd "$CASE_DIR"

echo "Compiling ${SRC_FILE_LOCAL} ..."
QCC_FLAGS=(-I../../src-local -O2 -Wall -disable-dimensions)
if [[ $USE_OPENMP -eq 1 ]]; then
  QCC_FLAGS+=(-fopenmp)
fi
if ! qcc "${QCC_FLAGS[@]}" "$SRC_FILE_LOCAL" -o "$EXECUTABLE_NAME" -lm; then
  if [[ $USE_OPENMP -eq 1 ]]; then
    echo "ERROR: OpenMP build failed. Re-run with --threads 1 for serial mode." >&2
  fi
  exit 1
fi
echo "Compilation successful: $EXECUTABLE_NAME"
echo ""

if [[ -f "restart" ]]; then
  echo "Restart file found - simulation will resume from checkpoint."
fi

if [[ $USE_OPENMP -eq 1 ]]; then
  echo "Running: OMP_NUM_THREADS=${OMP_THREADS} ./${EXECUTABLE_NAME} case.params"
  if OMP_NUM_THREADS="$OMP_THREADS" ./"$EXECUTABLE_NAME" case.params; then
    EXIT_CODE=0
  else
    EXIT_CODE=$?
  fi
else
  echo "Running (serial): ./${EXECUTABLE_NAME} case.params"
  if ./"$EXECUTABLE_NAME" case.params; then
    EXIT_CODE=0
  else
    EXIT_CODE=$?
  fi
fi

echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "Simulation completed successfully."
  echo "Output location: simulationCases/${CASE_NO}/"
  if [[ -f "${CASE_LOG_FILE}" ]]; then
    echo "Log file: simulationCases/${CASE_NO}/${CASE_LOG_FILE}"
  fi
else
  echo "Simulation failed with exit code: $EXIT_CODE"
fi

exit "$EXIT_CODE"
