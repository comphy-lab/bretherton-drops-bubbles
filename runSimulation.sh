#!/bin/bash
# runSimulation.sh
#
# Run a single bretherton-drops-bubbles simulation from the repository root.
# The script creates simulationCases/<CaseNo>/, copies the parameter file and
# source file, compiles the selected case, and runs it.
#
# Usage:
#   bash runSimulation.sh [params_file] [OPTIONS]

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
  --openmp       Force an OpenMP build when --threads is 1
  --ranks N      Build and run with N MPI ranks
  --rankfile FILE
                 Open MPI rankfile; mutually exclusive with --pe-list
  --pe-list LIST Explicit comma-separated core indices, one per MPI rank
  --mpi-timeout SEC
                 Ask mpirun to terminate the MPI job after SEC seconds
  --build-only   Compile the case-local executable without running it
  --no-build     Run a matching case-local executable without compiling
  -h, --help     Show this help message

Environment:
  OUTPUT_ROOT    Parent directory for numbered case directories
  MPIEXEC        MPI launcher executable or wrapper (default: mpirun)
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

sha256_file() {
  "${SHA256_COMMAND[@]}" "$1" | awk '{print $1}'
}

sha256_stream() {
  "${SHA256_COMMAND[@]}" | awk '{print $1}'
}

select_sha256_backend() {
  if command -v sha256sum >/dev/null 2>&1; then
    SHA256_COMMAND=(sha256sum --)
    SHA256_BACKEND_NAME="sha256sum"
  elif command -v shasum >/dev/null 2>&1; then
    SHA256_COMMAND=(shasum -a 256 --)
    SHA256_BACKEND_NAME="shasum -a 256"
  else
    echo "ERROR: SHA-256 utility not found (need sha256sum or shasum)." >&2
    exit 1
  fi
}

command_version() {
  local executable="$1"
  local version_output
  if version_output="$("$executable" --version 2>&1)"; then
    printf '%s\n' "${version_output%%$'\n'*}"
  else
    printf '%s\n' "unavailable"
  fi
}

build_input_sha256() {
  local source_file="$1"
  local header
  local LC_ALL=C
  local headers=("${SCRIPT_DIR}"/src-local/*.h)
  {
    printf '%s  %s\n' "$(sha256_file "$source_file")" "simulationCases/$(basename "$source_file")"
    for header in "${headers[@]}"; do
      [[ -f "$header" ]] || continue
      printf '%s  %s\n' "$(sha256_file "$header")" "${header#"${SCRIPT_DIR}/"}"
    done
    if [[ -f "${SCRIPT_DIR}/basilisk/.comphy-lock" ]]; then
      printf '%s  %s\n' "$(sha256_file "${SCRIPT_DIR}/basilisk/.comphy-lock")" \
        "basilisk/.comphy-lock"
    fi
  } | sha256_stream
}

# Defaults
EXEC_CODE="bretherton.c"
PARAM_FILE="default.params"
PARAM_FILE_SET=0
OMP_THREADS=1
FORCE_OPENMP=0
MPI_RANKS=""
MPI_RANKFILE=""
MPI_PE_LIST=""
MPI_TIMEOUT=""
BUILD_ONLY=0
NO_BUILD=0
MPIEXEC="${MPIEXEC:-mpirun}"

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
    --openmp)
      FORCE_OPENMP=1
      shift
      ;;
    --ranks)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: $1 requires a positive integer value." >&2
        usage
        exit 1
      fi
      MPI_RANKS="$2"
      shift 2
      ;;
    --ranks=*)
      MPI_RANKS="${1#*=}"
      shift
      ;;
    --rankfile)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: $1 requires a file path." >&2
        usage
        exit 1
      fi
      MPI_RANKFILE="$2"
      shift 2
      ;;
    --rankfile=*)
      MPI_RANKFILE="${1#*=}"
      shift
      ;;
    --pe-list)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: $1 requires explicit comma-separated core indices." >&2
        usage
        exit 1
      fi
      MPI_PE_LIST="$2"
      shift 2
      ;;
    --pe-list=*)
      MPI_PE_LIST="${1#*=}"
      shift
      ;;
    --mpi-timeout)
      if [[ -z "${2:-}" ]]; then
        echo "ERROR: $1 requires a positive integer value." >&2
        usage
        exit 1
      fi
      MPI_TIMEOUT="$2"
      shift 2
      ;;
    --mpi-timeout=*)
      MPI_TIMEOUT="${1#*=}"
      shift
      ;;
    --build-only)
      BUILD_ONLY=1
      shift
      ;;
    --no-build)
      NO_BUILD=1
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

if [[ -n "$MPI_RANKS" && ! "$MPI_RANKS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --ranks must be a positive integer, got: $MPI_RANKS" >&2
  exit 1
fi

if [[ -n "$MPI_RANKS" && "$OMP_THREADS" -gt 1 ]]; then
  echo "ERROR: --ranks cannot be combined with --threads greater than 1." >&2
  exit 1
fi

if [[ -n "$MPI_RANKS" && $FORCE_OPENMP -eq 1 ]]; then
  echo "ERROR: --ranks cannot be combined with --openmp." >&2
  exit 1
fi

if [[ -n "$MPI_RANKFILE" && -n "$MPI_PE_LIST" ]]; then
  echo "ERROR: --rankfile and --pe-list are mutually exclusive." >&2
  exit 1
fi

if [[ -n "$MPI_RANKS" && -z "$MPI_RANKFILE" && -z "$MPI_PE_LIST" ]]; then
  echo "ERROR: either --rankfile or --pe-list is required with --ranks." >&2
  exit 1
fi

if [[ -z "$MPI_RANKS" && ( -n "$MPI_RANKFILE" || -n "$MPI_PE_LIST" ) ]]; then
  echo "ERROR: --rankfile and --pe-list require --ranks." >&2
  exit 1
fi

if [[ -n "$MPI_TIMEOUT" && ! "$MPI_TIMEOUT" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --mpi-timeout must be a positive integer, got: $MPI_TIMEOUT" >&2
  exit 1
fi

if [[ -z "$MPI_RANKS" && -n "$MPI_TIMEOUT" ]]; then
  echo "ERROR: --mpi-timeout requires --ranks." >&2
  exit 1
fi

if [[ $BUILD_ONLY -eq 1 && $NO_BUILD -eq 1 ]]; then
  echo "ERROR: --build-only and --no-build are mutually exclusive." >&2
  exit 1
fi

if [[ -n "$MPI_PE_LIST" ]]; then
  if [[ ! "$MPI_PE_LIST" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
    echo "ERROR: --pe-list must contain explicit comma-separated non-negative integers." >&2
    exit 1
  fi
  IFS=',' read -r -a MPI_PE_CORES <<< "$MPI_PE_LIST"
  if [[ ${#MPI_PE_CORES[@]} -ne $MPI_RANKS ]]; then
    echo "ERROR: --pe-list must contain exactly ${MPI_RANKS} core indices." >&2
    exit 1
  fi
  MPI_PE_LIST_NORMALIZED=""
  MPI_PE_CORES_SEEN=()
  for core_index in "${MPI_PE_CORES[@]}"; do
    if [[ ! "$core_index" =~ ^(0|[1-9][0-9]*)$ ]]; then
      echo "ERROR: --pe-list indices must use canonical decimal notation: $core_index" >&2
      exit 1
    fi
    for seen_core in "${MPI_PE_CORES_SEEN[@]}"; do
      if [[ "$core_index" == "$seen_core" ]]; then
        echo "ERROR: --pe-list contains duplicate core index: $core_index" >&2
        exit 1
      fi
    done
    MPI_PE_CORES_SEEN+=("$core_index")
    MPI_PE_LIST_NORMALIZED+="${MPI_PE_LIST_NORMALIZED:+,}${core_index}"
  done
  MPI_PE_LIST="$MPI_PE_LIST_NORMALIZED"
fi

USE_OPENMP=0
if [[ "$OMP_THREADS" -gt 1 || $FORCE_OPENMP -eq 1 ]]; then
  USE_OPENMP=1
fi

USE_MPI=0
if [[ -n "$MPI_RANKS" ]]; then
  USE_MPI=1
fi

if [[ $USE_MPI -eq 1 ]]; then
  if [[ -n "$MPI_RANKFILE" ]]; then
    if [[ ! -f "$MPI_RANKFILE" || ! -r "$MPI_RANKFILE" ]]; then
      echo "ERROR: MPI rankfile is not a readable regular file: $MPI_RANKFILE" >&2
      exit 1
    fi
    MPI_RANKFILE="$(cd "$(dirname "$MPI_RANKFILE")" && pwd)/$(basename "$MPI_RANKFILE")"
  fi
fi

if [[ -z "$EXEC_CODE" || "$EXEC_CODE" == */* || "$EXEC_CODE" == -* ]]; then
  echo "ERROR: --exec must be a basename within simulationCases/: $EXEC_CODE" >&2
  exit 1
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

if [[ $NO_BUILD -eq 0 ]]; then
  if ! command -v qcc >/dev/null 2>&1; then
    echo "ERROR: qcc not found in PATH." >&2
    echo "Hint: source your Basilisk environment or provide .project_config." >&2
    exit 1
  fi
  QCC_PATH="$(command -v qcc)"
fi

if [[ $USE_MPI -eq 1 && $NO_BUILD -eq 0 ]]; then
  if ! command -v mpicc >/dev/null 2>&1; then
    echo "ERROR: mpicc not found in PATH." >&2
    exit 1
  fi
  MPICC_PATH="$(command -v mpicc)"
fi

if [[ $USE_MPI -eq 1 && $BUILD_ONLY -eq 0 ]]; then
  if ! MPIEXEC_PATH="$(command -v "$MPIEXEC")"; then
    echo "ERROR: MPIEXEC is not an executable command: $MPIEXEC" >&2
    exit 1
  fi
fi

select_sha256_backend

if [[ ! -f "$PARAM_FILE" ]]; then
  echo "ERROR: Parameter file not found: $PARAM_FILE" >&2
  exit 1
fi

# Run data is a separate surface from source. OUTPUT_ROOT defaults to the
# in-tree location so existing invocations are unchanged, but a dispatched
# campaign points it at an organised run directory on a registered volume.
OUTPUT_ROOT="${OUTPUT_ROOT:-${SCRIPT_DIR}/simulationCases}"
if [[ ! -d "$OUTPUT_ROOT" ]]; then
  mkdir -p "$OUTPUT_ROOT" || {
    echo "ERROR: cannot create output root: $OUTPUT_ROOT" >&2
    exit 1
  }
fi
OUTPUT_ROOT="$(cd "$OUTPUT_ROOT" && pwd)"

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

CASE_DIR="${OUTPUT_ROOT}/${CASE_NO}"
SRC_FILE_LOCAL="${EXEC_CODE}"
EXECUTABLE_NAME="${EXEC_CODE%.c}"
CASE_LOG_FILE="c${CASE_NO}-log"
BUILD_METADATA="${CASE_DIR}/build.meta"

if [[ $USE_MPI -eq 1 ]]; then
  BUILD_MODE="mpi"
elif [[ $USE_OPENMP -eq 1 ]]; then
  BUILD_MODE="openmp"
else
  BUILD_MODE="serial"
fi

echo "========================================="
echo "bretherton-drops-bubbles - Case Runner"
echo "========================================="
echo "Source file: ${EXEC_CODE}"
echo "Parameter file: ${PARAM_FILE}"
echo "CaseNo: ${CASE_NO}"
echo "Case directory: ${CASE_DIR}"
if [[ $USE_MPI -eq 1 ]]; then
  echo "Run mode: MPI (ranks=${MPI_RANKS})"
  if [[ $BUILD_ONLY -eq 0 ]]; then
    echo "MPI launcher: ${MPIEXEC_PATH}"
  fi
  if [[ -n "$MPI_RANKFILE" ]]; then
    echo "Rankfile: ${MPI_RANKFILE}"
  else
    echo "PE list: ${MPI_PE_LIST}"
  fi
elif [[ $USE_OPENMP -eq 1 ]]; then
  echo "Run mode: OpenMP (threads=${OMP_THREADS})"
else
  echo "Run mode: Serial"
fi
echo "Expected log file: ${CASE_DIR}/${CASE_LOG_FILE}"
echo "========================================="
echo ""

mkdir -p "$CASE_DIR"

SOURCE_SHA256="$(sha256_file "$SRC_FILE_ORIG")"
BUILD_INPUT_SHA256="$(build_input_sha256 "$SRC_FILE_ORIG")"
if [[ $NO_BUILD -eq 1 ]]; then
  if [[ ! -f "$BUILD_METADATA" ]]; then
    echo "ERROR: --no-build requires case-local metadata: $BUILD_METADATA" >&2
    exit 1
  fi
  RECORDED_MODE="$(get_param_value "mode" "$BUILD_METADATA")"
  RECORDED_SOURCE_SHA256="$(get_param_value "source_sha256" "$BUILD_METADATA")"
  RECORDED_BUILD_INPUT_SHA256="$(get_param_value "build_input_sha256" "$BUILD_METADATA")"
  RECORDED_EXECUTABLE_SHA256="$(get_param_value "executable_sha256" "$BUILD_METADATA")"
  if [[ "$RECORDED_MODE" != "$BUILD_MODE" ]]; then
    echo "ERROR: --no-build mode mismatch: requested $BUILD_MODE, built $RECORDED_MODE" >&2
    exit 1
  fi
  if [[ "$RECORDED_SOURCE_SHA256" != "$SOURCE_SHA256" ]]; then
    echo "ERROR: --no-build source mismatch for $SRC_FILE_ORIG" >&2
    exit 1
  fi
  if [[ "$RECORDED_BUILD_INPUT_SHA256" != "$BUILD_INPUT_SHA256" ]]; then
    echo "ERROR: --no-build input mismatch in source, local headers, or Basilisk lock." >&2
    exit 1
  fi
  if [[ ! -f "$CASE_DIR/$SRC_FILE_LOCAL" ||
        "$(sha256_file "$CASE_DIR/$SRC_FILE_LOCAL")" != "$RECORDED_SOURCE_SHA256" ]]; then
    echo "ERROR: case-local source does not match its build metadata." >&2
    exit 1
  fi
  if [[ ! -x "$CASE_DIR/$EXECUTABLE_NAME" ||
        "$(sha256_file "$CASE_DIR/$EXECUTABLE_NAME")" != "$RECORDED_EXECUTABLE_SHA256" ]]; then
    echo "ERROR: case-local executable does not match its build metadata." >&2
    exit 1
  fi
fi

if [[ "$PARAM_FILE" != "$CASE_DIR/case.params" ]]; then
  cp "$PARAM_FILE" "$CASE_DIR/case.params"
fi
if [[ $NO_BUILD -eq 0 ]]; then
  cp "$SRC_FILE_ORIG" "$CASE_DIR/$SRC_FILE_LOCAL"
fi

cd "$CASE_DIR"

if [[ $NO_BUILD -eq 0 ]]; then
  echo "Compiling ${SRC_FILE_LOCAL} ..."
  QCC_FLAGS=(-I"${SCRIPT_DIR}/src-local" -O2 -Wall -disable-dimensions)
  QCC_VERSION="$(command_version "$QCC_PATH")"
  if [[ $USE_MPI -eq 1 ]]; then
    GENERATED_SOURCE="_${SRC_FILE_LOCAL}"
    BUILD_RECIPE="qcc -I<src-local> -O2 -Wall -disable-dimensions -source -D_MPI=1 <source>; mpicc -O2 -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 <generated-source> -o <executable> -lm"
    MPICC_VERSION="$(command_version "$MPICC_PATH")"
    if ! qcc "${QCC_FLAGS[@]}" -source -D_MPI=1 "$SRC_FILE_LOCAL"; then
      echo "ERROR: qcc MPI source generation failed." >&2
      exit 1
    fi
    if ! mpicc -O2 -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 \
        "$GENERATED_SOURCE" -o "$EXECUTABLE_NAME" -lm; then
      echo "ERROR: MPI build failed." >&2
      exit 1
    fi
  else
    if [[ $USE_OPENMP -eq 1 ]]; then
      QCC_FLAGS+=(-fopenmp)
      BUILD_RECIPE="qcc -I<src-local> -O2 -Wall -disable-dimensions -fopenmp <source> -o <executable> -lm"
    else
      BUILD_RECIPE="qcc -I<src-local> -O2 -Wall -disable-dimensions <source> -o <executable> -lm"
    fi
    if ! qcc "${QCC_FLAGS[@]}" "$SRC_FILE_LOCAL" -o "$EXECUTABLE_NAME" -lm; then
      if [[ $USE_OPENMP -eq 1 ]]; then
        echo "ERROR: OpenMP build failed. Re-run with --threads 1 for serial mode." >&2
      fi
      exit 1
    fi
  fi
  EXECUTABLE_SHA256="$(sha256_file "$EXECUTABLE_NAME")"
  METADATA_TMP="${BUILD_METADATA}.tmp.$$"
  {
    printf 'mode=%s\n' "$BUILD_MODE"
    printf 'source_sha256=%s\n' "$SOURCE_SHA256"
    printf 'build_input_sha256=%s\n' "$BUILD_INPUT_SHA256"
    printf 'executable_sha256=%s\n' "$EXECUTABLE_SHA256"
    printf 'sha256_backend=%s\n' "$SHA256_BACKEND_NAME"
    printf 'build_recipe=%s\n' "$BUILD_RECIPE"
    printf 'qcc=%s\n' "$QCC_PATH"
    printf 'qcc_version=%s\n' "$QCC_VERSION"
    printf 'qcc_sha256=%s\n' "$(sha256_file "$QCC_PATH")"
    if [[ $USE_MPI -eq 1 ]]; then
      printf 'mpicc=%s\n' "$MPICC_PATH"
      printf 'mpicc_version=%s\n' "$MPICC_VERSION"
      printf 'mpicc_sha256=%s\n' "$(sha256_file "$MPICC_PATH")"
      printf 'generated_source_sha256=%s\n' "$(sha256_file "$GENERATED_SOURCE")"
    fi
  } > "$METADATA_TMP"
  mv "$METADATA_TMP" "$BUILD_METADATA"
  echo "Compilation successful: $EXECUTABLE_NAME"
  echo "Build metadata: $BUILD_METADATA"
  echo ""
else
  echo "Using existing checksum-matched ${BUILD_MODE} executable: $EXECUTABLE_NAME"
  echo ""
fi

if [[ $BUILD_ONLY -eq 1 ]]; then
  echo "Build-only mode complete."
  exit 0
fi

if [[ -f "restart" ]]; then
  echo "Restart file found - simulation will resume from checkpoint."
fi

if [[ $USE_MPI -eq 1 ]]; then
  RUN_COMMAND=("$MPIEXEC_PATH" --np "$MPI_RANKS")
  if [[ -n "$MPI_RANKFILE" ]]; then
    RUN_COMMAND+=(--map-by "rankfile:file=${MPI_RANKFILE}")
  else
    RUN_COMMAND+=(--map-by "pe-list=${MPI_PE_LIST}:ordered")
  fi
  RUN_COMMAND+=(--bind-to core --report-bindings)
  if [[ -n "$MPI_TIMEOUT" ]]; then
    RUN_COMMAND+=(--timeout "$MPI_TIMEOUT")
  fi
  RUN_COMMAND+=("./${EXECUTABLE_NAME}" case.params)
  echo "Running MPI (${MPI_RANKS} ranks; OpenMP disabled)"
  if OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE "${RUN_COMMAND[@]}"; then
    EXIT_CODE=0
  else
    EXIT_CODE=$?
  fi
elif [[ $USE_OPENMP -eq 1 ]]; then
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
  echo "Output location: ${CASE_DIR}/"
  if [[ -f "${CASE_LOG_FILE}" ]]; then
    echo "Log file: ${CASE_DIR}/${CASE_LOG_FILE}"
  fi
else
  echo "Simulation failed with exit code: $EXIT_CODE"
fi

exit "$EXIT_CODE"
