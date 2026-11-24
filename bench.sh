#!/bin/bash
# Run the Enzyme Fortran benchmarks directly from their lit Output directories.
#
# Usage: bench.sh [options] [benchmark ...]
#
# Options:
#   -r, --repetitions N   repetition count passed to each benchmark (default: 9)
#   -c, --case N          which compiled case to run: 1 = flang IR piped through
#                         opt+Enzyme, 2 = flang with the FlangEnzyme plugin
#                         (default: both)
#   -b, --build-dir DIR   CMake build directory (default: ./build)
#   -l, --list            list the discovered benchmarks and exit
#   -h, --help            show this help
#
# Positional arguments filter the benchmarks by name: each is matched as a
# substring against the benchmark file names (e.g. "bratu" selects all three
# bratu cases, "square" the square benchmark). With no arguments, every
# discovered benchmark is run.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"

REPETITIONS=9
CASES="1 2"
BUILD_DIR="${REPO_ROOT}/build"
LIST_ONLY=0
FILTERS=()

usage() {
  sed -n '2,17p' "${BASH_SOURCE[0]}"
}

while [ $# -gt 0 ]; do
  case "$1" in
    -r|--repetitions)
      REPETITIONS="${2:?--repetitions requires a value}"
      shift 2
      ;;
    -r=*|--repetitions=*)
      REPETITIONS="${1#*=}"
      shift
      ;;
    -c|--case)
      CASES="${2:?--case requires a value}"
      shift 2
      ;;
    -c=*|--case=*)
      CASES="${1#*=}"
      shift
      ;;
    -b|--build-dir)
      BUILD_DIR="${2:?--build-dir requires a value}"
      shift 2
      ;;
    -b=*|--build-dir=*)
      BUILD_DIR="${1#*=}"
      shift
      ;;
    -l|--list)
      LIST_ONLY=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [ $# -gt 0 ]; do
        FILTERS+=("$1")
        shift
      done
      ;;
    -*)
      echo "unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      FILTERS+=("$1")
      shift
      ;;
  esac
done

# The lit exec root mirrors the source layout: each benchmark file
# benchmarks/<group>/<name>.f90 gets its compiled binaries under
# ${BUILD_DIR}/benchmarks/benchmarks/<group>/Output/<name>.f90.tmp{1,2}.
EXEC_ROOT="${BUILD_DIR}/benchmarks/benchmarks"
if [ ! -d "${EXEC_ROOT}" ]; then
  echo "No benchmark outputs under ${EXEC_ROOT}." >&2
  echo "Run the lit suite first, e.g. ./test.sh (or cmake --build ${BUILD_DIR} --target bench-enzyme-fortran)." >&2
  exit 1
fi

# Discover benchmarks as "<group>/<name>" pairs from the lit Output dirs.
BENCHMARKS=()
while IFS= read -r script; do
  name="$(basename "${script}" .f90.script)"
  group="$(basename "$(dirname "$(dirname "${script}")")")"
  BENCHMARKS+=("${group}/${name}")
done < <(find "${EXEC_ROOT}" -name '*.f90.script' | sort)

if [ ${#BENCHMARKS[@]} -eq 0 ]; then
  echo "No benchmarks discovered under ${EXEC_ROOT}." >&2
  exit 1
fi

matches_filter() {
  local bench="$1"
  if [ ${#FILTERS[@]} -eq 0 ]; then
    return 0
  fi
  local filter
  for filter in "${FILTERS[@]}"; do
    case "${bench}" in
      *"${filter}"*) return 0 ;;
    esac
  done
  return 1
}

if [ "${LIST_ONLY}" -eq 1 ]; then
  for bench in "${BENCHMARKS[@]}"; do
    if matches_filter "${bench}"; then
      echo "${bench}"
    fi
  done
  exit 0
fi

RUN=0
for bench in "${BENCHMARKS[@]}"; do
  if ! matches_filter "${bench}"; then
    continue
  fi
  group="${bench%%/*}"
  name="${bench#*/}"
  output_dir="${EXEC_ROOT}/${group}/Output"
  for CASE in ${CASES}; do
    binary="${output_dir}/${name}.f90.tmp${CASE}"
    if [ ! -x "${binary}" ]; then
      echo "skipping ${bench}, case ${CASE}: ${binary} not found (was that RUN line executed?)"
      continue
    fi
    echo "${bench}, case ${CASE} (repetitions: ${REPETITIONS})"
    echo
    "${binary}" "${REPETITIONS}"
    echo
    RUN=$((RUN + 1))
  done
done

if [ "${RUN}" -eq 0 ]; then
  echo "No benchmarks matched. Available benchmarks:" >&2
  for bench in "${BENCHMARKS[@]}"; do
    echo "  ${bench}" >&2
  done
  exit 1
fi
