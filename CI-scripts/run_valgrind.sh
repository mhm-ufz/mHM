#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Run the curated serial memory checks and one routing memory profile.

Usage:
  run_valgrind.sh [--log-dir DIR] DRIVER REINITIALIZE_TEST

DRIVER and REINITIALIZE_TEST must be built with the same serial Debug toolchain.
Memcheck covers routing, land hydrology/output, lake restart write/read, and
repeated initialization. Massif profiles routing only. Restart pairs share a run
directory; other checks use separate staged roots. Logs default to build/valgrind-logs.
EOF
}

log_dir="build/valgrind-logs"
executables=()
while (($# > 0)); do
  case "$1" in
    --log-dir)
      if (($# < 2)); then
        usage >&2
        exit 2
      fi
      log_dir="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      executables+=("$@")
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      executables+=("$1")
      shift
      ;;
  esac
done
if ((${#executables[@]} != 2)); then
  usage >&2
  exit 2
fi
for i in "${!executables[@]}"; do
  if [[ ! -x "${executables[i]}" ]]; then
    echo "Executable not found or not executable: ${executables[i]}" >&2
    exit 2
  fi
  executables[i]=$(cd "$(dirname "${executables[i]}")" && pwd)/$(basename "${executables[i]}")
done
driver="${executables[0]}"
reinitialize_test="${executables[1]}"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/.." && pwd)
cd "${repo_root}"
# shellcheck source=CI-scripts/v6_namelist_helpers.sh
source "${script_dir}/v6_namelist_helpers.sh"
mkdir -p "${log_dir}" build/valgrind-runs
log_dir=$(cd "${log_dir}" && pwd)
run_base="${repo_root}/build/valgrind-runs"

# Host uprobes in clearenv break Valgrind's libc cleanup; application cleanup still runs.
memcheck_options=(
  --tool=memcheck --error-exitcode=99 --leak-check=full
  "--errors-for-leak-kinds=definite,indirect,possible" --show-leak-kinds=all
  --track-origins=yes --run-libc-freeres=no
)
valgrind --version | tee "${log_dir}/version.log"

stage_run() {
  local label="$1"
  run_dir=$(mktemp -d "${run_base}/${label}.XXXXXX")
  cmake -DSOURCE_DIR="${repo_root}" -DRUN_DIR="${run_dir}" -P "${script_dir}/stage_run.cmake"
  echo "Run root for ${label}: ${run_dir}"
}

run_memcheck() {
  local label="$1"
  shift
  local status start_seconds=${SECONDS}
  local report="${log_dir}/${label}.memcheck.log"
  local model_log="${log_dir}/${label}.model.log"
  echo "Running Memcheck: ${label}"
  if (cd "${run_dir}" && valgrind "${memcheck_options[@]}" --log-file="${report}" "$@") >"${model_log}" 2>&1; then
    status=0
  else
    status=$?
  fi
  if [[ -f "${report}" ]]; then
    tail -n 12 "${report}"
  fi
  if ((status != 0)); then
    echo "Memcheck failed for ${label} (exit ${status}); see ${report}" >&2
    tail -n 50 "${model_log}" >&2
    exit "${status}"
  fi
  echo "Completed ${label} in $((SECONDS - start_seconds)) s"
}

run_namelist() {
  local name="$1"
  local param_file out_file
  param_file=$(parameter_file_for_namelist "${name}")
  out_file=$(output_file_for_namelist "${name}")
  run_memcheck "${name}" "${driver}" -n "test_nml/${name}.nml" -p "${param_file}" -o "${out_file}"
  (cd "${run_dir}" && post_run_check "${name}")
}

stage_run routing
routing_run_dir="${run_dir}"
run_namelist mrm_minimal
stage_run land
run_namelist mhm_output_minimal
stage_run lake_restart
run_namelist mlm_lake_case1_restart_write
run_namelist mlm_lake_case1_restart_read
stage_run reinitialize
run_memcheck test_parameter_reinitialize "${reinitialize_test}"

echo "Profiling routing with Massif"
if (cd "${routing_run_dir}" && valgrind --tool=massif --stacks=yes --time-unit=i --run-libc-freeres=no \
    --massif-out-file="${log_dir}/massif.out" --log-file="${log_dir}/massif.log" \
    "${driver}" -n test_nml/mrm_minimal.nml -p test_nml/mhm_parameter_v6_routing.nml \
    -o mhm-output-template.nml) >"${log_dir}/massif.model.log" 2>&1; then
  # Packaged ms_print can retain a build-time shebang; use the available Perl.
  perl "$(command -v ms_print)" "${log_dir}/massif.out" >"${log_dir}/massif.txt"
  cat "${log_dir}/massif.txt"
else
  status=$?
  tail -n 50 "${log_dir}/massif.model.log" >&2
  exit "${status}"
fi
