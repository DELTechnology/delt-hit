#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DELI_ROOT="${ROOT}/other_tools/DELi"
DELI_VENV="${ROOT}/other_tools/deli/.venv"

export VIRTUAL_ENV="${DELI_VENV}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/deli-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/deli-cache}"

run_deli() {
  uv run --active --no-project deli "$@"
}

for experiment_name in synthetic_2cycle synthetic_3cycle synthetic_4cycle; do
  experiment_dir="${DELI_ROOT}/${experiment_name}"
  output_dir="${experiment_dir}/output_v02"
  prefix="${experiment_name}_v02"

  rm -rf "${output_dir}"
  mkdir -p "${output_dir}"

  run_deli --config-file "${experiment_dir}/deli_config" --disable-logging \
    decode run "${experiment_dir}/decode_synthetic.yaml" \
    --decode-settings-file "${experiment_dir}/decode_settings_v02.yaml" \
    --out-dir "${output_dir}" \
    --prefix "${prefix}" \
    --save-failed \
    --skip-report

  run_deli --config-file "${experiment_dir}/deli_config" --disable-logging \
    decode collect "${output_dir}/${prefix}_decoded.tsv" \
    --out-loc "${output_dir}/${prefix}_collected.ndjson"

  run_deli --config-file "${experiment_dir}/deli_config" --disable-logging \
    decode count "${output_dir}/${prefix}_collected.ndjson" \
    --out-loc "${output_dir}/${prefix}_counts.parquet" \
    --keep-raw-count \
    --output-format parquet

  run_deli --config-file "${experiment_dir}/deli_config" --disable-logging \
    decode summarize \
    "${output_dir}/${prefix}_counts.parquet" \
    "${output_dir}/${prefix}_decode_statistics.json" \
    --out-loc "${output_dir}/${prefix}_summary.json"
done
