#!/usr/bin/env bash

set -euo pipefail

PREFIX="${1:-lane-1}"

time pixi run delt-hit init --excel_path "$PREFIX.xlsx"

CONFIG_PATH=$PREFIX/config.yaml

time pixi run delt-hit demultiplex prepare --config_path="$CONFIG_PATH"
time bash "$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh"

time pixi run delt-hit demultiplex report --config_path="$CONFIG_PATH"
time pixi run delt-hit demultiplex qc --config_path="$CONFIG_PATH"

time pixi run delt-hit demultiplex process --config_path="$CONFIG_PATH"
time pixi run delt-hit demultiplex process --config_path="$CONFIG_PATH" --as_files=True
