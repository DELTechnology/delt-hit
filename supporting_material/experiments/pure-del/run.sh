#!/usr/bin/env bash

set -euo pipefail

# Run from this directory after: conda activate delt-hit

PREFIX="${1:-lane-1}"

time delt-hit init --excel_path "$PREFIX.xlsx"

CONFIG_PATH=$PREFIX/config.yaml

time delt-hit demultiplex prepare --config_path="$CONFIG_PATH"
time bash "$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh"

time delt-hit demultiplex report --config_path="$CONFIG_PATH"
time delt-hit demultiplex qc --config_path="$CONFIG_PATH"

time delt-hit demultiplex process --config_path="$CONFIG_PATH"
time delt-hit demultiplex process --config_path="$CONFIG_PATH" --as_files=True
