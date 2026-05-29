#!/usr/bin/env bash

set -euo pipefail

cd supporting_material/experiments/pure-del || exit

PREFIX="${1:-lane-1}"

pixi run delt-hit init --excel_path "$PREFIX.xlsx"

CONFIG_PATH=$PREFIX/config.yaml

pixi run delt-hit demultiplex prepare --config_path=$CONFIG_PATH
$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh

pixi run delt-hit demultiplex report --config_path=$CONFIG_PATH
pixi run delt-hit demultiplex qc --config_path=$CONFIG_PATH

pixi run delt-hit demultiplex process --config_path=$CONFIG_PATH
pixi run delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True
