#!/usr/bin/env bash

# Example Slurm submissions:
# sbatch --job-name=favalli-lane-1-fasta --mem=32G --time=04:00:00 --cpus-per-task=12 --output="$HOME/logs/%j.out" --wrap="bash supporting_material/experiments/favalli/run.sh"

set -euo pipefail

cd supporting_material/experiments/favalli || exit

PREFIX="lane-1"

time pixi run delt-hit init --excel_path "$PREFIX.xlsx"

CONFIG_PATH=$PREFIX/config.yaml

time pixi run delt-hit demultiplex prepare --config_path="$CONFIG_PATH"
time bash "$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh"

time pixi run delt-hit demultiplex report --config_path="$CONFIG_PATH"
time pixi run delt-hit demultiplex qc --config_path="$CONFIG_PATH"

time pixi run delt-hit demultiplex process --config_path="$CONFIG_PATH"
time pixi run delt-hit demultiplex process --config_path="$CONFIG_PATH" --as_files=True
