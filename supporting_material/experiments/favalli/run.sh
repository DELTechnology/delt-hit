#!/usr/bin/env bash

# Example Slurm submissions:
# sbatch --job-name=favalli-lane-1 --mem=32G --time=04:00:00 --cpus-per-task=12 --output="$HOME/logs/%j.out" --wrap="bash supporting_material/experiments/favalli/run.sh lane-1"
# sbatch --job-name=favalli-lane-2-fasta --mem=32G --time=04:00:00 --cpus-per-task=12 --output="$HOME/logs/%j.out" --wrap="bash supporting_material/experiments/favalli/run.sh lane-2-fasta"

set -euo pipefail

cd supporting_material/experiments/favalli || exit

delt-hit init --excel_path "lane-1.xlsx"

CONFIG_PATH=lane-1/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
lane-1/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True
