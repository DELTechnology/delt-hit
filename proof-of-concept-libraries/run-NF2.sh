#!/usr/bin/env bash

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/libraries/251021_NF2_library_rechecked.xlsx

delt-hit library enumerate --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml
delt-hit library properties --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml
delt-hit library represent --method=morgan --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml
delt-hit library represent --method=bert --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml

delt-hit demultiplex prepare --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml
/Volumes/T7/experiments/experiment-NF2/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml
delt-hit demultiplex qc --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml

delt-hit demultiplex process --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml --as_files=True

delt-hit dashboard \
  --config_path=/Volumes/T7/experiments/experiment-NF2/config.yaml \
  --counts_path=/Volumes/T7/experiments/experiment-NF2/selections/AG24_1_counts.txt
