#!/usr/bin/env bash

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/libraries/260115_NF2_library_rechecked_with_precursors.xlsx
CONFIG_PATH=/Volumes/T7/experiments/experiment-NF2-with-precursors/config.yaml

delt-hit library enumerate --config_path=$CONFIG_PATH --graph_only=True
delt-hit library enumerate --config_path=$CONFIG_PATH
delt-hit library properties --config_path=$CONFIG_PATH
delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
/Volumes/T7/experiments/experiment-NF2/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True

delt-hit dashboard \
  --config_path=$CONFIG_PATH \
  --counts_path=/Volumes/T7/experiments/experiment-NF2-with-precursors/selections/AG24_1_counts.txt
