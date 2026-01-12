#!/usr/bin/env bash

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/template.xlsx
delt-hit demultiplex prepare --config_path=/Volumes/T7/experiments/experiment-GB/config.yaml
/Volumes/T7/experiments/experiment-GB/demultiplex/cutadapt_input_files/demultiplex.sh
delt-hit demultiplex process --config_path=/Volumes/T7/experiments/experiment-GB/config.yaml --as_files=True
delt-hit demultiplex process --config_path=/Volumes/T7/experiments/experiment-GB/config.yaml --as_files=False
delt-hit demultiplex process --config_path=/Volumes/T7/experiments/experiment-GB/config.yaml --as_files=False --sort_by_counts=False

delt-hit demultiplex report --config_path=/Volumes/T7/experiments/experiment-GB/config.yaml
delt-hit demultiplex qc --config=/Volumes/T7/experiments/experiment-GB/config.yaml