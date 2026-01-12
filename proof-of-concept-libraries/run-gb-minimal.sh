#!/usr/bin/env bash

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/template_minimal.xlsx

delt-hit demultiplex prepare --config_path=/Volumes/T7/experiments/experiment-GB-minimal/config.yaml

/Volumes/T7/experiments/experiment-GB-minimal/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=/Volumes/T7/experiments/experiment-GB-minimal/config.yaml
delt-hit demultiplex qc --config_path=/Volumes/T7/experiments/experiment-GB-minimal/config.yaml

delt-hit demultiplex process --config_path=/Volumes/T7/experiments/experiment-GB-minimal/config.yaml --as_files=True

