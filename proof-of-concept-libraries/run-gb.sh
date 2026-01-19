#!/usr/bin/env bash

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/template.xlsx
CONFIG_PATH=/Volumes/T7/experiments/experiment-GB/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
/Volumes/T7/experiments/experiment-GB/demultiplex/cutadapt_input_files/demultiplex.sh
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=False
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=False --sort_by_counts=False

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit library enumerate --config_path=$CONFIG_PATH --graph_only=True
delt-hit library enumerate --config_path=$CONFIG_PATH --debug=invalid

delt-hit analyse enrichment \
  --config_path=/Volumes/T7/experiments/experiment-GB/analysis.yaml \
  --name=analysis-1 \
  --method=counts
Rscript --vanilla /Volumes/T7/experiments/experiment-GB/analysis/analysis-1/counts/enrichment_counts.R

delt-hit analyse enrichment \
  --config_path=/Volumes/T7/experiments/experiment-GB/analysis.yaml \
  --name=analysis-1 \
  --method=edgeR
Rscript --vanilla /Volumes/T7/experiments/experiment-GB/analysis/analysis-1/edgeR/enrichment_edgeR.R