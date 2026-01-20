#!/usr/bin/env bash

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/proof-of-concept-libraries/template.xlsx
CONFIG_PATH=/Volumes/T7/experiments/template/config.yaml

delt-hit library enumerate --config_path=$CONFIG_PATH --graph_only=True || exit
delt-hit library enumerate --config_path=$CONFIG_PATH --debug=invalid || exit

delt-hit demultiplex prepare --config_path=$CONFIG_PATH || exit
/Volumes/T7/experiments/template/demultiplex/cutadapt_input_files/demultiplex.sh || exit
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True --sort_by_counts=False || exit
delt-hit demultiplex process --config_path=$CONFIG_PATH || exit

delt-hit demultiplex report --config_path=$CONFIG_PATH || exit
delt-hit demultiplex qc --config_path=$CONFIG_PATH || exit

delt-hit dashboard \
  --config_path=$CONFIG_PATH \
  --counts_path=/Volumes/T7/experiments/template/selections/AG24_4/counts.txt

delt-hit analyse enrichment \
  --config_path=/Volumes/T7/experiments/template/analysis.yaml \
  --name=protein_vs_no_protein \
  --method=counts
Rscript --vanilla /Volumes/T7/experiments/template/analysis/protein_vs_no_protein/counts/enrichment_counts.R

delt-hit analyse enrichment \
  --config_path=/Volumes/T7/experiments/template/analysis.yaml \
  --name=protein_vs_no_protein \
  --method=edgeR
Rscript --vanilla /Volumes/T7/experiments/template/analysis/protein_vs_no_protein/edgeR/enrichment_edgeR.R