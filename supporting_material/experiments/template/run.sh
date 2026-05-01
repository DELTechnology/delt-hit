#!/usr/bin/env bash

delt-hit init --excel_path template.xlsx

CONFIG_PATH=experiments/template/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
experiments/template/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True
# delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True --sort_by_counts=False

delt-hit visualize enumerate --config_path=$CONFIG_PATH

delt-hit library enumerate \
  --config_path=$CONFIG_PATH \
  --counts_path=experiments/template/selections/AG24_4/counts.txt \
  --top_n=1000 \
  --library_name=AG24_4_top_hits

delt-hit visualize library --config_path=$CONFIG_PATH --library_name=AG24_4_top_hits
delt-hit library properties --config_path=$CONFIG_PATH --library_name=AG24_4_top_hits

delt-hit dashboard \
  --config_path=$CONFIG_PATH \
  --counts_path=experiments/template/selections/AG24_4/counts.txt

delt-hit analyse enrichment \
  --config_path=analysis.yaml \
  --name=protein_vs_no_protein \
  --method=counts

Rscript --vanilla experiments/template/analysis/protein_vs_no_protein/counts/enrichment_counts.R

delt-hit analyse enrichment \
  --config_path=analysis.yaml \
  --name=protein_vs_no_protein \
  --method=edgeR

Rscript --vanilla experiments/template/analysis/protein_vs_no_protein/edgeR/enrichment_edgeR.R

# full library enumeration, usually only needed for ML tasks
delt-hit library enumerate --config_path=$CONFIG_PATH
delt-hit library properties --config_path=$CONFIG_PATH

delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH
