#!/usr/bin/env bash

cd /Users/adrianomartinelli/projects/delt/delt-core/supporting_material
delt-hit init --excel_path template.xlsx

CONFIG_PATH=experiments/template/config.yaml

delt-hit library enumerate --config_path=$CONFIG_PATH --graph_only=True  # produce reaction graph only for inspection
delt-hit library enumerate --config_path=$CONFIG_PATH

delt-hit library properties --config_path=$CONFIG_PATH

delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
experiments/template/demultiplex/cutadapt_input_files/demultiplex.sh
delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True --sort_by_counts=False

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

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