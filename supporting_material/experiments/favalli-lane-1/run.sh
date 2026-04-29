#!/usr/bin/env bash

cd /users/amarti51/projects/delt-hit/supporting_material/experiments/favalli-lane-1
CONFIG_PATH=config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
demultiplex/cutadapt_input_files/demultiplex.sh
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