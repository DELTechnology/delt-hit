#!/usr/bin/env bash

cd /Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli || exit

PREFIX=lane-1
#PREFIX=lane-2

delt-hit init --excel_path $PREFIX.xlsx

CONFIG_PATH=$PREFIX/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True

delt-hit visualize enumerate --config_path=$CONFIG_PATH

delt-hit library enumerate \
  --config_path=$CONFIG_PATH \
  --counts_path=$PREFIX/selections/AG24_4/counts.txt \
  --top_n=1000 \
  --library_name=AG24_4_top_hits

delt-hit library properties --config_path=$CONFIG_PATH --library_name=AG24_4_top_hits

delt-hit dashboard \
  --config_path=$CONFIG_PATH \
  --counts_path=$PREFIX/selections/AG24_4/counts.txt

delt-hit analyse enrichment \
  --config_path=analysis.yaml \
  --name=protein_vs_no_protein \
  --method=counts

Rscript --vanilla $PREFIX/analysis/protein_vs_no_protein/counts/enrichment_counts.R

delt-hit analyse enrichment \
  --config_path=analysis.yaml \
  --name=protein_vs_no_protein \
  --method=edgeR

Rscript --vanilla $PREFIX/analysis/protein_vs_no_protein/edgeR/enrichment_edgeR.R

# full library enumeration, usually only needed for ML tasks
delt-hit library enumerate --config_path=$CONFIG_PATH
delt-hit library properties --config_path=$CONFIG_PATH

delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH