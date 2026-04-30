#!/usr/bin/env bash

delt-hit init --excel_path NFv3_JP_1_c0.1percent.xlsx

CONFIG_PATH=JP-1/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
JP-1/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True

delt-hit visualize enumerate --config_path=$CONFIG_PATH

delt-hit library enumerate \
  --config_path=$CONFIG_PATH \
  --counts_path=experiments/template/selections/AG24_4/counts.txt \
  --top_n=1000 \
  --library_name=AG24_4_top_hits

delt-hit library properties --config_path=$CONFIG_PATH --library_name=AG24_4_top_hits

delt-hit dashboard \
  --config_path=$CONFIG_PATH \
  --counts_path=experiments/template/selections/AG24_4/counts.txt

for ANALYSIS in his_pure_protein_vs_no_protein;
#for ANALYSIS in his_pure_protein_vs_no_protein dyna_protein_vs_no_protein;
do
  delt-hit analyse enrichment \
    --config_path=analysis.yaml \
    --name=$ANALYSIS \
    --method=counts

  Rscript --vanilla JP-1/analysis/$ANALYSIS/counts/enrichment_counts.R

  delt-hit analyse enrichment \
    --config_path=analysis.yaml \
    --name=$ANALYSIS \
    --method=edgeR

  Rscript --vanilla JP-1/analysis/$ANALYSIS/edgeR/enrichment_edgeR.R
done

# full library enumeration, usually only needed for ML tasks
delt-hit library enumerate --config_path=$CONFIG_PATH
delt-hit library properties --config_path=$CONFIG_PATH

delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH