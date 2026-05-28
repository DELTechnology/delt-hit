#!/usr/bin/env bash

set -euo pipefail

cd supporting_material/experiments/example-single-stranded || exit

delt-hit init --excel_path example-single-stranded.xlsx

CONFIG_PATH=campaign/config.yaml
ANALYSIS_CONFIG_PATH=analysis.yaml
ANALYSIS_OUTPUT_ROOT=campaign/analysis

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
campaign/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True
# delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True --sort_by_counts=False

delt-hit visualize enumerate --config_path=$CONFIG_PATH

delt-hit library enumerate \
  --config_path=$CONFIG_PATH \
  --counts_path=campaign/selections/AG24_4_counts.txt \
  --top_n=1000 \
  --library_name=AG24_4_top_hits

delt-hit visualize library --config_path=$CONFIG_PATH --library_name=AG24_4_top_hits
delt-hit library properties --config_path=$CONFIG_PATH --library_name=AG24_4_top_hits

delt-hit dashboard \
  --config_path=$CONFIG_PATH \
  --counts_path=campaign/selections/AG24_4_counts.txt

delt-hit analyse enrichment \
  --analysis_config=$ANALYSIS_CONFIG_PATH \
  --name=protein_vs_no_protein \
  --method=counts \
  --save_dir=$ANALYSIS_OUTPUT_ROOT/protein_vs_no_protein

Rscript --vanilla campaign/analysis/protein_vs_no_protein/counts/enrichment_counts.R

delt-hit analyse enrichment \
  --analysis_config=$ANALYSIS_CONFIG_PATH \
  --name=protein_vs_no_protein \
  --method=edgeR \
  --save_dir=$ANALYSIS_OUTPUT_ROOT/protein_vs_no_protein

Rscript --vanilla campaign/analysis/protein_vs_no_protein/edgeR/enrichment_edgeR.R

for selection in AG24_13 AG24_14 AG24_15; do
  delt-hit analyse enrichment \
    --config_path=$CONFIG_PATH \
    --counts=campaign/selections/$selection/counts.txt \
    --method=z_score \
    --save_dir=$ANALYSIS_OUTPUT_ROOT/z_score/$selection

  Rscript --vanilla $ANALYSIS_OUTPUT_ROOT/z_score/$selection/z_score/enrichment_z_score.R
done

# full library enumeration, usually only needed for ML tasks
delt-hit library enumerate --config_path=$CONFIG_PATH
delt-hit library properties --config_path=$CONFIG_PATH

delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH
