#!/usr/bin/env bash

set -euo pipefail

# Run from this directory after: conda activate delt-hit

time delt-hit init --excel_path example-single-display.xlsx

CONFIG_PATH=campaign/config.yaml
ANALYSIS_CONFIG_PATH=analysis.yaml
ANALYSIS_OUTPUT_ROOT=campaign/analysis

time delt-hit demultiplex prepare --config_path="$CONFIG_PATH"
time bash campaign/demultiplex/cutadapt_input_files/demultiplex.sh

time delt-hit demultiplex report --config_path="$CONFIG_PATH"
time delt-hit demultiplex qc --config_path="$CONFIG_PATH"

time delt-hit demultiplex process --config_path="$CONFIG_PATH"
time delt-hit demultiplex process --config_path="$CONFIG_PATH" --as_files=True

time delt-hit visualize enumerate --config_path="$CONFIG_PATH"

time delt-hit library enumerate \
  --config_path=$CONFIG_PATH \
  --counts_path=campaign/selections/AG24_4_counts.txt \
  --top_n=1000 \
  --library_name=AG24_4_top_hits

time delt-hit visualize library --config_path="$CONFIG_PATH" --library_name=AG24_4_top_hits
time delt-hit library properties --config_path="$CONFIG_PATH" --library_name=AG24_4_top_hits

#time delt-hit dashboard \
#  --config_path=$CONFIG_PATH \
#  --counts_path=campaign/selections/AG24_4_counts.txt

time delt-hit analyse enrichment \
  --analysis_config=$ANALYSIS_CONFIG_PATH \
  --name=condition_vs_control \
  --method=counts \
  --save_dir=$ANALYSIS_OUTPUT_ROOT

time Rscript --vanilla "$ANALYSIS_OUTPUT_ROOT/counts/condition_vs_control/enrichment_counts.R"

time delt-hit analyse enrichment \
  --analysis_config=$ANALYSIS_CONFIG_PATH \
  --name=condition_vs_control \
  --method=edgeR \
  --save_dir=$ANALYSIS_OUTPUT_ROOT

time Rscript --vanilla "$ANALYSIS_OUTPUT_ROOT/edgeR/condition_vs_control/enrichment_edgeR.R"

for selection in AG24_13 AG24_14 AG24_15; do
  time delt-hit analyse enrichment \
    --config_path=$CONFIG_PATH \
    --counts=campaign/selections/$selection/counts.txt \
    --method=z_score \
    --name=$selection \
    --save_dir=$ANALYSIS_OUTPUT_ROOT

  time Rscript --vanilla "$ANALYSIS_OUTPUT_ROOT/z_score/$selection/enrichment_z_score.R"
done

# full library enumeration, usually only needed for ML tasks
time delt-hit library enumerate --config_path="$CONFIG_PATH"
time delt-hit library properties --config_path="$CONFIG_PATH"

time delt-hit library represent --method=morgan --config_path="$CONFIG_PATH"
time delt-hit library represent --method=bert --config_path="$CONFIG_PATH"
