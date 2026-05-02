#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit

PREFIX="${1:-lane-1}"

analysis_lane() {
  case "$1" in
    his_pure_up|dyna_up)
      echo "lane-2"
      ;;
    his_pure_sp|dyna_sp)
      echo "lane-1"
      ;;
    *)
      echo "Unknown analysis: $1" >&2
      return 1
      ;;
  esac
}

pixi run delt-hit init --excel_path "$PREFIX.xlsx"

CONFIG_PATH=$PREFIX/config.yaml

pixi run delt-hit demultiplex prepare --config_path=$CONFIG_PATH
$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh

pixi run delt-hit demultiplex report --config_path=$CONFIG_PATH
pixi run delt-hit demultiplex qc --config_path=$CONFIG_PATH

pixi run delt-hit demultiplex process --config_path=$CONFIG_PATH
pixi run delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True

if [[ ! -d lane-1/selections || ! -d lane-2/selections ]]; then
  echo "Skipping enrichment: processed selections are required for both lane-1 and lane-2." >&2
  echo "Run 'bash run.sh lane-1' and 'bash run.sh lane-2' first, then rerun either command to generate analysis outputs." >&2
  exit 0
fi

for ANALYSIS in his_pure_up his_pure_sp dyna_up dyna_sp
do
  ANALYSIS_LANE="$(analysis_lane "$ANALYSIS")"
  ANALYSIS_DIR="$(pwd)/$ANALYSIS_LANE/analysis/$ANALYSIS"

  pixi run delt-hit analyse enrichment \
    --config_path=analysis.yaml \
    --name=$ANALYSIS \
    --method=counts

  Rscript --vanilla "$ANALYSIS_DIR/counts/enrichment_counts.R"

  pixi run delt-hit analyse enrichment \
    --config_path=analysis.yaml \
    --name=$ANALYSIS \
    --method=edgeR

  Rscript --vanilla "$ANALYSIS_DIR/edgeR/enrichment_edgeR.R"

  pixi run python enrichment.py \
    --data-dir "$ANALYSIS_DIR" \
    --output-dir "$(pwd)/enrichment/$ANALYSIS"
done
