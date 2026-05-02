#!/usr/bin/env bash

# Example Slurm submissions:
# sbatch --job-name=favalli-lane-1-fasta --mem=32G --time=04:00:00 --cpus-per-task=12 --output="$HOME/logs/%j.out" --wrap="bash supporting_material/experiments/favalli/run.sh lane-1-fasta"
# sbatch --job-name=favalli-lane-2-fasta --mem=32G --time=04:00:00 --cpus-per-task=12 --output="$HOME/logs/%j.out" --wrap="bash supporting_material/experiments/favalli/run.sh lane-2-fasta"

set -euo pipefail

cd supporting_material/experiments/favalli || exit

PREFIX="${1:-lane-1}"
VALID_PREFIXES=("lane-1" "lane-2" "lane-1-fasta" "lane-2-fasta")

if [[ ! " ${VALID_PREFIXES[*]} " =~ [[:space:]]${PREFIX}[[:space:]] ]]; then
  echo "Unsupported prefix '${PREFIX}'. Expected one of: ${VALID_PREFIXES[*]}" >&2
  exit 1
fi

delt-hit init --excel_path "$PREFIX.xlsx"

CONFIG_PATH=$PREFIX/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
delt-hit demultiplex process --config_path=$CONFIG_PATH --as_files=True

# TODO: check with Jörg what selections to compare
delt-hit analyse enrichment \
  --config_path=analysis.yaml \
  --name=protein_vs_no_protein \
  --method=counts

Rscript --vanilla "$PREFIX/analysis/protein_vs_no_protein/counts/enrichment_counts.R"

delt-hit analyse enrichment \
  --config_path=analysis.yaml \
  --name=protein_vs_no_protein \
  --method=edgeR

Rscript --vanilla "$PREFIX/analysis/protein_vs_no_protein/edgeR/enrichment_edgeR.R"
