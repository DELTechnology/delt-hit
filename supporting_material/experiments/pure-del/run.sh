#!/usr/bin/env bash

cd /Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del || exit

PREFIX=lane-1
#PREFIX=lane-2

delt-hit init --excel_path $PREFIX.xlsx

CONFIG_PATH=$PREFIX/config.yaml

delt-hit demultiplex prepare --config_path=$CONFIG_PATH
$PREFIX/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=$CONFIG_PATH
delt-hit demultiplex qc --config_path=$CONFIG_PATH

delt-hit demultiplex process --config_path=$CONFIG_PATH
