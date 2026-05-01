#!/usr/bin/env bash

cd /Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/template || exit
delt-hit init --excel_path template-dual-display.xlsx

CONFIG_PATH=campaign-dual/config.yaml

delt-hit visualize enumerate --config_path=$CONFIG_PATH
delt-hit library enumerate --config_path=$CONFIG_PATH
delt-hit library properties --config_path=$CONFIG_PATH

delt-hit library represent --method=morgan --config_path=$CONFIG_PATH
delt-hit library represent --method=bert --config_path=$CONFIG_PATH
