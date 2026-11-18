#!/usr/bin/env bash

set -euo pipefail

# Run from this directory after: conda activate delt-hit

cd supporting_material/experiments/example-dual-display || exit

delt-hit init --excel_path example-dual-display.xlsx

CONFIG_PATH=campaign-dual-display/config.yaml

delt-hit visualize enumerate --config_path="$CONFIG_PATH"
delt-hit library enumerate --config_path="$CONFIG_PATH"

# those commands are not supported yet, since they are ill-defined for dual-display campaigns
# delt-hit visualize library --config_path="$CONFIG_PATH"
# delt-hit library properties --config_path="$CONFIG_PATH"
