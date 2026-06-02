#!/usr/bin/env bash

set -euo pipefail

cd supporting_material/experiments/example-dual-display || exit

pixi run delt-hit init --excel_path example-dual-display.xlsx

CONFIG_PATH=campaign-dual-display/config.yaml

pixi run delt-hit visualize enumerate --config_path="$CONFIG_PATH"
pixi run delt-hit library enumerate --config_path="$CONFIG_PATH"

# those commands are not supported yet, since they are ill-defined for dual-display campaigns
# pixi run delt-hit visualize library --config_path="$CONFIG_PATH"
# pixi run delt-hit library properties --config_path="$CONFIG_PATH"
