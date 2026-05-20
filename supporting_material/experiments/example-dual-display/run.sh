#!/usr/bin/env bash

set -euo pipefail

cd supporting_material/experiments/example-dual-display || exit

delt-hit init --excel_path example-dual-display.xlsx
#delt-hit init --excel_path example-dual-display-extended.xlsx

CONFIG_PATH=campaign/config.yaml
#CONFIG_PATH=campaign-extended/config.yaml

delt-hit visualize enumerate --config_path="$CONFIG_PATH"
delt-hit library enumerate --config_path="$CONFIG_PATH"
