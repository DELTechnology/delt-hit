#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit

delt-hit init --excel_path example-dual-display.xlsx

CONFIG_PATH=campaign/config.yaml

delt-hit visualize enumerate --config_path="$CONFIG_PATH"
delt-hit library enumerate --config_path="$CONFIG_PATH"
