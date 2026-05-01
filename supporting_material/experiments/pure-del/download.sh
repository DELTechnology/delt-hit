#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit

download() {
  local url="$1"
  local out="$2"

  if [[ -s "$out" ]]; then
    echo "$out already exists, skipping"
    return
  fi

  echo "Downloading $out"
  curl -L --fail --continue-at - --output "$out" "$url"
}

download \
  "https://polybox.ethz.ch/index.php/s/H8YB6wdxpwWgzLB/download" \
  "328321_1-230906_MK_TSD501701_S1_R1_001_SP_fl1.fastq.gz"

download \
  "https://polybox.ethz.ch/index.php/s/wK6t8qjDTay8rAg/download" \
  "328321_2-230906_MK_TSD504704_S2_R1_001_UP_fl2.fastq.gz"

ZIP_PATH="Keller_Petrov_etal_Supplementary_Material2_Sequencing_Data.zip"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

download \
  "https://zenodo.org/records/11122901/files/Keller_Petrov_etal_Supplementary_Material2_Sequencing_Data.zip?download=1" \
  "$ZIP_PATH"

echo "Extracting published selection tables"
unzip -oq "$ZIP_PATH" -d "$TMP_DIR"

EXTRACTED_DIR="$(find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d ! -name '__MACOSX' | head -n 1)"

if [[ -z "$EXTRACTED_DIR" ]]; then
  echo "Could not locate extracted Zenodo directory" >&2
  exit 1
fi

mkdir -p published
find "$EXTRACTED_DIR" -maxdepth 1 -type f -name 'selection_*_.txt' -exec cp -f {} published/ \;
