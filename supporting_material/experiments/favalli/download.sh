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
  "https://polybox.ethz.ch/index.php/s/HbtDGTznwWGxwKr/download" \
  "20190812.A-1907_NF2GB2_s1_R1.fastq.gz"

download \
  "https://polybox.ethz.ch/index.php/s/9bX3YnzyDy8ez8g/download" \
  "20190812.A-1907_NF2GB2_s2_R1.fastq.gz"

mkdir -p published

download \
  "https://polybox.ethz.ch/index.php/s/cgxBWkHnpkT6rBi/download" \
  "published/1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt"
