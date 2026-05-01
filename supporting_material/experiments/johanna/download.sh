#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit

if [[ -s "407991_1-2601_JPAG_L1c_NF_501701_S3_R1_001.fastq.gz" ]]; then
  echo "407991_1-2601_JPAG_L1c_NF_501701_S3_R1_001.fastq.gz already exists, skipping"
  exit 0
fi

echo "Downloading 407991_1-2601_JPAG_L1c_NF_501701_S3_R1_001.fastq.gz"
curl -L --fail --continue-at - \
  --output "407991_1-2601_JPAG_L1c_NF_501701_S3_R1_001.fastq.gz" \
  "https://polybox.ethz.ch/index.php/s/bXTtkZ2666eZoDT/download"
