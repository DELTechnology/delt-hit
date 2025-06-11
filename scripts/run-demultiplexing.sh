#!/bin/bash

for EVAL_DIR in /work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/experiments/eval-*; do
  if [ -f "$EVAL_DIR/cutadapt_output_files/reads_with_adapters.gz" ]; then
    echo "Skipping $EVAL_DIR: reads_with_adapters.gz already exists."
    continue
  fi

  echo "Processing $EVAL_DIR"
  delt-cli demultiplex run "$EVAL_DIR"/config.yml > "$EVAL_DIR/output.log" 2>&1
  delt-cli qc report "$EVAL_DIR" > "$EVAL_DIR/qc_report.txt" 2>&1
  delt-cli qc plot "$EVAL_DIR"
  cat "$EVAL_DIR/qc_report.txt"
done