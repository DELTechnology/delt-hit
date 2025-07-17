#!/bin/bash

LOG_DIR="/work/FAC/FBM/DBC/mrapsoma/prometex/logs/adrianom"
mkdir -p "$LOG_DIR"

for EVAL_DIR in /work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/experiments/eval-*; do
  if [ -f "$EVAL_DIR/cutadapt_output_files/reads_with_adapters.gz" ]; then
    echo "Skipping $EVAL_DIR: reads_with_adapters.gz already exists."
    continue
  fi

  BASENAME=$(basename "$EVAL_DIR")
  echo "Submitting job for $EVAL_DIR"

  sbatch --job-name=delt-demultiplex \
         --cpus-per-task=48 \
         --mem=64G \
         --time=04:00:00 \
         --output="$LOG_DIR/${BASENAME}-%j.out" \
         --error="$LOG_DIR/${BASENAME}-%j.err" \
         <<EOF
#!/bin/bash
module load python/3.12.1
source /work/FAC/FBM/DBC/mrapsoma/prometex/envs/adrianom/delt-core/bin/activate

EVAL_DIR="${EVAL_DIR}"

echo "Processing \$EVAL_DIR"
delt-cli demultiplex run "\$EVAL_DIR/config.yml" > "\$EVAL_DIR/output.log" 2>&1
delt-cli qc report "\$EVAL_DIR" > "\$EVAL_DIR/qc_report.txt" 2>&1
delt-cli qc plot "\$EVAL_DIR"
cat "\$EVAL_DIR/qc_report.txt"
EOF

done