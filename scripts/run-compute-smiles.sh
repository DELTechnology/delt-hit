#!/bin/bash

LOG_DIR="/work/FAC/FBM/DBC/mrapsoma/prometex/logs/adrianom"
mkdir -p "$LOG_DIR"

LIBRARIES_DIR="/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries"
SMILES_DIR="$LIBRARIES_DIR/smiles"
mkdir -p "$SMILES_DIR"

LIBRARY_NAMES=('NF.xlsx' 'GB2.xlsx' 'NF-KS_corrected_codes3RC.xlsx')

for LIBRARY in "${LIBRARY_NAMES[@]}"; do
  BASENAME=$(basename "$LIBRARY" .xlsx)
  SAVE_PATH="$SMILES_DIR/${BASENAME}_smiles.txt.gz"

  if [ -f "$SAVE_PATH" ]; then
    echo "Skipping $SAVE_PATH: Smiles are already computed."
    continue
  fi

  LIBRARY_PATH="$LIBRARIES_DIR/$LIBRARY"
  echo "Submitting job for $LIBRARY_PATH"

  sbatch --job-name=delt-demultiplex \
         --cpus-per-task=32 \
         --mem=64G \
         --time=04:00:00 \
         --output="$LOG_DIR/${BASENAME}-%j.out" \
         --error="$LOG_DIR/${BASENAME}-%j.err" \
         <<EOF
#!/bin/bash
module load python/3.12.1
source /work/FAC/FBM/DBC/mrapsoma/prometex/envs/adrianom/delt-core/bin/activate

echo "Processing ${LIBRARY_PATH}"
delt-cli compute smiles "${LIBRARY_PATH}"
EOF

done