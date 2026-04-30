#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files"
ln -sf "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/328321_1-230906_MK_TSD501701_S1_R1_001_SP_fl1.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/0-S0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/0-S0.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/0-S0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/1-C0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/1-C0.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/1-C0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/2-B0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/2-B0.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/2-B0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/3-C1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/3-C1.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/3-C1.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/4-B1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/4-B1.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/4-B1.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/5-C2.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/5-C2.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/5-C2.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/6-B2.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/6-B2.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/6-B2.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/7-C3.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/7-C3.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/7-C3.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_input_files/8-S1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/8-S1.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/8-S1.cutadapt.log"

gzip -cd "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/out.fastq.gz" | awk 'NR % 4 == 1' | gzip -c > "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/pure-del/lane-1/demultiplex/cutadapt_output_files/reads_with_adapters.gz" || exit
