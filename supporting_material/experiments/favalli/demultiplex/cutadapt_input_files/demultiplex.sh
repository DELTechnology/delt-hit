#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files"
ln -sf "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/favalli.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/0-S0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/0-S0.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/0-S0.cutadapt.log"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/1-C0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/1-C0.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/1-C0.cutadapt.log"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/2-B0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/2-B0.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/2-B0.cutadapt.log"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/3-C1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/3-C1.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/3-C1.cutadapt.log"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/4-B1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/4-B1.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/4-B1.cutadapt.log"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/5-C2.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/5-C2.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/5-C2.cutadapt.log"

mv "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_input_files/6-S1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/6-S1.cutadapt.json" \
--cores=11 2>&1 | tee "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/6-S1.cutadapt.log"

zcat "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/out.fastq.gz" | sed -n '1~4p' | gzip -c > "/work/FAC/FBM/DBC/mrapsoma/prometex/projects/delt-hit/supporting_material/experiments/favalli/demultiplex/cutadapt_output_files/reads_with_adapters.gz" || exit
