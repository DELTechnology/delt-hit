#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files"
ln -sf "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/20190812.A-1907_NF2GB2_s2_R1.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/0-S0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/0-S0.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/0-S0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/1-C0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/1-C0.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/1-C0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/2-B0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/2-B0.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/2-B0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/3-C1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/3-C1.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/3-C1.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/4-B1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/4-B1.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/4-B1.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/5-C2.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/5-C2.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/5-C2.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/input.fasta.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_input_files/6-S1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/6-S1.cutadapt.json" \
--cores=11 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/6-S1.cutadapt.log"

gzip -cd "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/out.fasta.gz" | awk 'NR % 2 == 1' | gzip -c > "/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli/lane-2-fasta/demultiplex/cutadapt_output_files/reads_with_adapters.gz" || exit
