#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files"
ln -sf "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/20190812.A-1907_NF2GB2_s1_R1.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/0-S0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/0-S0.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/0-S0.cutadapt.log"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/1-C0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/1-C0.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/1-C0.cutadapt.log"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/2-B0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/2-B0.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/2-B0.cutadapt.log"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/3-C1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/3-C1.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/3-C1.cutadapt.log"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/4-B1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/4-B1.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/4-B1.cutadapt.log"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/5-C2.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/5-C2.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/5-C2.cutadapt.log"

mv "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_input_files/6-S1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/6-S1.cutadapt.json" \
--cores=10 2>&1 | tee "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/6-S1.cutadapt.log"

zgrep @ "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/out.fastq.gz" | gzip -c > "/Users/jorsch/Desktop/Nicholas_library_fastq/1907_NF2GB2MC/NF2_sels1907/DELT-hit/Evaluation_Delt-hit_NatChem_fl1/NF2_NatChem_flowlane1/demultiplex/cutadapt_output_files/reads_with_adapters.gz" || exit
