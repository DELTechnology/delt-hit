#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files"
ln -sf "/Users/adrianomartinelli/projects/delt-hit/benchmarks/data/synthetic_4cycle_100m/synthetic_4cycle_100m.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz"

mv "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_input_files/0-S0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/0-S0.cutadapt.json" \
--cores=1 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/0-S0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_input_files/1-B0.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/1-B0.cutadapt.json" \
--cores=1 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/1-B0.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_input_files/2-B1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/2-B1.cutadapt.json" \
--cores=1 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/2-B1.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_input_files/3-B2.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/3-B2.cutadapt.json" \
--cores=1 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/3-B2.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_input_files/4-B3.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/4-B3.cutadapt.json" \
--cores=1 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/4-B3.cutadapt.log"

mv "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz"

cutadapt "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/input.fastq.gz" \
-o "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" \
-e 0.0 --no-indels \
-g "^file:/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_input_files/5-S1.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json="/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/5-S1.cutadapt.json" \
--cores=1 2>&1 | tee "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/5-S1.cutadapt.log"

zgrep @ "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/out.fastq.gz" | gzip -c > "/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m/synthetic_4cycle_100m/demultiplex/cutadapt_output_files/reads_with_adapters.gz" || exit
