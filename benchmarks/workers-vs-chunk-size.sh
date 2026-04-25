cd /Users/adrianomartinelli/projects/delt-hit

# Generate fixtures once
./.venv/bin/python benchmarks/generate_reads_with_adapters.py \
  --output-path benchmarks/data/get_counts/reads_with_adapters_5m.gz \
  --n-reads 5000000

./.venv/bin/python benchmarks/generate_reads_with_adapters.py \
  --output-path benchmarks/data/get_counts/reads_with_adapters_50m.gz \
  --n-reads 50000000

# Run benchmarks sequentially
for input in \
  benchmarks/data/get_counts/reads_with_adapters_5m.gz \
  benchmarks/data/get_counts/reads_with_adapters_50m.gz
do
  for workers in 20 40
  do
    for chunk in 5000000 10000000
    do
      echo
      echo "=== input=${input} workers=${workers} chunk=${chunk} ==="
      ./.venv/bin/python tests/benchmark_get_counts.py \
        --input-path "${input}" \
        --num-workers "${workers}" \
        --chunk-size-bytes "${chunk}"
    done
  done
done
