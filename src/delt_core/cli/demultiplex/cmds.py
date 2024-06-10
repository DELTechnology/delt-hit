from pathlib import Path
import subprocess

from ... import demultiplex as d


def compute_counts_cli(
        fastq_file: str,
        struct_file: str,
        output_dir: str = None,
) -> None:
    fastq_file = Path(fastq_file).resolve()
    struct_file = Path(struct_file).resolve()
    dir = struct_file.parent
    
    if not output_dir:
        output_dir = dir / 'counts'
        output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = Path(output_dir).resolve()
    
    input_file = dir / 'cutadapt_input_files' / 'demultiplex.sh'
    output_file = dir / 'cutadapt_output_files' / 'reads_with_adapters.gz'
    
    d.generate_input_files(fastq_file, struct_file, input_file, output_file)
    subprocess.run(['bash', input_file])
    d.compute_counts(output_file, output_dir)

