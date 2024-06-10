from pathlib import Path
import subprocess
import typing as tp

from ... import compute as c


def compute_smiles_cli(
        input_path: tp.Tuple,
        output_path: str = None,
) -> None:
    if not output_path:
        output_path = Path(input_path[0]).parent / 'smiles.txt.gz'

    libraries = []
    for library in input_path:
        libraries += [c.load_data(library)]
    
    c.compute_smiles(libraries, Path(output_path))


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
    
    c.generate_input_files(fastq_file, struct_file, input_file, output_file)
    subprocess.run(['bash', input_file])
    c.compute_counts(output_file, output_dir)


def evaluate_cli():
    pass

