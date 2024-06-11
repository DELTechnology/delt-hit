from pathlib import Path
import subprocess

from ... import demultiplex as d


def convert(path_to_struct_file: Path) -> None:
    d.migrate_to_new_structure(path_to_struct_file)


def init() -> None:
    pass


def create_lists() -> None:
    pass


def create_cutadapt_input(
        struct_file: Path,
        fastq_file: Path = None,
) -> None:
    struct_file = Path(struct_file).resolve()
    if not fastq_file:
        fastq_file = struct_file.parent / 'input.fastq.gz'
    fastq_file = Path(fastq_file).resolve()
    d.generate_input_files(struct_file, fastq_file)


def compute_counts(
        input_file: Path,
        output_dir: Path,
) -> None:
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    d.compute_counts(input_file, output_dir)


def run(
        struct_file: Path,
        fastq_file: Path = None,
) -> None:
    struct_file = Path(struct_file).resolve()
    create_cutadapt_input(struct_file, fastq_file)
    input_file = struct_file.parent / 'cutadapt_input_files' / 'demultiplex.sh'
    subprocess.run(['bash', input_file])

