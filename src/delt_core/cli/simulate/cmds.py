from pathlib import Path

from ... import simulate as s


def init(
        *,
        root: Path,
        config_file: Path,
        experiment_name: str,
        selection_file: Path,
        library: Path,
        fastq_file: Path,
        output_file: Path,
        num_reads: int = 100,
) -> None:
    if not Path(library).exists():
        s.create_library_template(root, library)
    if not Path(selection_file).exists():
        s.create_selection_template(root, selection_file, library, fastq_file)
    s.create_config_file(
        root=root,
        config_file=config_file,
        experiment_name=experiment_name,
        selection_file=selection_file,
        fastq_file=fastq_file,
        library=library,
        output_file=output_file,
        num_reads=num_reads,
    )


def run(
        config_file: Path,
) -> None:
    s.run_simulation(config_file)

