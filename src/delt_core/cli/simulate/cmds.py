from pathlib import Path

from ... import simulate as s


def init(
        config_file: Path,
        struct_file: Path,
        output_file: Path = 'simulation.fastq.gz',
        num_reads: int = 100,
) -> None:
    s.create_config_file(config_file, struct_file, output_file, int(num_reads))


def run(
        config_file: Path,
) -> None:
    s.run_simulation(config_file)

