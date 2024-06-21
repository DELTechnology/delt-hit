from pathlib import Path

from ... import simulate as s


def init(
        config_simulation: Path,
        config_file: Path,
        output_file: Path,
        num_reads: int,
) -> None:
    s.create_config_file(config_simulation, config_file, output_file, int(num_reads))


def run(
        config_file: Path,
) -> None:
    s.run_simulation(config_file)

