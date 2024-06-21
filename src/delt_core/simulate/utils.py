from pathlib import Path

import yaml


def read_txt(
        path: Path,
) -> None:
    with open(path, 'r') as file:
        return file.readlines()


def read_yaml(
        path: Path,
) -> None:
    with open(path, 'r') as file:
        return yaml.safe_load(file)


def write_yaml(
        data: dict,
        path: Path,
) -> None:
    with open(path, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)


def create_config_file(
        config_simulation: Path,
        config_file: Path,
        output_file: Path,
        num_reads: int = 100,
) -> None:
    data = {
        'config_file': config_file,
        'output_file': output_file,
        'num_reads': num_reads,
        'errors': [],
    }
    write_yaml(data, config_simulation)

