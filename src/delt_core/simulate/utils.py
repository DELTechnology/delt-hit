import gzip
from pathlib import Path

import yaml


def read_txt(
        path: Path,
) -> None:
    with open(path, 'r') as file:
        return file.readlines()


def write_txt(
        data: list,
        path: Path,
) -> None:
    if Path(path).suffix == '.gz':
        with gzip.open(path, 'wt') as f:
            f.writelines(data)
    else:
        with open(path, 'w') as f:
            f.writelines(data)


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
        config_file: Path,
        struct_file: Path,
        output_file: Path,
        num_reads: int = 100,
) -> None:
    data = {
        'struct_file': struct_file,
        'output_file': output_file,
        'num_reads': num_reads,
        'errors': [],
    }
    write_yaml(data, config_file)

