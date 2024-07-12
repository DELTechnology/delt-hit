import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel

"""
Pydantic model for configuration.
"""


class Experiment(BaseModel):
    Name: str


class Selection(BaseModel):
    SelectionFile: Path
    FASTQFile: Path
    Library: Path


class StructureItem(BaseModel):
    MaxErrorRate: float
    Indels: bool


class ConfigSimulation(BaseModel):
    OutputFile: Path
    NumReads: int
    Errors: list


class Config(BaseModel):
    Root: Path
    Experiment: Experiment
    Selection: Selection
    Structure: dict[str, StructureItem]
    Simulation: Optional[ConfigSimulation] = None

    def to_dict(self):
        return json.loads(self.model_dump_json())

    def write_yaml(self, path: Path):
        yaml_path = path.with_suffix('.yml')
        with open(yaml_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)
        return yaml_path

    @classmethod
    def from_yaml(cls, path: Path):
        return cls(**read_yaml(path))


def init_config(
        structure: list,
        root: Path = Path.cwd(),
        experiment_name: str = None,
        selection_file: str = 'selections/selection.xlsx',
        fastq_file: str = 'fastq_files/input.fastq.gz',
        library: str = 'libraries/library.xlsx',
        simulation: dict = None,
) -> None:
    experiment_name = get_experiment_name(experiment_name)
    config = {
        'Root': root,
        'Experiment': {
            'Name': experiment_name,
        },
        'Selection': {
            'SelectionFile': selection_file,
            'FASTQFile': fastq_file,
            'Library': library,
        },
        'Structure': {},
    }
    max_error_rate = 0.0
    indels = 0
    for region in structure:
        config['Structure'][region] = {}
        config['Structure'][region]['MaxErrorRate'] = max_error_rate
        config['Structure'][region]['Indels'] = indels
    if simulation:
        config['Simulation'] = simulation
    config_file = Path(root) / 'experiments' / experiment_name / 'config.yml'
    config_file.parent.mkdir(parents=True, exist_ok=True)
    Config(**config).write_yaml(config_file)


def get_experiment_name(
        experiment_name: str = None,
) -> str:
    experiment_name = experiment_name or 'default'
    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    return f'{experiment_name}-{timestamp}'


def is_gz_file(
        path: Path,
) -> bool:
    with open(path, 'rb') as file:
        return file.read(2) == b'\x1f\x8b'


def read_yaml(
        path: Path,
) -> dict:
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def write_yaml(
        data: dict,
        path: Path,
) -> None:
    with open(path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def convert_dict(
        data: dict,
) -> tuple:
    if isinstance(data, dict):
        return tuple((key, convert_dict(value)) for key, value in data.items())
    else:
        return data


def hash_dict(
        data: dict,
) -> str:
    data_tuple = convert_dict(data)
    data_str = json.dumps(data_tuple)
    hash_object = hashlib.sha256()
    hash_object.update(data_str.encode())
    return hash_object.hexdigest()
