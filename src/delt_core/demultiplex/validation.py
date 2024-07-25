import json
from pathlib import Path
from typing import Optional

import pandas as pd
from pydantic import BaseModel, computed_field, ValidationError
import yaml

from .utils import read_yaml, get_experiment_name


class Region(BaseModel):
    name: str
    codons: list[str]
    max_error_rate: float
    indels: int
    position_in_construct: int = None

    @computed_field
    @property
    def region_id(self) -> str:
        return f'{self.position_in_construct}-{self.name}'


class SelectionFile(BaseModel):
    SelectionID: int
    Library: str
    FwdPrimer1: str
    RevPrimer1: str
    FASTQFile: str


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


def validate(
        data: pd.DataFrame,
        model: BaseModel,
) -> dict:
    records = data.to_dict(orient='records')
    for record in records:
        try:
            model(**record)
        except ValidationError as e:
            print(f'Validation error: {e.json()}')
    return data


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


