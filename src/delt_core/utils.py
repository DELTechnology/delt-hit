import hashlib
import json
from datetime import datetime
from pathlib import Path

import yaml


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


def hash_dict(
        data: dict,
) -> str:
    data_str = json.dumps(data, sort_keys=False)
    hash_object = hashlib.sha256()
    hash_object.update(data_str.encode())
    return hash_object.hexdigest()
