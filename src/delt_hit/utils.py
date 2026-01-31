import hashlib
import json
from datetime import datetime
from pathlib import Path

import yaml


def get_experiment_name(
        experiment_name: str = None,
) -> str:
    """Generate a timestamped experiment name.

    Args:
        experiment_name: Optional base name.

    Returns:
        A name with a timestamp suffix.
    """
    experiment_name = experiment_name or 'default'
    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    return f'{experiment_name}-{timestamp}'


def is_gz_file(
        path: Path,
) -> bool:
    """Check whether a file is gzip-compressed.

    Args:
        path: Path to the file.

    Returns:
        True if the file has a gzip magic header.
    """
    with open(path, 'rb') as file:
        return file.read(2) == b'\x1f\x8b'


def read_yaml(
        path: Path,
) -> dict:
    """Read a YAML file into a dictionary.

    Args:
        path: Path to the YAML file.

    Returns:
        Parsed YAML content.
    """
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def write_yaml(
        data: dict,
        path: Path,
) -> None:
    """Write a dictionary to a YAML file.

    Args:
        data: Data to serialize.
        path: Path to the output YAML file.
    """
    with open(path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def hash_dict(
        data: dict,
) -> str:
    """Hash a dictionary with SHA-256.

    Args:
        data: Dictionary to hash.

    Returns:
        Hex digest string.
    """
    data_str = json.dumps(data, sort_keys=False)
    hash_object = hashlib.sha256()
    hash_object.update(data_str.encode())
    return hash_object.hexdigest()
