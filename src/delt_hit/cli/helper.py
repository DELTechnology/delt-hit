from pathlib import Path

from delt_hit.utils import read_yaml


def get_experiment_dir(config_path: Path) -> Path:
    """Resolve the experiment output directory from a config file."""
    cfg = read_yaml(config_path)
    return Path(cfg["experiment"]["save_dir"]).expanduser().resolve() / cfg["experiment"]["name"]


def get_library_dir(config_path: Path) -> Path:
    """Resolve the experiment library output directory."""
    return get_experiment_dir(config_path) / "library"


def get_library_path(config_path: Path) -> Path:
    """Resolve the default library parquet path."""
    return get_library_dir(config_path) / "library.parquet"


def get_named_library_path(config_path: Path, library_name: str) -> Path:
    """Resolve a named library parquet path inside the experiment library dir."""
    return get_library_dir(config_path) / f"{library_name}.parquet"
