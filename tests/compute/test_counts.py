from pathlib import Path
import tempfile

import pytest

from delt_core.cli.compute.cmds import compute_counts_cli
from run_simulation import run_simulation
import utils


@pytest.fixture
def load_config():
    if 'exact':
        config_file = 'config/config-no-errors.json'
    # elif 'one_error':
    #     config_file = '../config/config.json'
    config = utils.read_json(config_file)
    return config


@pytest.fixture(scope='function')
def load_counts(load_config):
    def _load_counts():
        config = load_config
        output_path = Path(config['output_file']).parent
        counts_true, counts_pred = [], []
        files_true = sorted(list(Path(output_path / 'counts_true').iterdir()))
        files_pred = sorted(list(Path(output_path / 'counts').iterdir()))
        for file_true, file_pred in zip(files_true, files_pred):
            counts_true += [utils.read_txt(file_true)]
            counts_pred += [utils.read_txt(file_pred)]
        return counts_true, counts_pred
    yield _load_counts


def test_simulation(load_config, load_counts):
    config = load_config
    dir = Path(config['output_file']).parent
    with tempfile.TemporaryDirectory(dir=dir) as tmp:
        config['output_file'] = Path(tmp) / Path(config['output_file']).name
        run_simulation(config)
        compute_counts_cli(config['output_file'], config['struct_file'])
        counts_true, counts_pred = load_counts()
        assert counts_true == counts_pred

