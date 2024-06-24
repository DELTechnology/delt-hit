import os
from pathlib import Path
import shutil
import tempfile

import pytest

from delt_core.cli.init.cmds import init
from delt_core.cli.demultiplex.cmds import run as run_demultiplexing
from delt_core.cli.simulate.cmds import run as run_simulation
from delt_core.simulate.utils import read_yaml, read_txt, write_yaml, create_config_file


@pytest.fixture
def load_config():
    return 0
    config = 'config.yml'
    return read_yaml(config)


@pytest.fixture(scope='function')
def load_counts(load_config):
    def _load_counts():
        config = load_config
        output_path = Path.cwd()
        counts_true, counts_pred = [], []
        dirs_true = sorted(list((Path(output_path) / 'evaluations_true').iterdir()))
        dirs_pred = sorted(list((Path(output_path) / 'evaluations').iterdir()))
        for dir_true, dir_pred in zip(dirs_true, dirs_pred):
            files_true = sorted(list(dir_true.iterdir()))
            files_pred = sorted(list(dir_pred.iterdir()))
            for file_true, file_pred in zip(files_true, files_pred):
                counts_true += [read_txt(file_true)]
                counts_pred += [read_txt(file_pred)]
        return counts_true, counts_pred
    yield _load_counts


def test_simulation(load_config, load_counts):
    config = load_config
    dir = Path.cwd()
    with tempfile.TemporaryDirectory(dir=dir) as tmp:
        # tmp = Path(tmp)
        # for item in tmp.iterdir():
        #     if item.is_dir():
        #         shutil.copytree(item, tmp / item.name)
        #     else:
        #         shutil.copy(item, tmp / item.name)
        # config_file['output_file'] = str(Path(tmp) / Path(config_file['output_file']).name)
        # run_simulation(config_file)
        # run_demultiplexing(config['struct_file'])
        counts_true, counts_pred = load_counts()
        assert counts_true == counts_pred

