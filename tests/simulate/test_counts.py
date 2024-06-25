import os
from pathlib import Path
import tempfile

import pytest

from delt_core.cli.demultiplex.cmds import init as init_demultiplexing
from delt_core.cli.demultiplex.cmds import run as run_demultiplexing
from delt_core.cli.init.cmds import init
from delt_core.cli.simulate.cmds import init as init_simulation
from delt_core.cli.simulate.cmds import run as run_simulation
from delt_core.simulate.utils import read_txt


@pytest.fixture
def load_config():
    return {
        'config_file': 'config.yml',
        'selection_file': 'selections/selection.xlsx',
        'library': 'libraries/library.xlsx',
        'fastq_file': 'fastq_files/input.fastq.gz',
        'output_file': 'fastq_files/input.fastq.gz',
    }


@pytest.fixture(scope='function')
def load_counts():
    def _load_counts():
        output_path = Path.cwd()
        counts_true, counts_pred = [], []
        dirs_true = sorted(list((Path(output_path) / 'evaluations_true').iterdir()))
        dirs_pred = sorted(list((Path(output_path) / 'evaluations').iterdir()))
        for dir_true, dir_pred in zip(dirs_true, dirs_pred):
            file_true = sorted(list(dir_true.iterdir()))[0]
            file_pred = sorted(list(dir_pred.iterdir()))[0]
            counts_true += [read_txt(file_true)]
            counts_pred += [read_txt(file_pred)]
        return counts_true, counts_pred
    yield _load_counts


def test_simulation(load_config, load_counts):
    config = load_config
    dir = Path.cwd()
    with tempfile.TemporaryDirectory(dir=dir) as tmp:
        os.chdir(tmp)
        config['root'] = Path.cwd()
        init()
        init_simulation(**config)
        run_simulation(config['config_file'])
        init_demultiplexing(**config)
        run_demultiplexing(config['config_file'])
        counts_true, counts_pred = load_counts()
        assert counts_true == counts_pred

