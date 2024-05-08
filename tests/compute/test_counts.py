import json
from pathlib import Path
import typing as tp

import pytest

from delt_core.cli.compute.cmds import compute_counts_cli
from create_simulated_fastq import run_simulation


def read_txt(
        path: str,
) -> tp.Dict:
    with open(path, 'r') as file:
        return file.readlines()


def read_json(
        path: str,
) -> tp.Dict:
    with open(path, 'r') as file:
        return json.load(file)


@pytest.fixture
def load_config():
    if 'exact':
        config_file = '../config/config-no-errors.json'
    # elif 'one_error':
    #     config_file = '../config/config.json'
    config = read_json(config_file)
    return config


@pytest.fixture
def load_counts(load_config):
    config = load_config
    output_path = Path(config['output_file']).parent
    counts_true = read_txt(output_path / 'counts_true.txt')[1:]
    counts_pred = read_txt(output_path / 'counts.txt')[1:]
    return counts_true, counts_pred


def test_simulation(load_config, load_counts):
    config = load_config
    run_simulation(config)
    compute_counts_cli(config['output_file'], config['struct_file'])
    counts_true, counts_pred = load_counts
    assert counts_true == counts_pred


def test_num_codes(load_counts):
    counts_true, counts_pred = load_counts
    assert len(counts_true) == len(counts_pred)


def test_num_matches(load_counts):
    counts_true, counts_pred = load_counts
    sum_true, sum_pred = 0, 0
    for count_true, count_pred in zip(counts_true, counts_pred):
        sum_true += int(count_true.split('\t')[0])
        sum_pred += int(count_pred.split('\t')[0])
    assert sum_true == sum_pred

