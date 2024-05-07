import json
from pathlib import Path
import typing as tp

import pytest


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
def load_data():
    path = Path('../../data/simulation')
    counts_true = read_txt(path / 'counts_true.txt')[1:]
    counts_pred = read_txt(path / 'counts.txt')[1:]
    return counts_true, counts_pred


def test_num_codes(load_data):
    counts_true, counts_pred = load_data
    assert len(counts_true) == len(counts_pred)


def test_num_matches(load_data):
    counts_true, counts_pred = load_data
    sum_true, sum_pred = 0, 0
    for count_true, count_pred in zip(counts_true, counts_pred):
        sum_true += int(count_true.split('\t')[0])
        sum_pred += int(count_pred.split('\t')[0])
    assert sum_true == sum_pred


def test_match(load_data):
    counts_true, counts_pred = load_data
    assert counts_true == counts_pred



# counts_
# @pytest.mark.parametrize("true, pred", counts_true)
# def test_compare_elements(true, pred):
#     assert true == pred, f"Difference found: {true} != {pred}"
    


def main():

    path = Path('../../data/simulation')

    counts_true = read_txt(path / 'counts_true.txt')[1:]
    counts_pred = read_txt(path / 'counts.txt')[1:]

    structure = read_json(path / 'structure.json')
    codes1 = structure['B1']['Sequences']
    codes2 = structure['B2']['Sequences']

    reads = read_txt(path / 'simulation.fastq')
    # print(len(reads))
    # exit()

    errors = []
    num_reads = len(counts_pred)

    for i, (count_true, count_pred) in enumerate(zip(counts_true, counts_pred)):
        if count_true != count_pred:
            code1_true, code2_true = count_true.split('\t')[1:]
            code1_pred, code2_pred = count_pred.split('\t')[1:]
            errors += [i]
            print(count_true.strip())
            print(count_pred.strip())
            print('True:', codes1[int(code1_true)], codes1[int(code2_true)])
            print('Pred:', codes1[int(code1_pred)], codes1[int(code2_pred)])
            print('\n')


    print(f'{num_reads - len(errors)} / {num_reads}')
    print(errors)


if __name__ == '__main__':

    main()

