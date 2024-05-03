import json
from pathlib import Path
import typing as tp


def read_json(
        path: str,
) -> tp.Dict:
    with open(path, 'r') as file:
        return json.load(file)


def main():

    path = Path('../../data/simulation')

    counts_true = read_json(path / 'counts_true.json')
    counts_pred = read_json(path / 'counts.json')

    print(counts_true == counts_pred)


if __name__ == '__main__':

    main()

