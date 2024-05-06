import json
from pathlib import Path
import typing as tp

import numpy as np


def read_json(
        path: str,
) -> tp.Dict:
    with open(path, 'r') as file:
        return json.load(file)


def main():

    path = Path('../../data/simulation')

    counts_true = read_json(path / 'counts_true.json')
    counts_pred = read_json(path / 'counts.json')

    match = counts_true == counts_pred
    print(match)

    if not match:
        count = 0

        for key, _ in counts_true.items():
            if key in counts_pred:
                count += counts_pred[key]
        
        print(f'{count} / {sum(counts_true.values())}')
    
    # indices_true = np.load(path / 'indices_true.npy')
    # indices_pred = np.load(path / 'indices.npy')

    # print(indices_true)
    # print(indices_pred)
    # print((indices_true == indices_pred).all())


if __name__ == '__main__':

    main()

