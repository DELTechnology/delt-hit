import json
from pathlib import Path
import typing as tp


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


def write_json(
        data: tp.Dict,
        path: str,
) -> None:
    with open(path, 'w') as file:
        json.dump(data, file)


def generate_struct_file(
        path: str,
) -> None:
    
    path = Path(path)
    info = read_txt(path / 'structure.txt')[2:]
    info = sorted(info, key=lambda x: int(x.split('\t')[0]))

    count = {}
    structure = {}

    for line in info:
        _, _, component, file = line.split('\t')
        count[component] = count.get(component, 0) + 1
        component += str(count[component])
        structure[component] = {
            'Sequences': [],
            'Maximum Error Rate': 0.1,
            'Minimum Overlap': 3,
        }

        sequences = read_txt(path / file.strip())
        for sequence in sequences:
            structure[component]['Sequences'] += [sequence.strip()]

    write_json(structure, path / 'structure.json')



if __name__ == '__main__':

    path = 'data/sequences'
    generate_struct_file(path)

