from pathlib import Path

import yaml


def init() -> None:
    root = Path.cwd()
    dirs = ['libraries', 'fastq_files', 'selections']
    for dir in dirs:
        path = root / dir
        path.mkdir(exist_ok=True)
    config = {
        'Root': str(root),
        'Selection': {
            'SelectionFile': 'selection.xlsx',
            'FASTQFile': 'input.fastq.gz',
            'Library': 'library.xlsx'
        },
        'Structure': {},
    }
    max_error_rate = 0.0
    indels = 0
    structure = ['S1', 'C1', 'B1', 'C2', 'S2']
    for region in structure:
        config['Structure'][region] = {}
        config['Structure'][region]['MaxErrorRate'] = max_error_rate
        config['Structure'][region]['Indels'] = indels
    output_file = root / 'config.yml'
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

