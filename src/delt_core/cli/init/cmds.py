from pathlib import Path

import yaml


def init(
        root = None,
        selection_file: Path = 'selection.xlsx',
        fastq_file: Path = 'input.fastq.gz',
        library: Path = 'library.xlsx',
) -> None:
    if not root:
        root = Path.cwd()
    root = Path(root)
    dirs = ['libraries', 'fastq_files', 'selections']
    for dir in dirs:
        path = root / dir
        path.mkdir(exist_ok=True)
    config = {
        'Root': str(root),
        'Selection': {
            'SelectionFile': selection_file,
            'FASTQFile': fastq_file,
            'Library': library,
        },
        'Structure': {},
    }
    max_error_rate = 0.0
    indels = 0
    structure = ['S1', 'C1', 'B1', 'C2', 'B2', 'C3', 'S2']
    for region in structure:
        config['Structure'][region] = {}
        config['Structure'][region]['MaxErrorRate'] = max_error_rate
        config['Structure'][region]['Indels'] = indels
    output_file = root / 'config.yml'
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

