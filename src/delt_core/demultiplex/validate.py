import yaml
from pathlib import Path

def validate_config(path):
    config = yaml.safe_load(open(path))
    assert Path(config['experiment']['fastq_path']).exists()

    # TODO
    # check that all reactions are defined
    # check that all scaffolds are defined