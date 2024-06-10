import click
from pathlib import Path

@click.command()
@click.option('--from-fastq', type=click.Path(exists=True, path_type=Path))
@click.option('--from-old-structure', type=click.Path(exists=True, path_type=Path))
@click.option('--from-library', type=click.Path(exists=True, path_type=Path))
def init(from_fastq: Path = None, from_old_structure: Path = None, from_library: Path = None):
    if from_fastq:
        struct = get_default_struct(from_fastq)
    elif from_old_structure:
        struct = read_struct(from_old_structure)
    elif from_library:
        struct = read_struct_from_library(from_library)
    else:
        raise ValueError('Provide at least one of the following options: --from-fastq, --from-old-structure, --from-library')


def run():
    pass
