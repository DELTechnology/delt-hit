import click

from . import cmds
from pathlib import Path

@click.group()
def qc():
    pass


@qc.command()
@click.argument(
    'experiment_dir', type=click.Path(exists=True, path_type=Path)
)
def report(**kwargs):
    cmds.report(**kwargs)


@qc.command()
@click.argument(
    'experiment_dir', type=click.Path(exists=True, path_type=Path)
)
@click.option(
    '--output_dir',
    '-o',
    default=None,
    type=click.Path(writable=True),
)
def plot(**kwargs):
    cmds.plot(**kwargs)

