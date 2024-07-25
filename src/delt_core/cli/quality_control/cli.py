from pathlib import Path

import click

from . import cmds


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


@qc.command()
@click.argument(
    'config_file', type=click.Path(exists=True, path_type=Path)
)
@click.argument(
    'legacy_results_dir', type=click.Path(exists=True, path_type=Path)
)
def compare_with_legacy(**kwargs):
    cmds.compare_with_legacy(**kwargs)


@qc.command()
@click.argument(
    'config_file', type=click.Path(exists=True, path_type=Path)
)
def analyze_codons(**kwargs):
    cmds.analyze_codons(**kwargs)

