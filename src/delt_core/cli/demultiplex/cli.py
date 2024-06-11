import click

from . import cmds
from pathlib import Path

@click.group()
def demultiplex():
    pass


@demultiplex.command()
@click.argument('path_to_struct_file', type=click.Path(exists=True, path_type=Path))
def convert(path_to_struct_file: Path):
    cmds.convert(path_to_struct_file)


@demultiplex.command()
# @click.argument()
def init(**kwargs):
    cmds.init(**kwargs)


@demultiplex.command()
# @click.argument()
def create_lists(**kwargs):
    cmds.create_lists(**kwargs)


@demultiplex.command()
@click.argument(
    'struct_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    '--fastq_file',
    '-f',
    default=None,
    type=click.Path(writable=True),
)
def create_cutadapt_input(**kwargs):
    cmds.create_cutadapt_input(**kwargs)


@demultiplex.command()
@click.argument(
    'input_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
@click.argument(
    'output_dir',
    nargs=1,
    required=True,
    type=click.Path(exists=False),
)
def compute_counts(**kwargs):
    cmds.compute_counts(**kwargs)


@demultiplex.command()
@click.argument(
    'struct_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    '--fastq_file',
    '-f',
    default=None,
    type=click.Path(writable=True),
)
def run(**kwargs):
    cmds.run(**kwargs)

