import click

from . import cmds


@click.group()
def demultiplex():
    pass


@demultiplex.command()
@click.argument(
    'fastq_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
@click.argument(
    'struct_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    '--output_dir',
    '-o',
    default=None,
    type=click.Path(writable=True),
)
def counts(**kwargs):
    cmds.compute_counts_cli(**kwargs)

