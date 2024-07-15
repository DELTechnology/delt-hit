import click

from . import cmds


@click.group()
def simulate():
    pass


@simulate.command()
@click.option(
    '--root',
    '-r',
    default=None,
    type=click.Path(writable=True),
)
@click.option(
    '--experiment_name',
    '-e',
    default='',
    type=str,
)
@click.option(
    '--selection_file',
    '-s',
    default='selections/selection_template.xlsx',
    type=click.Path(writable=True),
)
@click.option(
    '--library',
    '-l',
    default='libraries/library_template.xlsx',
    type=click.Path(writable=True),
)
@click.option(
    '--fastq_file',
    '-f',
    default='fastq_files/input.fastq.gz',
    type=click.Path(writable=True),
)
@click.option(
    '--output_file',
    '-o',
    default='fastq_files/simulation.fastq.gz',
    type=click.Path(writable=True),
)
@click.option(
    '--num_reads',
    '-n',
    default=100,
    type=int,
)
def init(**kwargs):
    cmds.init(**kwargs)


@simulate.command()
@click.argument(
    'config_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def run(**kwargs):
    cmds.run(**kwargs)

