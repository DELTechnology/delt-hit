import click

from . import cmds


@click.group()
def simulate():
    pass


@simulate.command()
@click.option(
    '--config_simulation',
    '-s',
    default='config_simulation.yml',
    type=click.Path(writable=True),
)
@click.option(
    '--config_file',
    '-c',
    default='config.yml',
    type=click.Path(writable=True),
)
@click.option(
    '--output_file',
    '-o',
    default='simulation.fastq',
    type=click.Path(writable=True),
)
@click.option(
    '--num_reads',
    '-n',
    default=100,
    type=click.Path(writable=True),
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

