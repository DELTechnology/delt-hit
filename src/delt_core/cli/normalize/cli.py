import click

from . import cmds


@click.group()
def normalize():
    pass


@normalize.command()
@click.argument(
    'config_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
@click.argument(
    'target',
    nargs=1,
    required=True,
    type=str,
)
@click.argument(
    'control',
    nargs=1,
    required=True,
    type=str,
)
def run(**kwargs):
    cmds.run(**kwargs)

