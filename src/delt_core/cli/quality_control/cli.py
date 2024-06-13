import click

from . import cmds


@click.group()
def qc():
    pass


@qc.command()
@click.argument(
    'input_path',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def report(**kwargs):
    cmds.report(**kwargs)


@qc.command()
@click.argument(
    'input_dir',
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
def plot(**kwargs):
    cmds.plot(**kwargs)

