import click

from . import cmds


@click.group()
def compute():
    pass


@compute.command()
@click.argument(
    'input_path',
    nargs=-1,
    required=True,
    type=click.Path(exists=True),
)
def smiles(**kwargs):
    cmds.compute_smiles(**kwargs)


@compute.command()
@click.argument(
    'input_path',
    nargs=2,
    required=True,
    type=click.Path(exists=True),
)
def merge(**kwargs):
    cmds.merge_libraries(**kwargs)


@compute.command()
@click.argument(
    'input_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def properties(**kwargs):
    cmds.compute_properties(**kwargs)


@compute.command()
@click.argument(
    'input_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def plot(**kwargs):
    cmds.plot(**kwargs)

