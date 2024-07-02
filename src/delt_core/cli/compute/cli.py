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
@click.option(
    '--output_path',
    '-o',
    default=None,
    type=click.Path(writable=True),
)
def smiles(**kwargs):
    cmds.compute_smiles(**kwargs)


@compute.command()
@click.argument(
    'input_path',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def properties(**kwargs):
    cmds.compute_properties(**kwargs)

