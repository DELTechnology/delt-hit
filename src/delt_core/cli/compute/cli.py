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
    cmds.compute_smiles_cli(**kwargs)


@compute.command()
@click.argument(
    'input_file',
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
    '--output_file',
    '-o',
    default=None,
    type=click.Path(writable=True),
)
def counts(**kwargs):
    cmds.compute_counts_cli(**kwargs)


@click.group()
def evaluate():
    pass

