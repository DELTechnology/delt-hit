import click

from . import cmds


@click.group()
def demultiplex():
    pass


@demultiplex.command()
@click.argument(
    'struct_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def convert(**kwargs):
    cmds.convert(**kwargs)


@demultiplex.command()
@click.argument(
    'config_file',
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
def create_lists(**kwargs):
    cmds.create_lists(**kwargs)


@demultiplex.command()
@click.argument(
    'config_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def create_cutadapt_input(**kwargs):
    cmds.create_cutadapt_input(**kwargs)


@demultiplex.command()
@click.argument(
    'config_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
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
    'config_file',
    nargs=1,
    required=True,
    type=click.Path(exists=True),
)
def run(**kwargs):
    cmds.run(**kwargs)

