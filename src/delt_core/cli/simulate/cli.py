import click

from . import cmds


@click.group()
def simulate():
    pass


@simulate.command()
@click.argument(
    'config_file',
    nargs=-1,
    required=True,
    type=click.Path(exists=True),
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
# @click.option(
#     '--output_dir',
#     '-o',
#     default=None,
#     type=click.Path(writable=True),
# )
def run(**kwargs):
    cmds.run(**kwargs)

