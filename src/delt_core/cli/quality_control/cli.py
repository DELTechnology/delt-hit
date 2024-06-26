import click

from . import cmds


@click.group()
def qc():
    pass


@qc.command()
def report(**kwargs):
    cmds.report(**kwargs)


@qc.command()
@click.option(
    '--output_dir',
    '-o',
    default=None,
    type=click.Path(writable=True),
)
def plot(**kwargs):
    cmds.plot(**kwargs)

