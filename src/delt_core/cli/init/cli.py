import click

from . import cmds


@click.command()
def init(**kwargs):
    cmds.init(**kwargs)

