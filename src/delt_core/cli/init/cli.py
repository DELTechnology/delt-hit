import click

from . import cmds


@click.command()
@click.option(
    '--root',
    '-r',
    default=None,
    type=click.Path(writable=True),
)
def init(**kwargs):
    cmds.init(**kwargs)

