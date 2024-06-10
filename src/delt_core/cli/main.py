import click

from . import compute as c
from . import demultiplex as d


@click.group()
def cli():
    pass


cli.add_command(c.compute)
cli.add_command(d.demultiplex)


if __name__ == '__main__':
    cli()
