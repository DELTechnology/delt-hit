import click

from . import compute as c
from . import quality_control as q


@click.group()
def cli():
    pass


cli.add_command(c.compute)
cli.add_command(q.qc)


if __name__ == '__main__':
    cli()

