import click

from . import compute as c
from . import demultiplex as d
from . import quality_control as q
from . import simulate as s


@click.group()
def cli():
    pass


cli.add_command(c.compute)
cli.add_command(d.demultiplex)
cli.add_command(q.qc)
cli.add_command(s.simulate)


if __name__ == '__main__':
    cli()

