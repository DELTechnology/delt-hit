import click

from . import compute as c
from . import demultiplex as d
from . import init as i
from . import quality_control as qc
from . import simulate as s


@click.group()
def cli():
    pass


cli.add_command(c.compute)
cli.add_command(d.demultiplex)
cli.add_command(i.init)
cli.add_command(qc.qc)
cli.add_command(s.simulate)


if __name__ == '__main__':
    cli()

