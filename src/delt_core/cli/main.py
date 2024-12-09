import click

from . import compute as c
from . import demultiplex as d
from . import init as i
from . import normalize as n
from . import quality_control as q
from . import simulate as s
from . import visualize as v

@click.group()
def cli():
    pass


cli.add_command(c.compute)
cli.add_command(d.demultiplex)
cli.add_command(i.init)
cli.add_command(n.normalize)
cli.add_command(q.qc)
cli.add_command(s.simulate)
cli.add_command(v.viz)


if __name__ == '__main__':
    cli()

