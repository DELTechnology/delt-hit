import click

from . import compute as c


@click.group()
def cli():
    pass


cli.add_command(c.compute)
cli.add_command(c.evaluate)


if __name__ == '__main__':
    cli()
