import click

from . import compute as c


@click.group()
def cli():
    pass


cli.add_command(c.compute)


if __name__ == '__main__':
    cli()
