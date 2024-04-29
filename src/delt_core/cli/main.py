import click

from . import process as p


@click.group()
def cli():
    pass


cli.add_command(p.compute)
cli.add_command(p.evaluate)


if __name__ == '__main__':
    cli()
