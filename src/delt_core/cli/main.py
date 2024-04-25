import click

from . import evaluate


@click.group()
def cli():
    pass

cli.add_command(evaluate)


if __name__ == '__main__':
    cli()
