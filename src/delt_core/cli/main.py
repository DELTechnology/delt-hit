import click

from ..processor import process


@click.group()
def cli():
    pass

@cli.command()
def welcome():
    click.echo('Welcome to DEL Technology!')

cli.add_command(process)

if __name__ == '__main__':
    cli()
