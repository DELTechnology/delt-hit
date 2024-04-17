import click

@click.group()
def cli():
    pass

@cli.command()
def welcome():
    click.echo('Welcome to DEL Technology!')

# import command from other module
from ..processor import process
cli.add_command(process)

if __name__ == '__main__':
    cli()