import click

@click.command()
def cli():
    click.echo('Welcome to DEL Technology!')

if __name__ == '__main__':
    cli()