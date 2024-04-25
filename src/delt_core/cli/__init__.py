import click

from ..processor import compute_products


@click.command()
@click.argument('input_path', required=True, type=click.Path(writable=True))
@click.argument('output_path', required=False, type=click.Path(writable=True))
def evaluate(**kwargs):
    compute_products(**kwargs)

