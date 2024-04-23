import click

from .compute_products import compute_products


@click.group()
def process():
    print('Processing...')

process.add_command(compute_products)
