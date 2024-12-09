import click
from pathlib import Path
from . import cmds


@click.group()
def viz():
    pass


@viz.command()
@click.option(
    '--path_to_file',
    required=True,
    type=click.Path(exists=True, path_type=Path),
)
def show(path_to_file: Path):
    from ...visualizations.app import main
    import pandas as pd
    df = pd.read_csv(path_to_file, sep='\t')
    main(data=df)