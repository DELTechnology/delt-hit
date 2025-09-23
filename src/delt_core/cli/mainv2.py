from jsonargparse import CLI

from delt_core.cli.analyse.api import Analyse
from delt_core.cli.properties.api import Properties
from delt_core.cli.demultiplex.api import Demultiplex
from delt_core.cli.assembly.api import Assembly  # add/remove as your modules land

def cli() -> None:
    """
    Entry point for the delt-core CLI.

    This creates a two-level command tree:
      delt-core <group> <method> [--args]
    where <group> is one of the keys below (demultiplex, assembly),
    and <method> is any public method on that class (init, run, ...).
    """
    CLI(
        {
            "demultiplex": Demultiplex,
            "assembly": Assembly,
            "analyse": Analyse,
            "properties": Properties,
        },
        prog="delt-core",
        description="DEL-T core toolkit",
    )

if __name__ == "__main__":
    cli()