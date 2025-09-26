from jsonargparse import CLI

from delt_core.cli.init import init
from delt_core.cli.demultiplex.api import Demultiplex
from delt_core.cli.analyse.api import Analyse
from delt_core.cli.library.api import Library
from delt_core.cli.dashboard.api import dashboard

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
            "init": init,
            "library": Library,
            "demultiplex": Demultiplex,
            "analyse": Analyse,
            "dashboard": dashboard,
        },
        prog="delt-core",
        description="DEL-T core toolkit",
    )

if __name__ == "__main__":
    cli()