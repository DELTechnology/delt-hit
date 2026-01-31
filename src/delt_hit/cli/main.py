from jsonargparse import CLI

from delt_hit.cli.init import init
from delt_hit.cli.demultiplex.api import Demultiplex
from delt_hit.cli.analyse.api import Analyse
from delt_hit.cli.library.api import Library
from delt_hit.cli.dashboard.api import dashboard

def cli() -> None:
    """
    Entry point for the delt-hit CLI.

    This creates a two-level command tree:
      delt-hit <group> <method> [--args]
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
        prog="delt-hit",
        description="DELT Hit toolkit",
    )

if __name__ == "__main__":
    cli()