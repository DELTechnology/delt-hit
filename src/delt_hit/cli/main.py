from jsonargparse import CLI

from delt_hit.cli.init import init
from delt_hit.cli.demultiplex.api import Demultiplex
from delt_hit.cli.analyse.api import Analyse
from delt_hit.cli.library.api import Library
from delt_hit.cli.dashboard.api import dashboard

def cli() -> None:
    """Run the delt-hit CLI entrypoint.

    The command tree is: ``delt-hit <group> <method> [--args]``.
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
        add_config_file_arg=False,
    )

if __name__ == "__main__":
    cli()
