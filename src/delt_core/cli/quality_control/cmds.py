from pathlib import Path

from delt_core import quality_control as q
from delt_core.cli.demultiplex.cmds import create_lists
from delt_core.demultiplex.validation import Config


def report(experiment_dir: Path) -> None:
    q.print_report(experiment_dir)


def plot(
        experiment_dir: Path,
        output_dir: Path = None,
) -> None:
    output_dir = output_dir or experiment_dir / 'quality_control' / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    q.plot_hits(experiment_dir, output_dir)


def compare_with_legacy(config_file: Path, legacy_results_dir: Path):
    config = Config.from_yaml(config_file).model_dump()
    q.compare_counts_with_legacy(config, legacy_results_dir)


def analyze_codons(
        config_file: Path,
) -> None:
    output_dir = Path(config_file).parent / 'edit_distances'
    output_dir.mkdir(parents=True, exist_ok=True)
    structure = create_lists(config_file)
    q.analyze_codons(structure, output_dir)

