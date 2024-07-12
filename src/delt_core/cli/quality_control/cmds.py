from pathlib import Path

from ... import quality_control as q


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
    from ...demultiplex.utils import Config
    from ...quality_control.compare_output import compare_counts_with_legacy
    config = Config.from_yaml(config_file).model_dump()
    compare_counts_with_legacy(config, legacy_results_dir)
