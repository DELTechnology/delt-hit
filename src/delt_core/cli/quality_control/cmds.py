from pathlib import Path

from ... import quality_control as q


def report(
        input_path: Path,
) -> None:
    q.print_report(input_path)


def plot(
        input_dir: Path,
        output_dir: Path,
) -> None:
    q.plot_hits(input_dir, output_dir)

