from pathlib import Path

from ... import quality_control as q


def report() -> None:
    input_dir = Path.cwd() / 'cutadapt_output_files'
    q.print_report(input_dir)


def plot(
        output_dir: Path,
) -> None:
    input_dir = Path.cwd() / 'cutadapt_output_files'
    if not output_dir:
        output_dir = Path.cwd() / 'plots'
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    q.plot_hits(input_dir, output_dir)

