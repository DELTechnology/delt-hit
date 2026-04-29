#!/usr/bin/env python3
"""Plot benchmark runtimes for DELi and DELT-Hit across synthetic libraries."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TOOLS_ROOT = PROJECT_ROOT / "benchmarks" / "demultiplex" / "tools"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "benchmarks" / "demultiplex" / "plots"
DATASET_PATTERN = re.compile(
    r"^synthetic_(?P<cycles>\d+)cycle_(?P<bbpc>\d+)bbpc_(?P<depth>\d+)m(?:_err=\d+)?$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tools-root",
        type=Path,
        default=TOOLS_ROOT,
        help="Directory containing per-tool benchmark outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the plots are written.",
    )
    return parser.parse_args()


def parse_dataset_name(dataset_name: str) -> tuple[int, int, int] | None:
    match = DATASET_PATTERN.match(dataset_name)
    if match is None:
        return None
    return (
        int(match.group("cycles")),
        int(match.group("bbpc")),
        int(match.group("depth")) * 1_000_000,
    )


def load_timings(tools_root: Path) -> dict[tuple[int, int], dict[str, list[tuple[int, float, str]]]]:
    grouped: dict[tuple[int, int], dict[str, list[tuple[int, float, str]]]] = defaultdict(lambda: defaultdict(list))
    for timing_path in sorted(tools_root.glob("*/**/timing.json")):
        report = json.loads(timing_path.read_text())
        parsed = parse_dataset_name(report["dataset_name"])
        if parsed is None:
            continue
        cycles, bbpc, expected_reads = parsed
        grouped[(cycles, bbpc)][report["tool"]].append(
            (expected_reads, float(report["timings"]["total_s"]), report["dataset_name"])
        )
    return grouped


def format_reads(num_reads: int) -> str:
    if num_reads % 1_000_000 == 0:
        return f"{num_reads // 1_000_000}M"
    return str(num_reads)


def plot_cycle_group(
    *,
    cycles: int,
    bbpc: int,
    tool_points: dict[str, list[tuple[int, float, str]]],
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7, 4.5))

    styles = {
        "deli": {"label": "DELi", "color": "#1f77b4", "marker": "o"},
        "delt": {"label": "DELT-Hit", "color": "#d62728", "marker": "s"},
    }

    for tool_name in ("deli", "delt"):
        points = sorted(tool_points.get(tool_name, []), key=lambda item: item[0])
        if not points:
            continue
        x_values = [item[0] for item in points]
        y_values = [item[1] for item in points]
        style = styles[tool_name]
        axis.plot(
            x_values,
            y_values,
            label=style["label"],
            color=style["color"],
            marker=style["marker"],
            linewidth=2,
        )

    axis.set_xscale("log")
    axis.set_xlabel("Number of reads")
    axis.set_ylabel("Runtime (s)")
    axis.set_title(f"Synthetic {cycles}-cycle {bbpc} BB/cycle benchmark runtimes")
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    axis.legend()

    all_x_values = sorted({item[0] for points in tool_points.values() for item in points})
    if all_x_values:
        axis.set_xticks(all_x_values, [format_reads(value) for value in all_x_values])

    figure.tight_layout()
    output_path = output_dir / f"synthetic_{cycles}cycle_{bbpc}bbpc_runtime.png"
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


def ensure_expected_cycle_plots(
    grouped: dict[tuple[int, int], dict[str, list[tuple[int, float, str]]]]
) -> list[tuple[int, int]]:
    expected_groups = [(2, 10), (3, 10), (4, 10)]
    missing_groups = [group for group in expected_groups if group not in grouped]
    if missing_groups:
        missing_names = ", ".join(f"{cycles}-cycle/{bbpc}bbpc" for cycles, bbpc in missing_groups)
        raise FileNotFoundError(
            "Missing timing data for expected synthetic benchmark groups: "
            f"{missing_names}. Expected one plot each for the current 2/3/4-cycle libraries."
        )
    return expected_groups


def main(*, tools_root: Path, output_dir: Path) -> None:
    grouped = load_timings(tools_root.resolve())
    if not grouped:
        raise FileNotFoundError(f"No timing.json files found under {tools_root}")

    output_paths = []
    for cycles, bbpc in ensure_expected_cycle_plots(grouped):
        output_paths.append(
            plot_cycle_group(
                cycles=cycles,
                bbpc=bbpc,
                tool_points=grouped[(cycles, bbpc)],
                output_dir=output_dir.resolve(),
            )
        )

    for output_path in output_paths:
        print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(tools_root=cli_args.tools_root, output_dir=cli_args.output_dir)
