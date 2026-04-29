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


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TOOLS_ROOT = PROJECT_ROOT / "benchmarks" / "tools"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "benchmarks" / "plots"
DATASET_PATTERN = re.compile(r"^synthetic_(?P<cycles>\d+)cycle_(?P<depth>\d+)m$")


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


def parse_dataset_name(dataset_name: str) -> tuple[int, int] | None:
    match = DATASET_PATTERN.match(dataset_name)
    if match is None:
        return None
    return int(match.group("cycles")), int(match.group("depth")) * 1_000_000


def load_timings(tools_root: Path) -> dict[int, dict[str, list[tuple[int, float, str]]]]:
    grouped: dict[int, dict[str, list[tuple[int, float, str]]]] = defaultdict(lambda: defaultdict(list))
    for timing_path in sorted(tools_root.glob("*/**/timing.json")):
        report = json.loads(timing_path.read_text())
        parsed = parse_dataset_name(report["dataset_name"])
        if parsed is None:
            continue
        cycles, expected_reads = parsed
        grouped[cycles][report["tool"]].append((expected_reads, float(report["timings"]["total_s"]), report["dataset_name"]))
    return grouped


def format_reads(num_reads: int) -> str:
    if num_reads % 1_000_000 == 0:
        return f"{num_reads // 1_000_000}M"
    return str(num_reads)


def plot_cycle_group(*, cycles: int, tool_points: dict[str, list[tuple[int, float, str]]], output_dir: Path) -> Path:
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
    axis.set_title(f"Synthetic {cycles}-cycle benchmark runtimes")
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    axis.legend()

    all_x_values = sorted({item[0] for points in tool_points.values() for item in points})
    if all_x_values:
        axis.set_xticks(all_x_values, [format_reads(value) for value in all_x_values])

    figure.tight_layout()
    output_path = output_dir / f"synthetic_{cycles}cycle_runtime.png"
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


def main(*, tools_root: Path, output_dir: Path) -> None:
    grouped = load_timings(tools_root.resolve())
    if not grouped:
        raise FileNotFoundError(f"No timing.json files found under {tools_root}")

    output_paths = []
    for cycles, tool_points in sorted(grouped.items()):
        output_paths.append(plot_cycle_group(cycles=cycles, tool_points=tool_points, output_dir=output_dir.resolve()))

    for output_path in output_paths:
        print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(tools_root=cli_args.tools_root, output_dir=cli_args.output_dir)
