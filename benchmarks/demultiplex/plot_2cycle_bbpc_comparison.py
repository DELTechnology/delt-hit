#!/usr/bin/env python3
"""Plot 2-cycle runtime comparisons across 10, 100, and 1000 BB/cycle."""

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
EXPECTED_BBPCS = (10, 100, 1000)


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
        help="Directory where the plot is written.",
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


def load_timings(tools_root: Path) -> dict[int, dict[str, list[tuple[int, float]]]]:
    grouped: dict[int, dict[str, list[tuple[int, float]]]] = defaultdict(lambda: defaultdict(list))
    for timing_path in sorted(tools_root.glob("*/**/timing.json")):
        report = json.loads(timing_path.read_text())
        parsed = parse_dataset_name(report["dataset_name"])
        if parsed is None:
            continue
        cycles, bbpc, expected_reads = parsed
        if cycles != 2 or bbpc not in EXPECTED_BBPCS:
            continue
        grouped[bbpc][report["tool"]].append((expected_reads, float(report["timings"]["total_s"])))
    return grouped


def format_reads(num_reads: int) -> str:
    if num_reads % 1_000_000 == 0:
        return f"{num_reads // 1_000_000}M"
    return str(num_reads)


def main(*, tools_root: Path, output_dir: Path) -> None:
    grouped = load_timings(tools_root.resolve())
    if not grouped:
        raise FileNotFoundError(
            f"No 2-cycle timing.json files found under {tools_root} for BB/cycle values {EXPECTED_BBPCS}."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(8, 5))

    tool_styles = {
        "deli": {"label": "DELi", "color": "#1f77b4", "marker": "o"},
        "delt": {"label": "DELT-Hit", "color": "#d62728", "marker": "s"},
    }
    bbpc_styles = {
        10: {"linestyle": "-", "suffix": "10 BB/cycle"},
        100: {"linestyle": "--", "suffix": "100 BB/cycle"},
        1000: {"linestyle": ":", "suffix": "1000 BB/cycle"},
    }

    for tool_name in ("deli", "delt"):
        for bbpc in EXPECTED_BBPCS:
            points = sorted(grouped[bbpc].get(tool_name, []), key=lambda item: item[0])
            if not points:
                continue
            x_values = [item[0] for item in points]
            y_values = [item[1] for item in points]
            tool_style = tool_styles[tool_name]
            bbpc_style = bbpc_styles[bbpc]
            axis.plot(
                x_values,
                y_values,
                label=f"{tool_style['label']} {bbpc_style['suffix']}",
                color=tool_style["color"],
                marker=tool_style["marker"],
                linestyle=bbpc_style["linestyle"],
                linewidth=2,
            )

    all_x_values = sorted(
        {
            expected_reads
            for bbpc_group in grouped.values()
            for tool_points in bbpc_group.values()
            for expected_reads, _runtime in tool_points
        }
    )
    if all_x_values:
        axis.set_xticks(all_x_values, [format_reads(value) for value in all_x_values])

    axis.set_xscale("log")
    axis.set_xlabel("Number of reads")
    axis.set_ylabel("Runtime (s)")
    axis.set_title("Synthetic 2-cycle benchmark runtimes across BB/cycle settings")
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    if axis.lines:
        axis.legend(ncol=2)

    figure.tight_layout()
    output_path = output_dir / "synthetic_2cycle_bbpc_comparison_runtime.png"
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(tools_root=cli_args.tools_root, output_dir=cli_args.output_dir)
