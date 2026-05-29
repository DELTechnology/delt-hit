#!/usr/bin/env python3
"""Build Nature Protocols supplementary submission assets.

The script keeps generated submission assets reproducible:

* RGB JPEG copies of all Supplementary Figures.
* A single Supplementary Data workbook with one tab per figure source dataset
  where practical, plus benchmark tables used in the Supplementary Notes.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from statistics import median
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import Font
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
SUPP_DIR = ROOT / "paper" / "supplementary"
DEMUX_DIR = ROOT / "benchmarks" / "demultiplex"
DEMUX_TOOLS_DIR = DEMUX_DIR / "tools"
NF2_DIR = ROOT / "supporting_material" / "experiments" / "favalli" / "benchmark_nf2"
WORKBOOK_PATH = SUPP_DIR / "Supplementary_Data_1.xlsx"

DATASET_PATTERN = re.compile(
    r"^synthetic_(?P<cycles>\d+)cycle_(?P<bbpc>\d+)bbpc_(?P<depth>\d+)m(?:_err=\d+)?$"
)


def dataset_metadata(dataset_name: str) -> dict[str, int] | None:
    match = DATASET_PATTERN.match(dataset_name)
    if match is None:
        return None
    return {
        "cycles": int(match.group("cycles")),
        "bb_per_cycle": int(match.group("bbpc")),
        "reads": int(match.group("depth")) * 1_000_000,
    }


def convert_supplementary_figures() -> list[Path]:
    figure_paths = [
        DEMUX_DIR / "plots" / "synthetic_2cycle_10bbpc_runtime.png",
        DEMUX_DIR / "plots" / "synthetic_3cycle_10bbpc_runtime.png",
        DEMUX_DIR / "plots" / "synthetic_4cycle_10bbpc_runtime.png",
        DEMUX_DIR / "plots" / "synthetic_2cycle_bbpc_comparison_runtime.png",
        DEMUX_DIR / "plots" / "synthetic_2cycle_10bbpc_peak_rss.png",
        DEMUX_DIR / "plots" / "synthetic_3cycle_10bbpc_peak_rss.png",
        DEMUX_DIR / "plots" / "synthetic_4cycle_10bbpc_peak_rss.png",
        SUPP_DIR / "figures" / "nf2-recall-benchmark.png",
    ]
    written: list[Path] = []
    for path in figure_paths:
        destination = path.with_suffix(".jpg")
        with Image.open(path) as image:
            rgb_image = image.convert("RGB")
            rgb_image.save(destination, format="JPEG", quality=95, dpi=(300, 300))
        written.append(destination)
    return written


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def demux_timing_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for timing_path in sorted(DEMUX_TOOLS_DIR.glob("*/**/timing.json")):
        report = json.loads(timing_path.read_text())
        metadata = dataset_metadata(report["dataset_name"])
        if metadata is None:
            continue
        row = {
            "tool": report["tool"],
            "dataset": report["dataset_name"],
            **metadata,
            "total_runtime_s": report["timings"]["total_s"],
            "demultiplex_s": report["timings"].get("demultiplex_s"),
            "count_aggregation_s": report["timings"].get("count_aggregation_s"),
            "counts_match": report.get("counts_match"),
            "observed_rows": report.get("observed_rows"),
            "observed_total_reads": report.get("observed_total_reads"),
        }
        rows.append(row)
    return rows


def demux_memory_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stats_path in sorted(DEMUX_TOOLS_DIR.glob("*/**/job-stats.json")):
        parts = stats_path.relative_to(DEMUX_TOOLS_DIR).parts
        tool, dataset = parts[0], parts[1]
        metadata = dataset_metadata(dataset)
        if metadata is None:
            continue
        stats = json.loads(stats_path.read_text())
        peak_rss = stats.get("peak_rss_bytes")
        rows.append(
            {
                "tool": tool,
                "dataset": dataset,
                **metadata,
                "job_state": stats.get("state"),
                "elapsed_s": stats.get("elapsed_s"),
                "ncpus": stats.get("ncpus"),
                "peak_rss_bytes": peak_rss,
                "peak_rss_gib": None if peak_rss is None else peak_rss / 1024**3,
                "exit_code": stats.get("exit_code"),
            }
        )
    return rows


def recovery_summary_rows(strict_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    methods = [
        ("DELT-Hit counts", "delt_hit_counts_rank"),
        ("DELT-Hit z-score", "delt_hit_zscore_rank"),
        ("DELi NSC", "deli_nsc_rank"),
        ("DELi MLE", "deli_mle_rank"),
        ("DELT-Hit edgeR", "delt_hit_edger_rank"),
        ("DELi norm-z", "deli_normz_rank"),
    ]
    summary: list[dict[str, Any]] = []
    for method, column in methods:
        ranks = [int(row[column]) for row in strict_rows if row[column]]
        n = len(ranks)
        summary.append(
            {
                "method": method,
                "n": n,
                "median_rank": median(ranks),
                "worst_rank": max(ranks),
                "recall_at_10": sum(rank <= 10 for rank in ranks) / n,
                "recall_at_50": sum(rank <= 50 for rank in ranks) / n,
                "recall_at_100": sum(rank <= 100 for rank in ranks) / n,
            }
        )
    return summary


def recall_figure_rows(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in summary_rows:
        family = "DELT-Hit" if row["method"].startswith("DELT-Hit") else "DELi"
        for cutoff, column in (
            (10, "recall_at_10"),
            (50, "recall_at_50"),
            (100, "recall_at_100"),
        ):
            rows.append(
                {
                    "tool_family": family,
                    "method": row["method"],
                    "rank_cutoff": cutoff,
                    "recall": row[column],
                    "n": row["n"],
                }
            )
    return rows


def strict_truth_rows(strict_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in strict_rows:
        key = (row["paper_target"], row["paper_ab"], row["evidence"])
        entry = grouped.setdefault(
            key,
            {
                "target": row["paper_target"],
                "ab_pair": row["paper_ab"],
                "evidence": row["evidence"],
                "conditions": [],
                "instances": 0,
            },
        )
        entry["conditions"].append(row["strand"])
        entry["instances"] += 1
    return [
        {
            **entry,
            "conditions": ", ".join(entry["conditions"]),
        }
        for entry in grouped.values()
    ]


def reanalysis_example_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    favalli = [
        {"code_0": 268, "code_1": 278, "legacy": 62, "delt": 62, "identical": True},
        {"code_0": 67, "code_1": 278, "legacy": 59, "delt": 59, "identical": True},
        {"code_0": 67, "code_1": 60, "legacy": 53, "delt": 53, "identical": True},
        {"code_0": 3, "code_1": 278, "legacy": 51, "delt": 51, "identical": True},
        {"code_0": 268, "code_1": 60, "legacy": 49, "delt": 49, "identical": True},
        {"code_0": 204, "code_1": 16, "legacy": 47, "delt": 47, "identical": True},
        {"code_0": 67, "code_1": 304, "legacy": 43, "delt": 43, "identical": True},
        {"code_0": 193, "code_1": 16, "legacy": 42, "delt": 42, "identical": True},
        {"code_0": 67, "code_1": 230, "legacy": 42, "delt": 42, "identical": True},
        {"code_0": 268, "code_1": 230, "legacy": 39, "delt": 39, "identical": True},
    ]
    puredel = [
        {"code_0": 9, "code_1": 1, "code_2": 11, "legacy": 82237, "delt": 82237, "identical": True},
        {"code_0": 9, "code_1": 4, "code_2": 11, "legacy": 75428, "delt": 75428, "identical": True},
        {"code_0": 9, "code_1": 5, "code_2": 11, "legacy": 69178, "delt": 69178, "identical": True},
        {"code_0": 9, "code_1": 7, "code_2": 11, "legacy": 66951, "delt": 66951, "identical": True},
        {"code_0": 8, "code_1": 10, "code_2": 11, "legacy": 63161, "delt": 63161, "identical": True},
        {"code_0": 9, "code_1": 10, "code_2": 11, "legacy": 63147, "delt": 63147, "identical": True},
        {"code_0": 9, "code_1": 8, "code_2": 11, "legacy": 58942, "delt": 58942, "identical": True},
        {"code_0": 8, "code_1": 4, "code_2": 11, "legacy": 55386, "delt": 55386, "identical": True},
        {"code_0": 8, "code_1": 7, "code_2": 11, "legacy": 53387, "delt": 53387, "identical": True},
        {"code_0": 9, "code_1": 2, "code_2": 11, "legacy": 49307, "delt": 49307, "identical": True},
    ]
    return favalli, puredel


def workflow_comparison_rows() -> list[dict[str, str]]:
    return [
        {
            "feature": "Input files",
            "delt_hit": "Excel workbook compiled to config.yaml",
            "deli": "Library JSON, building-block CSV files, and decode-settings YAML",
        },
        {
            "feature": "User-friendliness",
            "delt_hit": "Single workbook keeps library design, barcode definitions, and analyses in one file",
            "deli": "More modular and explicit, but usually requires editing several structured input files",
        },
        {"feature": "Raw read decoding", "delt_hit": "yes", "deli": "yes"},
        {"feature": "Error-tolerant barcode matching", "delt_hit": "yes", "deli": "yes"},
        {"feature": "Library enumeration", "delt_hit": "yes", "deli": "yes"},
        {
            "feature": "Enrichment / analysis methods",
            "delt_hit": "Counts-based enrichment, normalized z-score, and edgeR",
            "deli": "MLE, NSC/SD_min, normalized z-score, log-space z-score, PolyO, overlap analyses",
        },
        {
            "feature": "Analysis reports",
            "delt_hit": "Organized analysis folders and static outputs",
            "deli": "Automated HTML report with figures and CSV exports",
        },
        {"feature": "Dual-Display", "delt_hit": "yes", "deli": ""},
        {"feature": "Counts-driven filtered enumeration", "delt_hit": "yes", "deli": ""},
        {"feature": "Interactive count-table dashboard", "delt_hit": "yes", "deli": ""},
        {"feature": "Single-compound structure lookup", "delt_hit": "", "deli": "yes"},
        {"feature": "Molecular property and representation generation", "delt_hit": "yes", "deli": ""},
    ]


def add_sheet(workbook: Workbook, title: str, rows: list[dict[str, Any]]) -> None:
    sheet = workbook.create_sheet(title)
    if not rows:
        return
    headers = list(rows[0].keys())
    sheet.append(headers)
    for cell in sheet[1]:
        cell.font = Font(bold=True)
    for row in rows:
        sheet.append([row.get(header) for header in headers])


def build_workbook() -> Path:
    strict_rows = read_csv(NF2_DIR / "strict14_method_ranks.csv")
    summary_rows = recovery_summary_rows(strict_rows)
    favalli_rows, puredel_rows = reanalysis_example_rows()

    workbook = Workbook()
    readme = workbook.active
    readme.title = "README"
    readme.append(["Supplementary Data 1"])
    readme.append(
        [
            "Source data for Supplementary Figures and reviewer-relevant "
            "benchmark tables in the DELT-Hit Nature Protocols supplementary notes."
        ]
    )
    readme.append(["Generated by", "paper/supplementary/build_submission_assets.py"])
    readme.append(["NF2 rank source", str(NF2_DIR / "strict14_method_ranks.csv")])
    readme.append(["Demultiplex source", str(DEMUX_TOOLS_DIR)])
    for cell in readme["1:1"]:
        cell.font = Font(bold=True)

    timing_rows = demux_timing_rows()
    memory_rows = demux_memory_rows()
    add_sheet(
        workbook,
        "SuppFig1_runtime_cycles",
        [
            row
            for row in timing_rows
            if row["bb_per_cycle"] == 10 and row["cycles"] in (2, 3, 4)
        ],
    )
    add_sheet(
        workbook,
        "SuppFig2_runtime_bbpc",
        [
            row
            for row in timing_rows
            if row["cycles"] == 2 and row["bb_per_cycle"] in (10, 100, 1000)
        ],
    )
    add_sheet(
        workbook,
        "SuppFig3_peak_memory",
        [
            row
            for row in memory_rows
            if row["bb_per_cycle"] == 10 and row["cycles"] in (2, 3, 4)
        ],
    )
    add_sheet(workbook, "SuppFig4_NF2_recall", recall_figure_rows(summary_rows))
    add_sheet(workbook, "SuppTable1_Favalli", favalli_rows)
    add_sheet(workbook, "SuppTable2_PureDEL", puredel_rows)
    add_sheet(workbook, "SuppTable3_Workflow", workflow_comparison_rows())
    add_sheet(workbook, "SuppTable4_NF2_truth", strict_truth_rows(strict_rows))
    add_sheet(workbook, "SuppTable5_NF2_summary", summary_rows)
    add_sheet(workbook, "Strict14_method_ranks", strict_rows)
    add_sheet(workbook, "Zscore_benchmark", read_csv(NF2_DIR / "zscore_outputs" / "benchmark_summary.csv"))

    workbook.save(WORKBOOK_PATH)
    return WORKBOOK_PATH


def main() -> None:
    figure_paths = convert_supplementary_figures()
    workbook_path = build_workbook()
    for figure_path in figure_paths:
        print(f"Wrote {figure_path}")
    print(f"Wrote {workbook_path}")


if __name__ == "__main__":
    main()
