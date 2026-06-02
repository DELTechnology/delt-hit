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
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import Font
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
SUPP_DIR = ROOT / "paper" / "supplementary"
DEMUX_DIR = ROOT / "benchmarks" / "demultiplex"
DEMUX_TOOLS_DIR = DEMUX_DIR / "tools"
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


def template_reaction_rows() -> list[dict[str, str]]:
    return [
        {"name": "MTT deprotection", "smirks": "[N:1]C(c1ccccc1)(c2ccccc2)c3ccc(C)cc3>>[N:1]"},
        {"name": "Nvoc deprotection", "smirks": "[N:1]C(OCc1cc(OC)c(OC)cc1[N+]([O-])=O)=O>>[N:1]"},
        {"name": "Fmoc deprotection", "smirks": "[N:1]C(OCC1c2c(c3c1cccc3)cccc2)=O>>[N:1]"},
        {"name": "Staudinger reduction", "smirks": "[#6:1][$([NX2-][NX2+]#[NX1]),$([NX2]=[NX2+]=[NX1-])]>>[#6:1][N;H2]"},
        {"name": "Azido transfer", "smirks": "[#6:1][N;H2:2]>>[#6:1][N:2]=[N+]=[N-]"},
        {"name": "Nitro reduction", "smirks": "[#6:1][N+]([O-])=O>>[#6:1][N;H2]"},
        {
            "name": "Ester cleavage",
            "smirks": "[CX3:1](=[O:2])[OX2;H0][$([C;H3]),$([C;H2][C;H3]),$([C;H0]([C;H3])([C;H3])[C;H3])]>>[CX3:1](=[O:2])[OX2;H1]",
        },
        {"name": "Reductive amination", "smirks": "[C$(C(=O)([CX4,c])([CX4,c])),C$([CH](=O)([CX4,c])):1](=[O]).[N:2]>>[C:1][N:2]"},
        {"name": "Heck", "smirks": "[cX3:1][Br,I].[CX3:2]=[CX3;H2:3]>>[cX3:1]-[CX3:3]=[CX3:2]"},
        {"name": "Thioether formation", "smirks": "[#6:1][S;H1:2].[C:3][Cl,Br,I]>>[#6:1][S:2][C:3]"},
        {"name": "Amide bond formation", "smirks": "[CX3:1](=[O:2])[OX2;H1].[N;$([N;H1,H2]),$([N;H3]):4]>>[CX3:1](=[O:2])[N:4]"},
        {"name": "Sonogashira", "smirks": "[cX3:1][Br,I].[CX2:2]#[CX2;H1:3]>>[cX3:1]-[CX2:3]#[CX2:2]"},
        {"name": "Suzuki", "smirks": "[cX3:1][Br,I].[#6:2][BX3]>>[cX3:1][#6:2]"},
        {"name": "C-N coupling reactions", "smirks": "[cX3:1][Cl,Br,I].[N;H1,H2;!$(NC=O);!$(NS=O):2]>>[cX3:1][N:2]"},
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
    favalli_rows, puredel_rows = reanalysis_example_rows()

    workbook = Workbook()
    readme = workbook.active
    readme.title = "README"
    readme.append(["Supplementary Data 1"])
    readme.append(
        [
            "Source data for the Supplementary Figures, supplementary comparison and benchmark tables, and additional reaction templates."
        ]
    )
    readme.append(["Workflow data archive (Zenodo)", "https://doi.org/10.5281/zenodo.20447074"])
    readme.append(["Code release archive (Zenodo)", "https://doi.org/10.5281/zenodo.20497548"])
    readme.append(["Supporting-material repository", "https://doi.org/10.6084/m9.figshare.31198468"])
    readme.append(["Example single-display direct download", "https://ndownloader.figshare.com/files/61487743"])
    readme.append(["Download links used by supporting_material/experiments/*/download.sh", None])
    readme.append([
        "example-single-display/download.sh",
        "https://api.figshare.com/v2/articles/31198468 ; https://ndownloader.figshare.com/files/61487743",
    ])
    readme.append([
        "favalli/download.sh",
        "https://polybox.ethz.ch/index.php/s/9xZkcW5jWtLHiFo/download ; https://polybox.ethz.ch/index.php/s/cgxBWkHnpkT6rBi/download",
    ])
    readme.append([
        "pure-del/download.sh",
        "https://polybox.ethz.ch/index.php/s/H8YB6wdxpwWgzLB/download ; https://polybox.ethz.ch/index.php/s/wK6t8qjDTay8rAg/download ; https://zenodo.org/records/11122901/files/Keller_Petrov_etal_Supplementary_Material2_Sequencing_Data.zip?download=1",
    ])
    readme.append(["Generated by", "paper/supplementary/build_submission_assets.py"])
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
    add_sheet(workbook, "SuppTable1_Favalli", favalli_rows)
    add_sheet(workbook, "SuppTable2_PureDEL", puredel_rows)
    add_sheet(workbook, "SuppTable3_Workflow", workflow_comparison_rows())
    add_sheet(workbook, "SuppTable6_ReactionTemplates", template_reaction_rows())

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
