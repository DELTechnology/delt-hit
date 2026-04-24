#!/usr/bin/env python3
"""Run prepared DELi and DELT-Hit benchmarks for one dataset and split timings by phase."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import time
from pathlib import Path
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS_ROOT = PROJECT_ROOT / "benchmarks"
DATA_ROOT = BENCHMARKS_ROOT / "data"
TOOLS_ROOT = BENCHMARKS_ROOT / "tools"
DEFAULT_DATASET_NAME = "synthetic_2cycle_1m"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "target" / "benchmarks"
DELT_HIT_PYTHON = PROJECT_ROOT / ".venv" / "bin" / "python"
DELI_BIN = PROJECT_ROOT / "other_tools" / "DELi" / ".venv" / "bin" / "deli"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-name", required=True, help="Dataset name under benchmarks/data and benchmarks/tools.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory for logs and summary. Defaults to target/benchmarks/<dataset>/split_timing.",
    )
    return parser.parse_args()


def read_manifest(dataset_name: str) -> dict:
    return json.loads((DATA_ROOT / dataset_name / "manifest.json").read_text())


def require_path(path: Path, description: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def benchmark_environment(output_dir: Path, tool_env_bin: Path | None = None) -> dict[str, str]:
    env = os.environ.copy()
    mpl_dir = output_dir / ".mplconfig"
    xdg_dir = output_dir / ".xdg-cache"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    xdg_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    env["XDG_CACHE_HOME"] = str(xdg_dir)
    if tool_env_bin is not None:
        env["PATH"] = f"{tool_env_bin}:{env['PATH']}"
    return env


def run_command(cmd: list[str | Path], *, env: dict[str, str], cwd: Path, log_path: Path) -> float:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with log_path.open("w") as handle:
        handle.write("$ " + " ".join(str(item) for item in cmd) + "\n\n")
        handle.flush()
        subprocess.run([str(item) for item in cmd], cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    return time.perf_counter() - start


def read_expected_counts(path: Path) -> dict[tuple[int, ...], int]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    code_cols = [col for col in rows[0].keys() if col.startswith("code_")]
    return {tuple(int(row[col]) for col in code_cols): int(row["count"]) for row in rows}


def read_delt_counts(path: Path) -> dict[tuple[int, ...], int]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    code_cols = [col for col in rows[0].keys() if col.startswith("code_")]
    return {tuple(int(row[col]) + 1 for col in code_cols): int(row["count"]) for row in rows}


def read_deli_counts(path: Path) -> dict[tuple[int, ...], int]:
    rows = pd.read_parquet(path).to_dict("records")
    counts: dict[tuple[int, ...], int] = {}
    for row in rows:
        codes = tuple(int(part.rsplit("_", 1)[-1]) for part in row["bb_ids"].split(","))
        counts[codes] = int(row["count"])
    return counts


def compare_counts(expected: dict[tuple[int, ...], int], observed: dict[tuple[int, ...], int]) -> dict[str, object]:
    mismatches = []
    all_keys = sorted(set(expected) | set(observed))
    for key in all_keys:
        if expected.get(key, 0) != observed.get(key, 0):
            mismatches.append(
                {
                    "codes": list(key),
                    "expected": expected.get(key, 0),
                    "observed": observed.get(key, 0),
                }
            )
            if len(mismatches) >= 10:
                break
    return {
        "counts_match": not mismatches and len(expected) == len(observed),
        "expected_rows": len(expected),
        "observed_rows": len(observed),
        "expected_total_reads": sum(expected.values()),
        "observed_total_reads": sum(observed.values()),
        "mismatch_examples": mismatches,
    }


def clean_delt_outputs(config_path: Path) -> None:
    config = yaml.safe_load(config_path.read_text())
    save_dir = Path(config["experiment"]["save_dir"])
    experiment_name = config["experiment"]["name"]
    for rel_path in [
        Path(experiment_name) / "demultiplex" / "cutadapt_output_files",
        Path(experiment_name) / "selections",
    ]:
        target = save_dir / rel_path
        if target.exists():
            shutil.rmtree(target)


def run_deli(dataset_name: str, output_dir: Path, expected: dict[tuple[int, ...], int]) -> dict[str, object]:
    dataset_root = require_path(TOOLS_ROOT / "deli" / dataset_name, "DELi prepared dataset")
    deli_bin = require_path(DELI_BIN, "DELi CLI")
    selection_file = require_path(dataset_root / "decode_synthetic.yaml", "DELi selection file")
    decode_settings = require_path(dataset_root / "decode_settings_v02.yaml", "DELi decode settings")
    config_file = require_path(dataset_root / "deli_config", "DELi config file")

    run_root = output_dir / "deli"
    if run_root.exists():
        shutil.rmtree(run_root)
    run_root.mkdir(parents=True, exist_ok=True)
    env = benchmark_environment(run_root, DELI_BIN.parent)

    prefix = dataset_name
    timings = {
        "decode_run_s": run_command(
            [deli_bin, "--config-file", config_file, "decode", "run", selection_file, "--decode-settings-file", decode_settings, "--out-dir", run_root, "--prefix", prefix, "--skip-report"],
            env=env,
            cwd=PROJECT_ROOT,
            log_path=run_root / "logs" / "decode_run.log",
        ),
        "decode_collect_s": run_command(
            [deli_bin, "--config-file", config_file, "decode", "collect", run_root / f"{prefix}_decoded.tsv", "--out-loc", run_root / f"{prefix}_collected.ndjson"],
            env=env,
            cwd=PROJECT_ROOT,
            log_path=run_root / "logs" / "decode_collect.log",
        ),
        "decode_count_s": run_command(
            [deli_bin, "--config-file", config_file, "decode", "count", run_root / f"{prefix}_collected.ndjson", "--out-loc", run_root / f"{prefix}_counts.parquet", "--output-format", "parquet"],
            env=env,
            cwd=PROJECT_ROOT,
            log_path=run_root / "logs" / "decode_count.log",
        ),
    }
    timings["total_s"] = timings["decode_run_s"] + timings["decode_collect_s"] + timings["decode_count_s"]
    observed = read_deli_counts(run_root / f"{prefix}_counts.parquet")
    comparison = compare_counts(expected, observed)
    return {
        "tool": "deli",
        "output_dir": str(run_root),
        "timings": timings,
        **comparison,
    }


def run_delt(dataset_name: str, output_dir: Path, expected: dict[tuple[int, ...], int]) -> dict[str, object]:
    dataset_root = require_path(TOOLS_ROOT / "delt" / dataset_name, "DELT-Hit prepared dataset")
    tool_python = require_path(DELT_HIT_PYTHON, "DELT-Hit Python executable")
    config_path = require_path(dataset_root / "config.yaml", "DELT-Hit config")

    config = yaml.safe_load(config_path.read_text())
    experiment_name = config["experiment"]["name"]
    demultiplex_script = require_path(
        dataset_root / experiment_name / "demultiplex" / "cutadapt_input_files" / "demultiplex.sh",
        "prepared DELT-Hit demultiplex script",
    )
    counts_path = dataset_root / experiment_name / "selections" / "synthetic_selection" / "counts.txt"

    clean_delt_outputs(config_path)
    run_root = output_dir / "delt"
    run_root.mkdir(parents=True, exist_ok=True)
    env = benchmark_environment(run_root, DELT_HIT_PYTHON.parent)

    timings = {
        "demultiplex_s": run_command(
            ["bash", "-e", demultiplex_script],
            env=env,
            cwd=PROJECT_ROOT,
            log_path=run_root / "logs" / "demultiplex.log",
        ),
        "count_aggregation_s": run_command(
            [tool_python, "-m", "delt_hit.cli.main", "demultiplex", "process", "--config_path", config_path],
            env=env,
            cwd=PROJECT_ROOT,
            log_path=run_root / "logs" / "count_aggregation.log",
        ),
    }
    timings["total_s"] = timings["demultiplex_s"] + timings["count_aggregation_s"]
    observed = read_delt_counts(require_path(counts_path, "DELT-Hit counts file"))
    comparison = compare_counts(expected, observed)
    return {
        "tool": "delt",
        "output_dir": str(dataset_root / experiment_name),
        "timings": timings,
        **comparison,
    }


def print_summary(dataset_name: str, manifest: dict, deli: dict[str, object], delt: dict[str, object]) -> None:
    print(f"Dataset: {dataset_name}")
    print(f"Expected compounds: {manifest['expected_compounds']}")
    print(f"Expected reads: {manifest['expected_reads']}")
    print("")
    print("DELi")
    print(f"  decode run:      {deli['timings']['decode_run_s']:.3f}s")
    print(f"  decode collect:  {deli['timings']['decode_collect_s']:.3f}s")
    print(f"  decode count:    {deli['timings']['decode_count_s']:.3f}s")
    print(f"  total:           {deli['timings']['total_s']:.3f}s")
    print(f"  counts match:    {deli['counts_match']}")
    print("")
    print("DELT-Hit")
    print(f"  demultiplex:     {delt['timings']['demultiplex_s']:.3f}s")
    print(f"  aggregation:     {delt['timings']['count_aggregation_s']:.3f}s")
    print(f"  total:           {delt['timings']['total_s']:.3f}s")
    print(f"  counts match:    {delt['counts_match']}")


def main(*, dataset_name: str, output_dir: Path | None = None) -> None:
    dataset_dir = require_path(DATA_ROOT / dataset_name, "dataset directory")
    require_path(TOOLS_ROOT / "deli" / dataset_name, "DELi prepared tool directory")
    require_path(TOOLS_ROOT / "delt" / dataset_name, "DELT-Hit prepared tool directory")

    output_dir = (output_dir or (DEFAULT_OUTPUT_ROOT / dataset_name / "split_timing")).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = read_manifest(dataset_name)
    expected = read_expected_counts(require_path(dataset_dir / "expected_counts.tsv", "expected counts table"))
    deli_result = run_deli(dataset_name, output_dir, expected)
    delt_result = run_delt(dataset_name, output_dir, expected)

    summary = {
        "dataset_name": dataset_name,
        "dataset_dir": str(dataset_dir),
        "expected_compounds": manifest["expected_compounds"],
        "expected_reads": manifest["expected_reads"],
        "deli": deli_result,
        "delt_hit": delt_result,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print_summary(dataset_name, manifest, deli_result, delt_result)
    print(f"\nWrote summary to {summary_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(dataset_name=cli_args.dataset_name, output_dir=cli_args.output_dir)
