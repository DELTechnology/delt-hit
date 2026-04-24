#!/usr/bin/env python3
"""Run synthetic DELi and DELT-Hit benchmark workflows and compare outputs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DELT_HIT_PYTHON = ROOT / ".venv" / "bin" / "python"
DELI_BIN = ROOT / "other_tools" / "DELi" / ".venv" / "bin" / "deli"
DEFAULT_OUTPUT_DIR = ROOT / "target" / "benchmarks" / "synthetic"
DEFAULT_TOOLS = ("delt-hit", "deli")
DEFAULT_DELI_CYCLES = (2, 3, 4)
BB_COUNTS_BY_CYCLE = {1: 100, 2: 10, 3: 10, 4: 10}
EXPERIMENT_NAME = "synthetic_2cycle"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where generated fixtures, logs, and benchmark summaries are written.",
    )
    parser.add_argument(
        "--num-reads-per-compound",
        type=int,
        default=1,
        help="Number of replicated reads to emit per synthetic compound combination.",
    )
    parser.add_argument(
        "--tools",
        nargs="+",
        choices=("delt-hit", "deli"),
        default=list(DEFAULT_TOOLS),
        help="Benchmark DELT-Hit, DELi, or both.",
    )
    parser.add_argument(
        "--deli-cycles",
        nargs="+",
        type=int,
        default=list(DEFAULT_DELI_CYCLES),
        help="DELi cycle counts to benchmark. DELT-Hit currently only supports the synthetic 2-cycle fixture.",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Reuse an existing output directory instead of deleting it first.",
    )
    args = parser.parse_args()
    if args.num_reads_per_compound < 1:
        parser.error("num_reads_per_compound must be at least 1")
    unsupported = sorted(set(args.deli_cycles) - set(DEFAULT_DELI_CYCLES))
    if unsupported:
        parser.error(f"unsupported DELi cycle counts: {unsupported}; choose from {list(DEFAULT_DELI_CYCLES)}")
    return args


def benchmark_environment(run_root: Path) -> dict[str, str]:
    env = os.environ.copy()
    mpl_dir = run_root / ".mplconfig"
    xdg_dir = run_root / ".xdg-cache"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    xdg_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    env["XDG_CACHE_HOME"] = str(xdg_dir)
    env["PATH"] = f"{ROOT / '.venv' / 'bin'}:{ROOT / 'other_tools' / 'DELi' / '.venv' / 'bin'}:{env['PATH']}"
    return env


def run_command(
    cmd: list[str],
    *,
    env: dict[str, str],
    cwd: Path,
    log_path: Path,
) -> float:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n\n")
        handle.flush()
        subprocess.run(
            cmd,
            cwd=cwd,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=True,
        )
    return time.perf_counter() - start


def expected_compounds(cycle_count: int) -> int:
    total = 1
    for cycle in range(1, cycle_count + 1):
        total *= BB_COUNTS_BY_CYCLE[cycle]
    return total


def generate_deli_inputs(*, output_dir: Path, num_reads_per_compound: int, env: dict[str, str]) -> None:
    run_command(
        [
            str(DELT_HIT_PYTHON),
            str(ROOT / "scripts" / "generate_synthetic_deli_fixture.py"),
            "--output-root",
            str(output_dir),
            "--num-reads-per-compound",
            str(num_reads_per_compound),
        ],
        cwd=ROOT,
        env=env,
        log_path=output_dir / "logs" / "generate_deli_inputs.log",
    )


def generate_delt_hit_inputs(*, output_dir: Path, num_reads_per_compound: int, env: dict[str, str]) -> Path:
    run_command(
        [
            str(DELT_HIT_PYTHON),
            str(ROOT / "experiments" / "build_synthetic_2cycle_delt_hit_inputs.py"),
            "--output-dir",
            str(output_dir),
            "--num-reads-per-compound",
            str(num_reads_per_compound),
        ],
        cwd=ROOT,
        env=env,
        log_path=output_dir / "logs" / "generate_delt_hit_inputs.log",
    )
    return output_dir / f"{EXPERIMENT_NAME}.xlsx"


def summarize_deli_counts(path: Path) -> dict[str, object]:
    df = pd.read_parquet(path)
    total_molecules = int(df["count"].sum())
    return {
        "counted_rows": int(len(df)),
        "total_molecules": total_molecules,
        "raw_total_reads": int(df["raw_count"].sum()) if "raw_count" in df.columns else total_molecules,
        "libraries": sorted(df["library_id"].unique().tolist()),
    }


def run_deli_benchmark(
    *,
    cycle_count: int,
    fixture_root: Path,
    run_root: Path,
    env: dict[str, str],
    num_reads_per_compound: int,
) -> dict[str, object]:
    experiment_name = f"synthetic_{cycle_count}cycle"
    experiment_dir = fixture_root / experiment_name
    selection_file = experiment_dir / "decode_synthetic.yaml"
    config_file = experiment_dir / "deli_config"
    decode_settings = experiment_dir / "decode_settings_v02.yaml"
    out_dir = run_root / experiment_name
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = experiment_name
    timings = {
        "decode_run_s": run_command(
            [
                str(DELI_BIN),
                "--config-file",
                str(config_file),
                "decode",
                "run",
                str(selection_file),
                "--decode-settings-file",
                str(decode_settings),
                "--out-dir",
                str(out_dir),
                "--prefix",
                prefix,
                "--skip-report",
            ],
            cwd=ROOT,
            env=env,
            log_path=out_dir / "logs" / "decode_run.log",
        ),
        "decode_collect_s": run_command(
            [
                str(DELI_BIN),
                "--config-file",
                str(config_file),
                "decode",
                "collect",
                str(out_dir / f"{prefix}_decoded.tsv"),
                "--out-loc",
                str(out_dir / f"{prefix}_collected.ndjson"),
            ],
            cwd=ROOT,
            env=env,
            log_path=out_dir / "logs" / "decode_collect.log",
        ),
        "decode_count_s": run_command(
            [
                str(DELI_BIN),
                "--config-file",
                str(config_file),
                "decode",
                "count",
                str(out_dir / f"{prefix}_collected.ndjson"),
                "--out-loc",
                str(out_dir / f"{prefix}_counts.parquet"),
                "--output-format",
                "parquet",
            ],
            cwd=ROOT,
            env=env,
            log_path=out_dir / "logs" / "decode_count.log",
        ),
    }
    timings["total_s"] = sum(timings.values())

    decode_stats = json.loads((out_dir / f"{prefix}_decode_statistics.json").read_text())
    counts_summary = summarize_deli_counts(out_dir / f"{prefix}_counts.parquet")
    return {
        "tool": "deli",
        "cycle_count": cycle_count,
        "selection_file": str(selection_file),
        "output_dir": str(out_dir),
        "timings": timings,
        "expected_compounds": expected_compounds(cycle_count),
        "expected_reads": expected_compounds(cycle_count) * num_reads_per_compound,
        "decode_statistics": decode_stats,
        "counts_summary": counts_summary,
    }


def summarize_delt_hit_counts(path: Path) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    return {
        "counted_rows": int(len(df)),
        "total_reads": int(df["count"].sum()),
        "min_code_1": int(df["code_1"].min()),
        "max_code_1": int(df["code_1"].max()),
        "min_code_2": int(df["code_2"].min()),
        "max_code_2": int(df["code_2"].max()),
    }


def run_delt_hit_benchmark(
    *,
    input_root: Path,
    run_root: Path,
    env: dict[str, str],
) -> dict[str, object]:
    workbook = input_root / f"{EXPERIMENT_NAME}.xlsx"
    prefix = EXPERIMENT_NAME
    config_dir = input_root / prefix
    config_path = config_dir / "config.yaml"

    timings = {
        "init_s": run_command(
            [str(DELT_HIT_PYTHON), "-m", "delt_hit.cli.main", "init", "--excel_path", str(workbook)],
            cwd=ROOT,
            env=env,
            log_path=run_root / "logs" / "delt_hit_init.log",
        ),
        "demultiplex_prepare_s": run_command(
            [str(DELT_HIT_PYTHON), "-m", "delt_hit.cli.main", "demultiplex", "prepare", "--config_path", str(config_path)],
            cwd=ROOT,
            env=env,
            log_path=run_root / "logs" / "delt_hit_prepare.log",
        ),
    }

    demultiplex_script = config_dir / "demultiplex" / "cutadapt_input_files" / "demultiplex.sh"
    timings["demultiplex_shell_s"] = run_command(
        ["bash", "-e", str(demultiplex_script)],
        cwd=ROOT,
        env=env,
        log_path=run_root / "logs" / "delt_hit_demultiplex.log",
    )
    timings["demultiplex_process_s"] = run_command(
        [str(DELT_HIT_PYTHON), "-m", "delt_hit.cli.main", "demultiplex", "process", "--config_path", str(config_path)],
        cwd=ROOT,
        env=env,
        log_path=run_root / "logs" / "delt_hit_process.log",
    )
    timings["total_s"] = sum(timings.values())

    counts_path = config_dir / "selections" / "synthetic_selection" / "counts.txt"
    return {
        "tool": "delt-hit",
        "cycle_count": 2,
        "config_path": str(config_path),
        "output_dir": str(config_dir),
        "timings": timings,
        "expected_compounds": expected_compounds(2),
        "expected_reads": summarize_delt_hit_counts(counts_path)["total_reads"],
        "counts_summary": summarize_delt_hit_counts(counts_path),
        "counts_path": str(counts_path),
    }


def normalized_deli_counts(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    bb_codes = df["bb_ids"].str.split(",", expand=True)
    return pd.DataFrame(
        {
            "code_1": bb_codes[0].str.extract(r"(\d+)$").astype(int)[0],
            "code_2": bb_codes[1].str.extract(r"(\d+)$").astype(int)[0],
            "count": df["count"].astype(int),
        }
    ).sort_values(["code_1", "code_2"]).reset_index(drop=True)


def normalized_delt_hit_counts(path: Path) -> pd.DataFrame:
    return (
        pd.read_csv(path, sep="\t")[["code_1", "code_2", "count"]]
        .astype(int)
        .sort_values(["code_1", "code_2"])
        .reset_index(drop=True)
    )


def normalized_expected_counts(path: Path) -> pd.DataFrame:
    return (
        pd.read_csv(path, sep="\t")[["code_1", "code_2", "count"]]
        .astype(int)
        .sort_values(["code_1", "code_2"])
        .reset_index(drop=True)
    )


def compare_normalized_counts(
    *,
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    left_name: str,
    right_name: str,
) -> dict[str, object]:
    equal = left_df.equals(right_df)
    result = {
        "matched": bool(equal),
        f"row_count_{left_name}": int(len(left_df)),
        f"row_count_{right_name}": int(len(right_df)),
    }
    if not equal:
        merged = left_df.merge(
            right_df,
            on=["code_1", "code_2"],
            how="outer",
            suffixes=(f"_{left_name}", f"_{right_name}"),
        ).fillna(0)
        mismatches = merged[merged[f"count_{left_name}"] != merged[f"count_{right_name}"]]
        result["mismatch_examples"] = mismatches.head(10).to_dict("records")
        result["mismatch_count"] = int(len(mismatches))
    return result


def compare_two_cycle_outputs(*, deli_counts_path: Path, delt_hit_counts_path: Path) -> dict[str, object]:
    return compare_normalized_counts(
        left_df=normalized_deli_counts(deli_counts_path),
        right_df=normalized_delt_hit_counts(delt_hit_counts_path),
        left_name="deli",
        right_name="delt_hit",
    )


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and not args.keep_existing:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    env = benchmark_environment(output_dir)

    results: dict[str, object] = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "num_reads_per_compound": args.num_reads_per_compound,
        "tools": {},
    }

    if "deli" in args.tools:
        fixture_root = output_dir / "deli_fixtures"
        run_root = output_dir / "deli_runs"
        generate_deli_inputs(output_dir=fixture_root, num_reads_per_compound=args.num_reads_per_compound, env=env)
        deli_results = {}
        for cycle_count in args.deli_cycles:
            deli_results[f"{cycle_count}_cycle"] = run_deli_benchmark(
                cycle_count=cycle_count,
                fixture_root=fixture_root,
                run_root=run_root,
                env=env,
                num_reads_per_compound=args.num_reads_per_compound,
            )
        results["tools"]["deli"] = deli_results

    if "delt-hit" in args.tools:
        input_root = output_dir / "delt_hit_inputs"
        run_root = output_dir / "delt_hit_run"
        generate_delt_hit_inputs(output_dir=input_root, num_reads_per_compound=args.num_reads_per_compound, env=env)
        results["tools"]["delt-hit"] = run_delt_hit_benchmark(
            input_root=input_root,
            run_root=run_root,
            env=env,
        )

    if "deli" in results["tools"] and "2_cycle" in results["tools"]["deli"] and "delt-hit" in results["tools"]:
        results["two_cycle_comparison"] = compare_two_cycle_outputs(
            deli_counts_path=output_dir / "deli_runs" / "synthetic_2cycle" / "synthetic_2cycle_counts.parquet",
            delt_hit_counts_path=output_dir / "delt_hit_inputs" / EXPERIMENT_NAME / "selections" / "synthetic_selection" / "counts.txt",
        )

    summary_path = output_dir / "benchmark_summary.json"
    summary_path.write_text(json.dumps(results, indent=2) + "\n")
    print(f"Wrote benchmark summary to {summary_path}")
    if "two_cycle_comparison" in results:
        comparison = results["two_cycle_comparison"]
        status = "matched" if comparison["matched"] else "mismatched"
        print(f"Two-cycle DELi vs DELT-Hit counts {status}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
