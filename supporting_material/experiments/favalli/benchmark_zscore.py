#!/usr/bin/env python3
"""Benchmark the delt-hit normalized z-score method on the Favalli NF2 dataset.

Runs delt-hit z-score enrichment on the NF2 benchmark experiments and looks up
the rank of known positive-control compounds from Favalli et al. 2021.

Quick start (from the delt-hit repository root):

    pixi run python supporting_material/experiments/favalli/benchmark_zscore.py

This uses the self-contained benchmark_nf2/ subdirectory (30 pre-computed count
tables, 8 experiments, 14 strict-affinity positive controls).

Reference: Faver et al., ACS Combinatorial Science 2019, 21(2), 75-82.
           DOI: 10.1021/acscombsci.8b00116
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

# Allow running from delt-hit repo root without pip install
_HERE = Path(__file__).resolve().parent
_SRC = _HERE.parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from delt_hit.cli.analyse.api import Analyse

# ---------------------------------------------------------------------------
# Mapping: experiment name in analysis.yaml → paper_target in truth_set.csv
# Covers both ss (single-strand) and ds (double-strand) conditions.
# ---------------------------------------------------------------------------
_EXPERIMENT_TO_TARGET = {
    "ss_caix":       "CAIX",
    "ds_caix":       "CAIX",
    "ca9_ss":        "CAIX",
    "ca9_ds":        "CAIX",
    "ss_crebbp":     "CREBBP",
    "ds_crebbp":     "CREBBP",
    "ss_wt_pik3":    "wt-PI3K",
    "ds_wt_pik3":    "wt-PI3K",
    "ss_h1047r_pik3":"H1047R-PI3K",
    "ds_h1047r_pik3":"H1047R-PI3K",
}


def parse_args() -> argparse.Namespace:
    default_config = _HERE / "benchmark_nf2" / "analysis.yaml"
    default_truth  = _HERE / "benchmark_nf2" / "truth_set.csv"

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config,
        help=f"Path to delt-hit analysis YAML (default: {default_config}).",
    )
    parser.add_argument(
        "--truth-csv",
        type=Path,
        default=default_truth,
        help=f"CSV with positive controls from truth_set.csv (default: {default_truth}). "
             "Expected columns: paper_target, code_1, code_2, label.",
    )
    parser.add_argument(
        "--label-filter",
        default="strict_affinity",
        help="Which 'label' rows to use as positive controls (default: strict_affinity). "
             "Use 'all' to include all rows.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_HERE / "benchmark_nf2" / "zscore_outputs",
        help="Directory to write z-score analysis outputs.",
    )
    return parser.parse_args()


def load_experiment_names(config_path: Path) -> list[str]:
    cfg = yaml.safe_load(config_path.read_text())
    return [exp["name"] for exp in cfg["experiments"]]


def find_zscore_stats(config_path: Path, experiment_name: str) -> Path:
    """Locate zscore/stats.csv written by Analyse.enrichment()."""
    cfg = yaml.safe_load(config_path.read_text())
    exp = next(e for e in cfg["experiments"] if e["name"] == experiment_name)
    save_dir = Path(exp["save_dir"]).expanduser()
    if not save_dir.is_absolute():
        save_dir = config_path.parent / save_dir
    return save_dir.resolve() / experiment_name / "zscore" / "stats.csv"


def run_zscore_for_all(config_path: Path) -> None:
    """Run z-score enrichment for every experiment in the config."""
    analyse = Analyse()
    experiment_names = load_experiment_names(config_path)
    print(f"Running z-score for {len(experiment_names)} experiment(s):")
    orig_dir = Path.cwd()
    os.chdir(config_path.parent)
    try:
        for name in experiment_names:
            print(f"  {name} ...", end=" ", flush=True)
            analyse.enrichment(config_path=config_path, name=name, method="zscore")
            print("done")
    finally:
        os.chdir(orig_dir)


def lookup_rank(stats_path: Path, code_0: int, code_1: int) -> tuple[int | None, float | None]:
    """Return (1-based rank by descending protein zscore, protein zscore) for a compound."""
    if not stats_path.exists():
        return None, None
    stats = (
        pd.read_csv(stats_path)
        .sort_values("protein", ascending=False, na_position="last")
        .reset_index(drop=True)
    )
    stats["_rank"] = range(1, len(stats) + 1)
    row = stats[(stats["code_0"] == code_0) & (stats["code_1"] == code_1)]
    if row.empty:
        return None, None
    rank = int(row["_rank"].iloc[0])
    zscore = float(row["protein"].iloc[0]) if "protein" in row.columns else float("nan")
    return rank, zscore


def build_truth_entries(
    truth_csv: Path,
    available_experiments: set[str],
    label_filter: str,
) -> list[dict]:
    """Build truth entries from truth_set.csv, mapping targets to experiments.

    truth_set.csv uses 1-based code_1/code_2 (paper numbering); delt-hit
    stats.csv uses 0-based code_0/code_1.  We subtract 1 from each code.
    """
    truth_df = pd.read_csv(truth_csv)
    if label_filter != "all":
        truth_df = truth_df[truth_df["label"] == label_filter]

    entries = []
    for _, row in truth_df.iterrows():
        target = row["paper_target"]
        code_0 = int(row["code_1"])  # already 0-based (matching counts.txt code_0)
        code_1 = int(row["code_2"])  # already 0-based (matching counts.txt code_1)
        for exp_name, exp_target in _EXPERIMENT_TO_TARGET.items():
            if exp_target == target and exp_name in available_experiments:
                entries.append({
                    "experiment": exp_name,
                    "paper_target": target,
                    "paper_ab": row["paper_ab"],
                    "code_0": code_0,
                    "code_1": code_1,
                    "label": row["label"],
                    "evidence": row.get("evidence", ""),
                })
    return entries


def build_summary(config_path: Path, truth_entries: list[dict]) -> pd.DataFrame:
    rows = []
    for entry in truth_entries:
        stats_path = find_zscore_stats(config_path, entry["experiment"])
        rank, zscore = lookup_rank(stats_path, entry["code_0"], entry["code_1"])
        rows.append({
            "experiment":   entry["experiment"],
            "paper_target": entry["paper_target"],
            "paper_ab":     entry["paper_ab"],
            "label":        entry["label"],
            "code_0":       entry["code_0"],
            "code_1":       entry["code_1"],
            "zscore_rank":  rank,
            "zscore_value": round(zscore, 3) if zscore is not None else None,
        })
    return pd.DataFrame(rows)


def print_summary(summary: pd.DataFrame) -> None:
    print("\n=== Z-score positive control benchmark ===")
    print(summary.to_string(index=False))
    ranked = summary.dropna(subset=["zscore_rank"])
    if len(ranked) > 0:
        median_rank = ranked["zscore_rank"].median()
        recall_10  = (ranked["zscore_rank"] <= 10).mean()
        recall_50  = (ranked["zscore_rank"] <= 50).mean()
        recall_100 = (ranked["zscore_rank"] <= 100).mean()
        print(f"\nn (ranked) = {len(ranked)}")
        print(f"Median rank : {median_rank:.1f}")
        print(f"Recall@10   : {recall_10:.1%}")
        print(f"Recall@50   : {recall_50:.1%}")
        print(f"Recall@100  : {recall_100:.1%}")


def main() -> None:
    args = parse_args()

    config_path = args.config.resolve()
    if not config_path.exists():
        print(f"ERROR: config not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    if not args.truth_csv.exists():
        print(f"ERROR: truth CSV not found: {args.truth_csv}", file=sys.stderr)
        sys.exit(1)

    args.output.mkdir(parents=True, exist_ok=True)

    available_experiments = set(load_experiment_names(config_path))
    truth_entries = build_truth_entries(
        args.truth_csv, available_experiments, args.label_filter
    )
    if not truth_entries:
        print(
            f"ERROR: no truth entries found for the available experiments.\n"
            f"Available: {sorted(available_experiments)}\n"
            f"label_filter: {args.label_filter}",
            file=sys.stderr,
        )
        sys.exit(1)

    print(
        f"Loaded {len(truth_entries)} positive control instances "
        f"(label={args.label_filter!r}) from {args.truth_csv.name}"
    )

    run_zscore_for_all(config_path)

    summary = build_summary(config_path, truth_entries)
    print_summary(summary)

    summary_path = args.output / "benchmark_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\nSummary written to {summary_path}")


if __name__ == "__main__":
    main()
