#!/usr/bin/env python3
"""Run counts and edgeR analyses on a synthetic enrichment benchmark dataset."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BENCHMARKS_ROOT = Path(__file__).resolve().parent
SRC_ROOT = PROJECT_ROOT / "src"
if str(BENCHMARKS_ROOT) not in sys.path:
    sys.path.insert(0, str(BENCHMARKS_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from delt_hit.cli.analyse.api import Analyse
from generate_synthetic_enrichment import DEFAULT_ANALYSIS_NAME, read_truth_table


DEFAULT_TOP_K = 100


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        required=True,
        help="Synthetic enrichment dataset directory produced by generate_synthetic_enrichment.py.",
    )
    parser.add_argument("--top-k", type=int, default=DEFAULT_TOP_K)
    parser.add_argument("--rscript-bin", default="Rscript")
    return parser.parse_args()


def read_analysis_config(config_path: Path) -> dict:
    return yaml.safe_load(config_path.read_text())


def find_analysis_root(config_path: Path, analysis_name: str) -> Path:
    config = read_analysis_config(config_path)
    experiment, = [item for item in config["experiments"] if item["name"] == analysis_name]
    return Path(experiment["save_dir"]).expanduser().resolve() / analysis_name


def run_r_script(*, rscript_bin: str, script_path: Path, log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        handle.write(f"$ {rscript_bin} --vanilla {script_path}\n\n")
        handle.flush()
        subprocess.run(
            [rscript_bin, "--vanilla", str(script_path)],
            cwd=PROJECT_ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run_analysis_method(
    *,
    analyse: Analyse,
    config_path: Path,
    analysis_name: str,
    method: str,
    analysis_root: Path,
    rscript_bin: str,
) -> Path:
    output_dir = analysis_root / method
    if output_dir.exists():
        shutil.rmtree(output_dir)

    analyse.enrichment(
        analysis_config=config_path,
        name=analysis_name,
        method=method,
        save_dir=analysis_root,
    )


    script_name = "enrichment_counts.R" if method == "counts" else "enrichment_edgeR.R"
    script_path = output_dir / script_name
    run_r_script(
        rscript_bin=rscript_bin,
        script_path=script_path,
        log_path=analysis_root / "logs" / f"{method}.log",
    )
    return output_dir


def load_counts_results(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    return table.rename(columns={"enrichment": "counts_score"})[["code_0", "code_1", "counts_score"]]


def load_edger_results(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    columns = ["code_0", "code_1", "logFC", "logCPM", "LR", "PValue", "FDR"]
    available = [column for column in columns if column in table.columns]
    renamed = table[available].rename(
        columns={
            "logFC": "edgeR_logFC",
            "logCPM": "edgeR_logCPM",
            "LR": "edgeR_LR",
            "PValue": "edgeR_PValue",
            "FDR": "edgeR_FDR",
        }
    )
    return renamed


def add_rank_column(table: pd.DataFrame, *, score_column: str, rank_column: str) -> pd.DataFrame:
    ranked = table.sort_values(by=score_column, ascending=False, na_position="last").reset_index(drop=True).copy()
    ranked[rank_column] = range(1, len(ranked) + 1)
    return ranked


def build_comparison_table(
    truth: pd.DataFrame,
    counts_results: pd.DataFrame,
    edger_results: pd.DataFrame,
) -> pd.DataFrame:
    comparison = truth.merge(counts_results, on=["code_0", "code_1"], how="left")
    comparison = comparison.merge(edger_results, on=["code_0", "code_1"], how="left")

    counts_ranked = add_rank_column(
        comparison[["id", "counts_score"]],
        score_column="counts_score",
        rank_column="counts_rank",
    )[["id", "counts_rank"]]
    edger_ranked = add_rank_column(
        comparison[["id", "edgeR_logFC"]],
        score_column="edgeR_logFC",
        rank_column="edgeR_rank",
    )[["id", "edgeR_rank"]]

    comparison = comparison.merge(counts_ranked, on="id", how="left")
    comparison = comparison.merge(edger_ranked, on="id", how="left")
    return comparison.sort_values(["code_0", "code_1"]).reset_index(drop=True)


def compute_method_metrics(
    comparison: pd.DataFrame,
    *,
    score_column: str,
    rank_column: str,
    top_k: int,
) -> dict[str, object]:
    scored = comparison.sort_values(by=score_column, ascending=False, na_position="last").reset_index(drop=True)
    positives = scored[scored["truth_tier"] != "negative"].copy()
    top_hits = scored.head(top_k)
    positive_in_top_k = int((top_hits["truth_tier"] != "negative").sum())
    positive_total = int(len(positives))

    metrics: dict[str, object] = {
        "score_column": score_column,
        "rank_column": rank_column,
        "top_k": int(top_k),
        "positive_total": positive_total,
        "positive_in_top_k": positive_in_top_k,
        "positive_recall_at_k": 0.0 if positive_total == 0 else positive_in_top_k / positive_total,
        "precision_at_k": float((top_hits["truth_tier"] != "negative").mean()),
        "negatives_in_top_k": int((top_hits["truth_tier"] == "negative").sum()),
        "median_positive_rank": float(positives[rank_column].median()),
        "mean_positive_rank": float(positives[rank_column].mean()),
        "tier_recall_at_k": {},
    }

    for tier in ("positive_low", "positive_medium", "positive_high"):
        tier_total = int((scored["truth_tier"] == tier).sum())
        tier_hits = int((top_hits["truth_tier"] == tier).sum())
        metrics["tier_recall_at_k"][tier] = 0.0 if tier_total == 0 else tier_hits / tier_total

    return metrics


def build_top_hits_table(comparison: pd.DataFrame, *, top_n: int) -> pd.DataFrame:
    counts_top = comparison.sort_values("counts_score", ascending=False, na_position="last").head(top_n).reset_index(drop=True)
    edger_top = comparison.sort_values("edgeR_logFC", ascending=False, na_position="last").head(top_n).reset_index(drop=True)

    side_by_side = pd.DataFrame(
        {
            "counts_rank": range(1, len(counts_top) + 1),
            "counts_id": counts_top["id"],
            "counts_tier": counts_top["truth_tier"],
            "counts_score": counts_top["counts_score"],
            "edgeR_rank": range(1, len(edger_top) + 1),
            "edgeR_id": edger_top["id"],
            "edgeR_tier": edger_top["truth_tier"],
            "edgeR_logFC": edger_top["edgeR_logFC"],
            "edgeR_FDR": edger_top.get("edgeR_FDR"),
        }
    )
    return side_by_side


def write_summary_markdown(
    path: Path,
    *,
    dataset_dir: Path,
    truth_metadata: dict[str, object],
    counts_metrics: dict[str, object],
    edger_metrics: dict[str, object],
    top_hits: pd.DataFrame,
) -> None:
    lines = [
        "# Synthetic enrichment comparison",
        "",
        f"- Dataset: `{dataset_dir}`",
        f"- Top K: `{counts_metrics['top_k']}`",
        f"- Baseline distribution: `{truth_metadata['baseline_distribution']}`",
        f"- Dispersion: `{truth_metadata['dispersion']}`",
        f"- Library size CV: `{truth_metadata['library_size_cv']}`",
        "",
        "## counts",
        "",
        f"- Precision at K: `{counts_metrics['precision_at_k']:.3f}`",
        f"- Median positive rank: `{counts_metrics['median_positive_rank']:.2f}`",
        f"- Mean positive rank: `{counts_metrics['mean_positive_rank']:.2f}`",
        f"- Tier recall at K: `{json.dumps(counts_metrics['tier_recall_at_k'], sort_keys=True)}`",
        "",
        "## edgeR",
        "",
        f"- Precision at K: `{edger_metrics['precision_at_k']:.3f}`",
        f"- Median positive rank: `{edger_metrics['median_positive_rank']:.2f}`",
        f"- Mean positive rank: `{edger_metrics['mean_positive_rank']:.2f}`",
        f"- Tier recall at K: `{json.dumps(edger_metrics['tier_recall_at_k'], sort_keys=True)}`",
        "",
        "## Top hits",
        "",
        "```text",
        top_hits.to_string(index=False),
        "```",
        "",
    ]
    path.write_text("\n".join(lines))


def compare_dataset(
    *,
    dataset_dir: Path,
    top_k: int,
    rscript_bin: str,
) -> Path:
    dataset_dir = dataset_dir.expanduser().resolve()
    config_path = dataset_dir / "config.yaml"
    truth_path = dataset_dir / "truth.tsv"
    output_dir = dataset_dir / "comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    truth_metadata, truth = read_truth_table(truth_path)
    analysis_name = DEFAULT_ANALYSIS_NAME
    analysis_root = find_analysis_root(config_path, analysis_name)
    analyse = Analyse()

    counts_dir = run_analysis_method(
        analyse=analyse,
        config_path=config_path,
        analysis_name=analysis_name,
        method="counts",
        analysis_root=analysis_root,
        rscript_bin=rscript_bin,
    )
    edger_dir = run_analysis_method(
        analyse=analyse,
        config_path=config_path,
        analysis_name=analysis_name,
        method="edgeR",
        analysis_root=analysis_root,
        rscript_bin=rscript_bin,
    )

    counts_results = load_counts_results(counts_dir / "stats.csv")
    edger_results = load_edger_results(edger_dir / "stats.csv")
    comparison = build_comparison_table(truth, counts_results, edger_results)

    counts_metrics = compute_method_metrics(
        comparison,
        score_column="counts_score",
        rank_column="counts_rank",
        top_k=top_k,
    )
    edger_metrics = compute_method_metrics(
        comparison,
        score_column="edgeR_logFC",
        rank_column="edgeR_rank",
        top_k=top_k,
    )

    top_hits = build_top_hits_table(comparison, top_n=min(20, top_k))
    comparison_path = output_dir / "comparison.tsv"
    top_hits_path = output_dir / "top_hits.tsv"
    summary_path = output_dir / "summary.md"
    comparison_json_path = output_dir / "comparison.json"

    comparison.to_csv(comparison_path, sep="\t", index=False)
    top_hits.to_csv(top_hits_path, sep="\t", index=False)
    write_summary_markdown(
        summary_path,
        dataset_dir=dataset_dir,
        truth_metadata=truth_metadata,
        counts_metrics=counts_metrics,
        edger_metrics=edger_metrics,
        top_hits=top_hits,
    )
    comparison_json_path.write_text(
        json.dumps(
            {
                "dataset_dir": str(dataset_dir),
                "analysis_root": str(analysis_root),
                "top_k": top_k,
                "truth_metadata": truth_metadata,
                "counts": counts_metrics,
                "edgeR": edger_metrics,
                "top_hits": top_hits.to_dict("records"),
                "files": {
                    "comparison_tsv": str(comparison_path),
                    "top_hits_tsv": str(top_hits_path),
                    "summary_md": str(summary_path),
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    print(f"Wrote comparison artifacts to {output_dir}")
    return output_dir


def main(*, dataset_dir: Path, top_k: int = DEFAULT_TOP_K, rscript_bin: str = "Rscript") -> Path:
    return compare_dataset(dataset_dir=dataset_dir, top_k=top_k, rscript_bin=rscript_bin)


if __name__ == "__main__":
    cli_args = parse_args()
    main(dataset_dir=cli_args.dataset_dir, top_k=cli_args.top_k, rscript_bin=cli_args.rscript_bin)
