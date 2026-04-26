#!/usr/bin/env python3
"""Generate synthetic count tables for counts-vs-edgeR enrichment benchmarks."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "target" / "synthetic_enrichment"
DEFAULT_EXPERIMENT_NAME = "synthetic_enrichment"
DEFAULT_ANALYSIS_NAME = "protein_vs_no_protein"
DEFAULT_BUILDING_BLOCKS_PER_CYCLE = 100
DEFAULT_POSITIVE_LOW = 50
DEFAULT_POSITIVE_MEDIUM = 25
DEFAULT_POSITIVE_HIGH = 10
DEFAULT_FC_LOW = 1.5
DEFAULT_FC_MEDIUM = 3.0
DEFAULT_FC_HIGH = 8.0
DEFAULT_BASELINE_DISTRIBUTION = "lognormal"
DEFAULT_BASELINE_LOGMEAN = 2.5
DEFAULT_BASELINE_LOGSIGMA = 1.0
DEFAULT_BASELINE_GAMMA_SHAPE = 2.0
DEFAULT_BASELINE_GAMMA_SCALE = 10.0
DEFAULT_DISPERSION = 0.20
DEFAULT_LIBRARY_SIZE_CV = 0.10
DEFAULT_SEED = 13
DEFAULT_NUM_REPLICATES = 3
TRUTH_METADATA_PREFIX = "# meta\t"

GROUPS = ("no_protein", "protein")
TIERS = ("negative", "positive_low", "positive_medium", "positive_high")


def build_compound_table(*, building_blocks_per_cycle: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for code_0 in range(building_blocks_per_cycle):
        for code_1 in range(building_blocks_per_cycle):
            rows.append(
                {
                    "code_0": code_0,
                    "code_1": code_1,
                    "id": f"{code_0}_{code_1}",
                    "building_block_id_1": f"BB1_{code_0 + 1:03d}",
                    "building_block_id_2": f"BB2_{code_1 + 1:03d}",
                }
            )
    return pd.DataFrame(rows)


def validate_positive_counts(
    *,
    building_blocks_per_cycle: int,
    positive_low: int,
    positive_medium: int,
    positive_high: int,
) -> None:
    total_compounds = building_blocks_per_cycle**2
    total_positives = positive_low + positive_medium + positive_high
    if total_positives >= total_compounds:
        raise ValueError(
            "Requested positive controls must leave at least one negative compound. "
            f"Got {total_positives} positives for {total_compounds} compounds."
        )


def sample_baseline_means(
    rng: np.random.Generator,
    *,
    size: int,
    baseline_distribution: str,
    baseline_logmean: float,
    baseline_logsigma: float,
    baseline_gamma_shape: float,
    baseline_gamma_scale: float,
) -> np.ndarray:
    if baseline_distribution == "lognormal":
        return rng.lognormal(mean=baseline_logmean, sigma=baseline_logsigma, size=size)
    if baseline_distribution == "gamma":
        return rng.gamma(shape=baseline_gamma_shape, scale=baseline_gamma_scale, size=size)
    raise ValueError(f"Unsupported baseline distribution: {baseline_distribution}")


def assign_truth_tiers(
    rng: np.random.Generator,
    *,
    num_compounds: int,
    positive_low: int,
    positive_medium: int,
    positive_high: int,
) -> tuple[np.ndarray, np.ndarray]:
    tiers = np.full(num_compounds, "negative", dtype=object)
    fold_changes = np.ones(num_compounds, dtype=float)
    shuffled = rng.permutation(num_compounds)

    boundaries = {
        "positive_high": positive_high,
        "positive_medium": positive_medium,
        "positive_low": positive_low,
    }
    start = 0
    for tier_name, count in boundaries.items():
        stop = start + count
        indices = shuffled[start:stop]
        tiers[indices] = tier_name
        start = stop

    fold_changes[tiers == "positive_low"] = DEFAULT_FC_LOW
    fold_changes[tiers == "positive_medium"] = DEFAULT_FC_MEDIUM
    fold_changes[tiers == "positive_high"] = DEFAULT_FC_HIGH
    return tiers, fold_changes


def sample_library_size_factors(
    rng: np.random.Generator,
    *,
    sample_names: list[str],
    library_size_cv: float,
) -> dict[str, float]:
    if library_size_cv == 0:
        return {sample_name: 1.0 for sample_name in sample_names}

    sigma_sq = math.log1p(library_size_cv**2)
    sigma = math.sqrt(sigma_sq)
    mean = -sigma_sq / 2.0
    sampled = rng.lognormal(mean=mean, sigma=sigma, size=len(sample_names))
    return {sample_name: float(value) for sample_name, value in zip(sample_names, sampled, strict=True)}


def negative_binomial_draw(
    rng: np.random.Generator,
    *,
    means: np.ndarray,
    dispersion: float,
) -> np.ndarray:
    means = np.asarray(means, dtype=float)
    if dispersion < 0:
        raise ValueError("dispersion must be non-negative")
    if dispersion == 0:
        return rng.poisson(lam=means).astype(np.int64)

    clipped_means = np.clip(means, a_min=1e-12, a_max=None)
    size = 1.0 / dispersion
    probability = size / (size + clipped_means)
    return rng.negative_binomial(size, probability).astype(np.int64)


def sample_counts_for_replicates(
    truth: pd.DataFrame,
    *,
    rng: np.random.Generator,
    dispersion: float,
    library_size_factors: dict[str, float],
) -> dict[str, pd.DataFrame]:
    result: dict[str, pd.DataFrame] = {}
    for sample_name, library_factor in library_size_factors.items():
        group = "protein" if sample_name.startswith("protein") else "no_protein"
        mean_column = "protein_mean" if group == "protein" else "no_protein_mean"
        replicate_means = truth[mean_column].to_numpy(dtype=float) * library_factor
        counts = negative_binomial_draw(rng, means=replicate_means, dispersion=dispersion)
        result[sample_name] = truth[["code_0", "code_1", "id"]].assign(count=counts)[
            ["code_0", "code_1", "count", "id"]
        ]
    return result


def build_truth_table(
    *,
    building_blocks_per_cycle: int,
    positive_low: int,
    positive_medium: int,
    positive_high: int,
    fc_low: float,
    fc_medium: float,
    fc_high: float,
    baseline_distribution: str,
    baseline_logmean: float,
    baseline_logsigma: float,
    baseline_gamma_shape: float,
    baseline_gamma_scale: float,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    validate_positive_counts(
        building_blocks_per_cycle=building_blocks_per_cycle,
        positive_low=positive_low,
        positive_medium=positive_medium,
        positive_high=positive_high,
    )

    rng = np.random.default_rng(seed)
    truth = build_compound_table(building_blocks_per_cycle=building_blocks_per_cycle)
    baseline_means = sample_baseline_means(
        rng,
        size=len(truth),
        baseline_distribution=baseline_distribution,
        baseline_logmean=baseline_logmean,
        baseline_logsigma=baseline_logsigma,
        baseline_gamma_shape=baseline_gamma_shape,
        baseline_gamma_scale=baseline_gamma_scale,
    )
    tiers, fold_changes = assign_truth_tiers(
        rng,
        num_compounds=len(truth),
        positive_low=positive_low,
        positive_medium=positive_medium,
        positive_high=positive_high,
    )

    fold_change_map = {
        "negative": 1.0,
        "positive_low": fc_low,
        "positive_medium": fc_medium,
        "positive_high": fc_high,
    }
    fold_changes = np.array([fold_change_map[tier] for tier in tiers], dtype=float)

    truth = truth.assign(
        truth_tier=tiers,
        expected_direction=np.where(tiers == "negative", "none", "protein_gt_no_protein"),
        baseline_mean=baseline_means,
        no_protein_mean=baseline_means,
        protein_mean=baseline_means * fold_changes,
        target_fold_change=fold_changes,
    )

    metadata = {
        "baseline_distribution": baseline_distribution,
        "baseline_logmean": baseline_logmean,
        "baseline_logsigma": baseline_logsigma,
        "baseline_gamma_shape": baseline_gamma_shape,
        "baseline_gamma_scale": baseline_gamma_scale,
        "building_blocks_per_cycle": building_blocks_per_cycle,
        "positive_low": positive_low,
        "positive_medium": positive_medium,
        "positive_high": positive_high,
        "fc_low": fc_low,
        "fc_medium": fc_medium,
        "fc_high": fc_high,
        "seed": seed,
    }
    return truth, metadata


def truth_metadata_lines(metadata: dict[str, object]) -> list[str]:
    lines = []
    for key in sorted(metadata):
        lines.append(f"{TRUTH_METADATA_PREFIX}{key}\t{json.dumps(metadata[key], sort_keys=True)}\n")
    return lines


def write_truth_table(path: Path, truth: pd.DataFrame, metadata: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.writelines(truth_metadata_lines(metadata))
        truth.to_csv(handle, sep="\t", index=False)


def read_truth_table(path: Path) -> tuple[dict[str, object], pd.DataFrame]:
    metadata: dict[str, object] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(TRUTH_METADATA_PREFIX):
                continue
            _, key, raw_value = line.rstrip("\n").split("\t", 2)
            metadata[key] = json.loads(raw_value)
    truth = pd.read_csv(path, sep="\t", comment="#")
    return metadata, truth


def selection_names() -> list[str]:
    selections: list[str] = []
    for group in GROUPS:
        for replicate in range(1, DEFAULT_NUM_REPLICATES + 1):
            selections.append(f"{group}_rep{replicate}")
    return selections


def selection_group(selection_name: str) -> str:
    return "protein" if selection_name.startswith("protein") else "no_protein"


def write_count_tables(path: Path, counts_by_selection: dict[str, pd.DataFrame]) -> dict[str, str]:
    selection_paths: dict[str, str] = {}
    for selection_name, table in counts_by_selection.items():
        counts_path = path / "selections" / selection_name / "counts.txt"
        counts_path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(counts_path, sep="\t", index=False)
        selection_paths[selection_name] = str(counts_path.resolve())
    return selection_paths


def write_config(
    path: Path,
    *,
    analysis_name: str,
    analysis_save_dir: Path,
    selection_paths: dict[str, str],
) -> None:
    config = {
        "experiments": [
            {
                "name": analysis_name,
                "save_dir": str(analysis_save_dir.resolve()),
                "selections": [
                    {
                        "name": selection_name,
                        "counts_path": selection_paths[selection_name],
                        "group": selection_group(selection_name),
                    }
                    for selection_name in selection_names()
                ],
            }
        ]
    }
    path.write_text(yaml.safe_dump(config, sort_keys=False))


def write_manifest(
    path: Path,
    *,
    experiment_dir: Path,
    experiment_name: str,
    analysis_name: str,
    config_path: Path,
    truth_path: Path,
    selection_paths: dict[str, str],
    library_size_factors: dict[str, float],
    generation_parameters: dict[str, object],
) -> None:
    manifest = {
        "experiment_name": experiment_name,
        "analysis_name": analysis_name,
        "output_dir": str(experiment_dir.resolve()),
        "config_path": str(config_path.resolve()),
        "truth_path": str(truth_path.resolve()),
        "selection_paths": selection_paths,
        "library_size_factors": library_size_factors,
        "generation_parameters": generation_parameters,
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--building-blocks-per-cycle",
        type=int,
        default=DEFAULT_BUILDING_BLOCKS_PER_CYCLE,
        help="Number of building blocks for each of the two cycles.",
    )
    parser.add_argument("--positive-low", type=int, default=DEFAULT_POSITIVE_LOW)
    parser.add_argument("--positive-medium", type=int, default=DEFAULT_POSITIVE_MEDIUM)
    parser.add_argument("--positive-high", type=int, default=DEFAULT_POSITIVE_HIGH)
    parser.add_argument("--fc-low", type=float, default=DEFAULT_FC_LOW)
    parser.add_argument("--fc-medium", type=float, default=DEFAULT_FC_MEDIUM)
    parser.add_argument("--fc-high", type=float, default=DEFAULT_FC_HIGH)
    parser.add_argument(
        "--baseline-distribution",
        choices=("lognormal", "gamma"),
        default=DEFAULT_BASELINE_DISTRIBUTION,
        help="Distribution used to sample heterogeneous per-compound background means.",
    )
    parser.add_argument("--baseline-logmean", type=float, default=DEFAULT_BASELINE_LOGMEAN)
    parser.add_argument("--baseline-logsigma", type=float, default=DEFAULT_BASELINE_LOGSIGMA)
    parser.add_argument("--baseline-gamma-shape", type=float, default=DEFAULT_BASELINE_GAMMA_SHAPE)
    parser.add_argument("--baseline-gamma-scale", type=float, default=DEFAULT_BASELINE_GAMMA_SCALE)
    parser.add_argument("--dispersion", type=float, default=DEFAULT_DISPERSION)
    parser.add_argument("--library-size-cv", type=float, default=DEFAULT_LIBRARY_SIZE_CV)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--experiment-name", default=DEFAULT_EXPERIMENT_NAME)
    return parser.parse_args()


def main(
    *,
    building_blocks_per_cycle: int = DEFAULT_BUILDING_BLOCKS_PER_CYCLE,
    positive_low: int = DEFAULT_POSITIVE_LOW,
    positive_medium: int = DEFAULT_POSITIVE_MEDIUM,
    positive_high: int = DEFAULT_POSITIVE_HIGH,
    fc_low: float = DEFAULT_FC_LOW,
    fc_medium: float = DEFAULT_FC_MEDIUM,
    fc_high: float = DEFAULT_FC_HIGH,
    baseline_distribution: str = DEFAULT_BASELINE_DISTRIBUTION,
    baseline_logmean: float = DEFAULT_BASELINE_LOGMEAN,
    baseline_logsigma: float = DEFAULT_BASELINE_LOGSIGMA,
    baseline_gamma_shape: float = DEFAULT_BASELINE_GAMMA_SHAPE,
    baseline_gamma_scale: float = DEFAULT_BASELINE_GAMMA_SCALE,
    dispersion: float = DEFAULT_DISPERSION,
    library_size_cv: float = DEFAULT_LIBRARY_SIZE_CV,
    seed: int = DEFAULT_SEED,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    experiment_name: str = DEFAULT_EXPERIMENT_NAME,
) -> Path:
    if building_blocks_per_cycle < 2:
        raise ValueError("building_blocks_per_cycle must be at least 2")
    if dispersion < 0:
        raise ValueError("dispersion must be non-negative")
    if library_size_cv < 0:
        raise ValueError("library_size_cv must be non-negative")

    experiment_dir = Path(output_dir) / experiment_name
    experiment_dir.mkdir(parents=True, exist_ok=True)

    truth, truth_metadata = build_truth_table(
        building_blocks_per_cycle=building_blocks_per_cycle,
        positive_low=positive_low,
        positive_medium=positive_medium,
        positive_high=positive_high,
        fc_low=fc_low,
        fc_medium=fc_medium,
        fc_high=fc_high,
        baseline_distribution=baseline_distribution,
        baseline_logmean=baseline_logmean,
        baseline_logsigma=baseline_logsigma,
        baseline_gamma_shape=baseline_gamma_shape,
        baseline_gamma_scale=baseline_gamma_scale,
        seed=seed,
    )

    sample_names = selection_names()
    sampling_rng = np.random.default_rng(seed + 1)
    library_size_factors = sample_library_size_factors(
        sampling_rng,
        sample_names=sample_names,
        library_size_cv=library_size_cv,
    )
    counts_by_selection = sample_counts_for_replicates(
        truth,
        rng=sampling_rng,
        dispersion=dispersion,
        library_size_factors=library_size_factors,
    )

    selection_paths = write_count_tables(experiment_dir, counts_by_selection)
    truth_path = experiment_dir / "truth.tsv"
    truth_metadata = {
        **truth_metadata,
        "dispersion": dispersion,
        "library_size_cv": library_size_cv,
    }
    write_truth_table(truth_path, truth, truth_metadata)

    config_path = experiment_dir / "config.yaml"
    analysis_save_dir = experiment_dir / "analysis"
    write_config(
        config_path,
        analysis_name=DEFAULT_ANALYSIS_NAME,
        analysis_save_dir=analysis_save_dir,
        selection_paths=selection_paths,
    )

    manifest_path = experiment_dir / "manifest.json"
    write_manifest(
        manifest_path,
        experiment_dir=experiment_dir,
        experiment_name=experiment_name,
        analysis_name=DEFAULT_ANALYSIS_NAME,
        config_path=config_path,
        truth_path=truth_path,
        selection_paths=selection_paths,
        library_size_factors=library_size_factors,
        generation_parameters={
            **truth_metadata,
            "num_replicates_per_group": DEFAULT_NUM_REPLICATES,
            "groups": list(GROUPS),
        },
    )

    print(f"Wrote synthetic enrichment dataset to {experiment_dir}")
    print(f"Config: {config_path}")
    print(f"Truth:  {truth_path}")
    print(f"Manifest: {manifest_path}")
    return experiment_dir


if __name__ == "__main__":
    cli_args = parse_args()
    main(
        building_blocks_per_cycle=cli_args.building_blocks_per_cycle,
        positive_low=cli_args.positive_low,
        positive_medium=cli_args.positive_medium,
        positive_high=cli_args.positive_high,
        fc_low=cli_args.fc_low,
        fc_medium=cli_args.fc_medium,
        fc_high=cli_args.fc_high,
        baseline_distribution=cli_args.baseline_distribution,
        baseline_logmean=cli_args.baseline_logmean,
        baseline_logsigma=cli_args.baseline_logsigma,
        baseline_gamma_shape=cli_args.baseline_gamma_shape,
        baseline_gamma_scale=cli_args.baseline_gamma_scale,
        dispersion=cli_args.dispersion,
        library_size_cv=cli_args.library_size_cv,
        seed=cli_args.seed,
        output_dir=cli_args.output_dir,
        experiment_name=cli_args.experiment_name,
    )
