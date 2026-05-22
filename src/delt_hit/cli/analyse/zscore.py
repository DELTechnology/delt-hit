"""Compound-level normalized z-score enrichment for the DELT-Hit CLI."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from delt_hit.analyse.enrichment import calculate_fold_enrichment, calculate_zscore
from delt_hit.utils import read_yaml


def get_code_columns(data: pd.DataFrame) -> list[str]:
    """Return DELT-Hit code columns in dataframe order."""
    return [col for col in data.columns if col.startswith("code_")]


def count_observed_compounds(data: pd.DataFrame, code_columns: list[str]) -> int:
    """Count distinct compound identities (code-column combinations) in a count table.

    This is the number of compounds that appeared at least once in sequencing.
    It is NOT the library diversity. Do not use this value as the z-score
    diversity denominator — use infer_diversity_from_cycle_bbs() or a configured
    library_size instead.
    """
    return int(data[code_columns].drop_duplicates().shape[0])


def infer_diversity_from_cycle_bbs(data: pd.DataFrame, code_columns: list[str]) -> int:
    """Estimate theoretical library diversity from per-cycle building-block counts.

    Computes the product of unique building-block identifiers observed per cycle:
    n_bb1 × n_bb2 × ... This is the Cartesian-product size of the designed library
    and is the correct null-hypothesis denominator for the Faver z-score.

    The estimate is exact when every designed building block appears in at least one
    sequenced compound. It underestimates when some building blocks are entirely
    absent from the count file. Set library_size explicitly in the experiment config
    whenever the designed diversity is known.
    """
    diversity = 1
    for col in code_columns:
        diversity *= int(data[col].nunique())
    return diversity


def infer_diversity_from_library_config(
    config: dict,
    code_columns: list[str] | None = None,
) -> int:
    """Calculate designed compound diversity from a DELT-Hit config.

    DELT-Hit demultiplex ``config.yaml`` files store the designed building-block
    families in ``library.building_blocks`` and their members in
    ``whitelists``. The compound-level theoretical diversity is therefore the
    Cartesian product of whitelist lengths for all configured building-block
    families.
    """
    building_blocks = list(config.get("library", {}).get("building_blocks", []))
    if not building_blocks:
        building_blocks = sorted(
            key for key in config.get("whitelists", {}) if key.startswith("B")
        )
    if not building_blocks:
        raise ValueError(
            "Cannot infer library diversity: config has no library.building_blocks "
            "or B* whitelists"
        )
    if code_columns is not None and len(building_blocks) != len(code_columns):
        raise ValueError(
            "Cannot infer library diversity: number of building-block families "
            f"({len(building_blocks)}) does not match count-file code columns "
            f"({len(code_columns)})"
        )

    diversity = 1
    whitelists = config.get("whitelists", {})
    for building_block in building_blocks:
        if building_block not in whitelists:
            raise ValueError(
                "Cannot infer library diversity: missing whitelist for "
                f"{building_block!r}"
            )
        n_building_blocks = len(whitelists[building_block])
        if n_building_blocks < 1:
            raise ValueError(
                "Cannot infer library diversity: empty whitelist for "
                f"{building_block!r}"
            )
        diversity *= n_building_blocks
    return diversity


def resolve_library_size_from_experiment(
    *,
    exp: dict,
    analysis_config_path: Path | None = None,
    code_columns: list[str] | None = None,
) -> int | None:
    """Resolve z-score library size from analysis or demultiplex config.

    Resolution order:

    1. Explicit ``library_size`` in the analysis experiment.
    2. Explicit demultiplex config path under ``library_config_path``,
       ``library_config``, ``demultiplex_config_path``, or
       ``delt_hit_config_path``.
    3. Auto-discovered DELT-Hit ``config.yaml`` in parent directories of the
       configured count files.
    4. ``None`` when no exact source is available; callers may use a fallback.
    """
    configured_size = exp.get("library_size")
    if configured_size is not None:
        return _coerce_library_size(configured_size)

    explicit_config = _get_explicit_library_config_path(exp)
    if explicit_config is not None:
        cfg_path = _resolve_path(explicit_config, analysis_config_path)
        return infer_diversity_from_library_config(
            read_yaml(cfg_path),
            code_columns=code_columns,
        )

    for cfg_path in _discover_library_config_paths(exp, analysis_config_path):
        try:
            return infer_diversity_from_library_config(
                read_yaml(cfg_path),
                code_columns=code_columns,
            )
        except (KeyError, TypeError, ValueError):
            continue
    return None


def compute_zscore(
    data: pd.DataFrame,
    library_size: int,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Compute normalized z-score and fold enrichment per compound/selection."""
    code_columns = get_code_columns(data)
    _validate_scoring_inputs(data=data, code_columns=code_columns, library_size=library_size)

    totals = data.groupby("name", as_index=False)["count"].sum().rename(
        columns={"count": "total_reads"}
    )
    zero_total = totals.loc[totals["total_reads"] <= 0, "name"].tolist()
    if zero_total:
        raise ValueError(
            "Each selection must have total read count > 0; zero total for "
            f"{', '.join(map(str, zero_total))}"
        )

    result = data.merge(totals, on="name", how="left")
    expected_fraction = 1.0 / library_size
    scored_parts = []

    for _, selection_df in result.groupby("name", sort=False):
        counts = selection_df["count"].to_numpy(dtype=float)
        total_reads = int(selection_df["total_reads"].iloc[0])
        zscore, ci_lower, ci_upper = calculate_zscore(
            counts=counts,
            total_reads=total_reads,
            diversity=library_size,
            alpha=alpha,
        )
        fold = calculate_fold_enrichment(
            counts=counts,
            total_reads=total_reads,
            diversity=library_size,
        )

        selection_result = selection_df.copy()
        selection_result["observed_fraction"] = counts / total_reads
        selection_result["expected_fraction"] = expected_fraction
        selection_result["expected_count"] = total_reads / library_size
        selection_result["fold_enrichment"] = fold
        selection_result["zscore"] = zscore
        selection_result["zscore_ci_lower"] = ci_lower
        selection_result["zscore_ci_upper"] = ci_upper
        scored_parts.append(selection_result)

    scored = pd.concat(scored_parts, axis=0, ignore_index=True)
    columns = code_columns + [
        "name",
        "group",
        "count",
        "total_reads",
        "observed_fraction",
        "expected_fraction",
        "expected_count",
        "fold_enrichment",
        "zscore",
        "zscore_ci_lower",
        "zscore_ci_upper",
    ]
    return scored[columns].sort_values(["name", "zscore"], ascending=[True, False])


def zscore_analysis(
    *,
    data: pd.DataFrame,
    samples: pd.DataFrame,
    library_size: int | None,
    save_dir: Path,
    hits_top_n: int = 100,
    alpha: float = 0.05,
) -> None:
    """Run compound-level z-score enrichment and write DELT-Hit CSV outputs."""
    _validate_analysis_inputs(data=data, samples=samples)

    code_columns = get_code_columns(data)
    n_observed = count_observed_compounds(data, code_columns)
    if library_size is None:
        library_size = infer_diversity_from_cycle_bbs(data, code_columns)
        logger.warning(
            "library_size not set in config; inferring theoretical diversity from "
            f"per-cycle building-block counts: {library_size}. "
            "Set library_size explicitly in the experiment config for best results."
        )
    _validate_library_size(library_size)
    if library_size < n_observed:
        raise ValueError(
            "library_size is smaller than the number of observed distinct compounds "
            f"({library_size} < {n_observed}). library_size must equal the total "
            "number of designed compounds in the library, not the observed count."
        )

    save_dir.mkdir(parents=True, exist_ok=True)
    complete_data = _complete_observed_compound_grid(
        data=data,
        samples=samples[["name", "group"]].drop_duplicates(),
        code_columns=code_columns,
    )
    scored = compute_zscore(complete_data, library_size=library_size, alpha=alpha)

    per_selection_columns = code_columns + [
        "count",
        "observed_fraction",
        "expected_fraction",
        "expected_count",
        "fold_enrichment",
        "zscore",
        "zscore_ci_lower",
        "zscore_ci_upper",
    ]
    for selection_name, selection_df in scored.groupby("name", sort=False):
        selection_df[per_selection_columns].reset_index(drop=True).to_csv(
            save_dir / f"{selection_name}.csv",
            index=False,
        )

    stats = _aggregate_group_stats(scored=scored, code_columns=code_columns)
    stats.to_csv(save_dir / "stats.csv", index=False)

    group_names = set(samples["group"].dropna().astype(str))
    sort_column = "protein" if "protein" in group_names else "enrichment"
    stats.sort_values(sort_column, ascending=False).head(hits_top_n).to_csv(
        save_dir / "hits.csv",
        index=False,
    )


def _complete_observed_compound_grid(
    *,
    data: pd.DataFrame,
    samples: pd.DataFrame,
    code_columns: list[str],
) -> pd.DataFrame:
    """Create the observed compound union across every selection."""
    counts = (
        data.groupby(code_columns + ["name"], as_index=False)["count"]
        .sum()
        .astype({"count": float})
    )

    compounds = data[code_columns].drop_duplicates().reset_index(drop=True)
    compound_key = compounds.assign(_merge_key=1)
    sample_key = samples[["name", "group"]].drop_duplicates().assign(_merge_key=1)
    grid = compound_key.merge(sample_key, on="_merge_key").drop(columns="_merge_key")

    complete = grid.merge(counts, on=code_columns + ["name"], how="left")
    complete["count"] = complete["count"].fillna(0.0)
    return complete[code_columns + ["name", "group", "count"]]


def _aggregate_group_stats(scored: pd.DataFrame, code_columns: list[str]) -> pd.DataFrame:
    """Average replicate z-scores by group and compute protein-control contrast."""
    group_mean = (
        scored.groupby(code_columns + ["group"], as_index=False)["zscore"]
        .mean()
        .pivot(index=code_columns, columns="group", values="zscore")
        .reset_index()
    )
    group_mean.columns.name = None
    group_mean = group_mean.fillna(0.0)

    for group in ("protein", "no_protein"):
        if group not in group_mean.columns:
            group_mean[group] = 0.0

    group_mean["enrichment"] = group_mean["protein"] - group_mean["no_protein"]
    preferred_groups = [
        group for group in ("protein", "no_protein") if group in group_mean.columns
    ]
    other_groups = [
        col
        for col in group_mean.columns
        if col not in set(code_columns + preferred_groups + ["enrichment"])
    ]
    return group_mean[code_columns + preferred_groups + other_groups + ["enrichment"]]


def _validate_analysis_inputs(data: pd.DataFrame, samples: pd.DataFrame) -> None:
    code_columns = get_code_columns(data)
    if not code_columns:
        raise ValueError("z-score analysis requires at least one code_* column")

    missing_data_columns = {"name", "count"} - set(data.columns)
    if missing_data_columns:
        raise ValueError(
            f"Count data is missing required columns: {sorted(missing_data_columns)}"
        )

    missing_sample_columns = {"name", "group"} - set(samples.columns)
    if missing_sample_columns:
        raise ValueError(
            f"Sample metadata is missing required columns: {sorted(missing_sample_columns)}"
        )

    unknown_samples = sorted(set(data["name"]) - set(samples["name"]))
    if unknown_samples:
        raise ValueError(
            "Count data contains selections absent from samples metadata: "
            f"{unknown_samples}"
        )

    duplicate_samples = samples.loc[samples["name"].duplicated(), "name"].tolist()
    if duplicate_samples:
        raise ValueError(f"Duplicate sample names in metadata: {duplicate_samples}")

    if data["count"].isna().any():
        raise ValueError("Count data contains missing count values")
    if (data["count"] < 0).any():
        raise ValueError("Count data contains negative count values")
    if data[code_columns].isna().any().any():
        raise ValueError("Count data contains missing code_* values")


def _validate_scoring_inputs(
    *,
    data: pd.DataFrame,
    code_columns: list[str],
    library_size: int,
) -> None:
    if not code_columns:
        raise ValueError("z-score scoring requires at least one code_* column")
    if {"name", "group", "count"} - set(data.columns):
        raise ValueError("z-score scoring requires name, group, and count columns")
    _validate_library_size(library_size)
    if data["count"].isna().any():
        raise ValueError("Count data contains missing count values")
    if (data["count"] < 0).any():
        raise ValueError("Count data contains negative count values")
    if data[code_columns].isna().any().any():
        raise ValueError("Count data contains missing code_* values")


def _validate_library_size(library_size: int) -> None:
    if not isinstance(library_size, (int, np.integer)):
        raise ValueError(f"library_size must be an integer, got {library_size!r}")
    if library_size < 2:
        raise ValueError(
            "library_size must be >= 2 for compound-level z-score scoring; "
            f"got {library_size}"
        )


def _coerce_library_size(value: object) -> int:
    if not isinstance(value, (int, np.integer)):
        raise ValueError(f"library_size must be an integer, got {value!r}")
    return int(value)


def _get_explicit_library_config_path(exp: dict) -> str | Path | None:
    for key in (
        "library_config_path",
        "library_config",
        "demultiplex_config_path",
        "delt_hit_config_path",
    ):
        value = exp.get(key)
        if value is not None:
            return value
    return None


def _resolve_path(path: str | Path, analysis_config_path: Path | None) -> Path:
    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()

    search_roots = []
    if analysis_config_path is not None:
        search_roots.append(Path(analysis_config_path).expanduser().resolve().parent)
    search_roots.append(Path.cwd())

    for root in search_roots:
        resolved = (root / candidate).resolve()
        if resolved.exists():
            return resolved
    return (search_roots[0] / candidate).resolve()


def _discover_library_config_paths(
    exp: dict,
    analysis_config_path: Path | None,
) -> list[Path]:
    paths = []
    seen: set[Path] = set()
    for selection in exp.get("selections", []):
        counts_path = selection.get("counts_path")
        if counts_path is None:
            continue
        resolved_counts_path = _resolve_path(counts_path, analysis_config_path)
        for parent in resolved_counts_path.parents:
            candidate = parent / "config.yaml"
            if candidate.exists() and candidate not in seen:
                seen.add(candidate)
                paths.append(candidate)
    return paths
