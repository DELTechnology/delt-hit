"""Normalized z-score enrichment scoring for DEL selections.

Implements the depth-independent z-score from Faver et al. (2019):

    z_n = (p_o - p_i) / sqrt(p_i * (1 - p_i))

where p_o = count / total_reads (observed fraction per selection) and
p_i = 1 / library_size (expected uniform fraction).

Reference:
    Faver et al., ACS Combinatorial Science 2019, 21(2), 75–82.
    DOI: 10.1021/acscombsci.8b00116
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


def derive_library_size(data: pd.DataFrame, code_columns: list[str]) -> int:
    """Count distinct compound combinations across all selections.

    Args:
        data: Merged counts dataframe with code_* and name columns.
        code_columns: List of code column names (e.g. ['code_0', 'code_1']).

    Returns:
        Number of unique compound combinations observed in data.
    """
    return int(data[code_columns].drop_duplicates().shape[0])


def compute_zscore(
    data: pd.DataFrame,
    library_size: int,
) -> pd.DataFrame:
    """Compute the normalized z-score for each compound in each selection.

    Args:
        data: Long-format dataframe with columns code_*, name, group, count.
        library_size: Total number of distinct compounds in the library.
            Used to derive the expected uniform fraction p_i = 1 / library_size.

    Returns:
        Dataframe with columns code_*, name, group, count,
        observed_fraction, expected_fraction, zscore.
        Rows are sorted by name then descending zscore.
    """
    code_columns = [c for c in data.columns if c.startswith("code_")]

    p_i = 1.0 / library_size
    denominator = math.sqrt(p_i * (1.0 - p_i))

    total_reads = (
        data.groupby("name")["count"]
        .sum()
        .rename("total_reads")
        .reset_index()
    )
    result = data.merge(total_reads, on="name", how="left")

    result = result.assign(
        observed_fraction=result["count"] / result["total_reads"],
        expected_fraction=p_i,
    )
    result = result.assign(
        zscore=(result["observed_fraction"] - p_i) / denominator,
    )

    col_order = code_columns + ["name", "group", "count", "observed_fraction", "expected_fraction", "zscore"]
    return result[col_order].sort_values(["name", "zscore"], ascending=[True, False]).reset_index(drop=True)


def zscore_analysis(
    *,
    data: pd.DataFrame,
    samples: pd.DataFrame,
    library_size: int | None,
    save_dir: Path,
    hits_top_n: int = 100,
) -> None:
    """Run z-score enrichment analysis and write outputs to save_dir.

    Outputs mirror the structure of the counts/ method:
      - stats.csv         — all compounds with per-group mean z-score and enrichment
      - hits.csv          — top hits_top_n rows sorted by descending protein zscore
      - <selection>.csv   — per-selection: code_*, count, observed_fraction, expected_fraction, zscore

    Args:
        data: Merged long-format counts (code_*, name, count columns).
        samples: Sample metadata with name and group columns.
        library_size: Override library size. If None, derived from observed data.
        save_dir: Directory to write outputs into.
        hits_top_n: Number of top hits to write to hits.csv.
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    data_with_group = data.merge(samples, on="name", how="left")
    code_columns = [c for c in data_with_group.columns if c.startswith("code_")]

    if library_size is None:
        library_size = derive_library_size(data_with_group, code_columns)

    scored = compute_zscore(data_with_group, library_size)

    # Per-selection exports
    for selection_name in scored["name"].unique():
        sel_df = scored[scored["name"] == selection_name][
            code_columns + ["count", "observed_fraction", "expected_fraction", "zscore"]
        ].reset_index(drop=True)
        sel_df.to_csv(save_dir / f"{selection_name}.csv", index=False)

    # Aggregate by compound and group: mean z-score across replicates
    group_mean = (
        scored.groupby(code_columns + ["group"])["zscore"]
        .mean()
        .reset_index()
        .pivot(index=code_columns, columns="group", values="zscore")
        .reset_index()
    )
    group_mean.columns.name = None

    # Ensure protein and no_protein columns exist (fill with 0 if a group is absent)
    for col in ("protein", "no_protein"):
        if col not in group_mean.columns:
            group_mean[col] = 0.0

    group_mean = group_mean.fillna(0.0)
    group_mean["enrichment"] = group_mean["protein"] - group_mean["no_protein"]

    # Reorder: code_* then protein, no_protein, enrichment (parallel to counts/stats.csv)
    present_groups = [g for g in ("protein", "no_protein") if g in group_mean.columns]
    other_groups = [g for g in group_mean.columns if g not in code_columns + present_groups + ["enrichment"]]
    col_order = code_columns + present_groups + other_groups + ["enrichment"]
    stats = group_mean[col_order]

    stats.to_csv(save_dir / "stats.csv", index=False)

    # Hits: top N by descending protein zscore
    sort_col = "protein" if "protein" in stats.columns else "enrichment"
    hits = stats.sort_values(sort_col, ascending=False).head(hits_top_n)
    hits.to_csv(save_dir / "hits.csv", index=False)
