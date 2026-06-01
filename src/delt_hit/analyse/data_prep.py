from __future__ import annotations

from pathlib import Path

import pandas as pd


def canonicalize_group_name(group: str) -> str:
    """Map legacy group labels onto generic control/condition labels."""
    mapping = {
        "protein": "condition",
        "no_protein": "control",
    }
    return mapping.get(group, group)


def prepare_replicate_analysis_data(*, exp: dict, data_path: Path, samples_path: Path) -> tuple[Path, Path]:
    """Compile counts and sample metadata for counts/edgeR analysis."""
    selections = exp["selections"]
    data = []
    for selection in selections:
        counts_path = Path(selection["counts_path"]).expanduser().resolve()
        if not counts_path.exists():
            raise FileNotFoundError(f"Counts file for selection {selection} not found at {counts_path}")

        counts = pd.read_csv(counts_path, delimiter="\t")
        counts["name"] = selection["name"]
        data.append(counts)

    samples = pd.DataFrame(selections)[["name", "group"]]
    samples["group"] = samples["group"].map(canonicalize_group_name)
    samples_path.parent.mkdir(parents=True, exist_ok=True)
    samples.to_csv(samples_path, index=False)

    merged = pd.concat(data, axis=0)
    data_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(data_path, index=False)
    return data_path, samples_path


def load_zscore_counts(*, counts_path: Path) -> pd.DataFrame:
    """Load and validate one counts table for z-score analysis."""
    counts_path = Path(counts_path).expanduser().resolve()
    if not counts_path.exists():
        raise FileNotFoundError(f"Counts file not found at {counts_path}")

    table = pd.read_csv(counts_path, sep="\t")
    code_columns = [column for column in table.columns if column.startswith("code_")]
    if not code_columns:
        raise ValueError("Counts file must contain at least one `code_*` column.")
    if "count" not in table.columns:
        raise ValueError("Counts file must contain a `count` column.")
    if float(table["count"].sum()) <= 0:
        raise ValueError("Counts file must have a positive total count.")
    return table
