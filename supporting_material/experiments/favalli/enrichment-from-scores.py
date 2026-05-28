from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from rich.console import Console
from rich.table import Table

DEFAULT_TOP_N = 10
DEFAULT_PLOT_TOP_NS = [10, 25, 50, 100]
METHOD_ORDER = ["counts", "edgeR", "z_score"]
MAX_CODES_PER_PLOT = 10
TOP_HITS_EXPORT_N = 10
PALETTE = {
    "counts": "#4C78A8",
    "edgeR": "#F58518",
    "z_score": "#54A24B",
}
QUERY_COMPOUNDS = [(159, 474), (357, 474), (68, 474)]
console = Console()


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_analysis_root = script_dir / "lane-1-fasta" / "analysis"
    default_analysis_config = script_dir / "analysis.yaml"

    parser = argparse.ArgumentParser(
        description="Plot top-hit code frequencies from counts, edgeR, and z-score score tables."
    )
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N, help="Number of top compounds to select per method.")
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=default_analysis_root,
        help="Analysis root containing counts/<name>/stats.csv, edgeR/<name>/stats.csv, and z_score/<selection>/stats.csv.",
    )
    parser.add_argument(
        "--analysis-config",
        type=Path,
        default=default_analysis_config,
        help="Analysis YAML used to infer the protein replicate selections for z-score.",
    )
    parser.add_argument("--name", default="ca9_ds", help="Experiment name for the counts and edgeR results.")
    parser.add_argument(
        "--zscore-selections",
        nargs="*",
        default=None,
        help="Optional explicit z-score selection names. Defaults to the protein selections for --name in the analysis config.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for plots and CSV summaries. Defaults to supporting_material/experiments/favalli/enrichment/<name>/.",
    )
    return parser.parse_args()


def validate_top_n(top_n: int) -> None:
    if top_n <= 0:
        raise ValueError("--top-n must be a positive integer.")


def get_code_columns(df: pd.DataFrame) -> list[str]:
    code_columns = [column for column in df.columns if column.startswith("code_")]
    if not code_columns:
        raise ValueError("Expected at least one code_* column.")
    return code_columns


def ensure_id_column(df: pd.DataFrame, code_columns: list[str]) -> pd.DataFrame:
    if "id" in df.columns:
        return df.copy()

    table = df.copy()
    table["id"] = table[code_columns].astype(int).astype(str).agg("_".join, axis=1)
    return table


def get_plot_top_ns(top_n: int) -> list[int]:
    return sorted({*DEFAULT_PLOT_TOP_NS, top_n})


def load_analysis_experiment(*, analysis_config: Path, name: str) -> dict:
    config = yaml.safe_load(Path(analysis_config).read_text())
    experiments = config.get("experiments", [])
    matches = [experiment for experiment in experiments if experiment.get("name") == name]
    if not matches:
        raise ValueError(f"Experiment {name} not found in analysis config {analysis_config}.")
    return matches[0]


def infer_zscore_selections(*, analysis_config: Path, name: str) -> list[str]:
    experiment = load_analysis_experiment(analysis_config=analysis_config, name=name)
    selections = [selection["name"] for selection in experiment["selections"] if selection.get("group") == "protein"]
    if not selections:
        raise ValueError(f"Experiment {name} has no protein selections in {analysis_config}.")
    return selections


def load_counts_scores(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    code_columns = get_code_columns(table)
    table = ensure_id_column(table, code_columns)
    required = {"enrichment"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Counts stats file {path} is missing required columns: {sorted(missing)}")
    ranked = table.sort_values(["enrichment"], ascending=[False]).reset_index(drop=True).copy()
    ranked["rank"] = ranked.index + 1
    ranked["method"] = "counts"
    ranked["score"] = ranked["enrichment"]
    return ranked


def load_edger_scores(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    code_columns = get_code_columns(table)
    table = ensure_id_column(table, code_columns)
    required = {"logFC", "FDR"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"edgeR stats file {path} is missing required columns: {sorted(missing)}")
    ranked = table.loc[(table["FDR"] < 0.05) & (table["logFC"] > 0)].copy()
    ranked = ranked.sort_values(["logFC", "FDR"], ascending=[False, True]).reset_index(drop=True)
    ranked["rank"] = ranked.index + 1
    ranked["method"] = "edgeR"
    ranked["score"] = ranked["logFC"]
    return ranked


def load_zscore_replicates(paths: dict[str, Path]) -> tuple[list[str], dict[str, pd.DataFrame]]:
    tables: dict[str, pd.DataFrame] = {}
    code_columns: list[str] | None = None
    for selection, path in paths.items():
        table = pd.read_csv(path)
        selection_code_columns = get_code_columns(table)
        table = ensure_id_column(table, selection_code_columns)
        required = {"count", "norm_z_score"}
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"z-score stats file {path} is missing required columns: {sorted(missing)}")
        if code_columns is None:
            code_columns = selection_code_columns
        elif selection_code_columns != code_columns:
            raise ValueError("All z-score replicate tables must use the same code columns.")

        ranked = table.sort_values(["norm_z_score", "count"], ascending=[False, False]).reset_index(drop=True).copy()
        ranked["rank"] = ranked.index + 1
        ranked["method"] = "z_score"
        ranked["score"] = ranked["norm_z_score"]
        tables[selection] = ranked

    assert code_columns is not None
    return code_columns, tables


def build_zscore_aggregate(zscore_tables: dict[str, pd.DataFrame], code_columns: list[str]) -> pd.DataFrame:
    merged: pd.DataFrame | None = None
    for selection, table in zscore_tables.items():
        subset = table[["id", *code_columns, "count", "norm_z_score"]].rename(
            columns={
                "count": f"count_{selection}",
                "norm_z_score": f"norm_z_score_{selection}",
            }
        )
        if merged is None:
            merged = subset
        else:
            merged = merged.merge(subset, on=["id", *code_columns], how="outer")

    assert merged is not None
    count_columns = [column for column in merged.columns if column.startswith("count_")]
    score_columns = [column for column in merged.columns if column.startswith("norm_z_score_")]

    merged[count_columns] = merged[count_columns].fillna(0.0)
    merged[score_columns] = merged[score_columns].fillna(0.0)
    merged["mean_count"] = merged[count_columns].mean(axis=1)
    merged["mean_norm_z_score"] = merged[score_columns].mean(axis=1)
    merged["sd_norm_z_score"] = merged[score_columns].std(axis=1, ddof=0)

    ranked = merged.sort_values(["mean_count", "mean_norm_z_score"], ascending=[False, False]).reset_index(drop=True).copy()
    ranked["rank"] = ranked.index + 1
    ranked["method"] = "z_score"
    ranked["score"] = ranked["mean_norm_z_score"]
    return ranked


def compute_frequency_table(ranked_tables: dict[str, pd.DataFrame], top_n: int, code_column: str) -> pd.DataFrame:
    frames = []
    all_codes: set[int] = set()

    for method in ("counts", "edgeR"):
        table = ranked_tables[method].head(top_n)
        all_codes.update(table[code_column].astype(int).tolist())

    zscore_replicates: dict[str, pd.DataFrame] = ranked_tables["z_score_replicates"]  # type: ignore[assignment]
    for table in zscore_replicates.values():
        all_codes.update(table.head(top_n)[code_column].astype(int).tolist())

    ordered_codes = sorted(all_codes)
    if not ordered_codes:
        return pd.DataFrame(columns=["code_value", "method", "frequency", "error", "code_type"])

    for method in ("counts", "edgeR"):
        table = ranked_tables[method].head(top_n)
        counts = table[code_column].value_counts().reindex(ordered_codes, fill_value=0)
        frames.append(
            pd.DataFrame(
                {
                    "code_value": ordered_codes,
                    "method": method,
                    "frequency": counts.to_numpy(dtype=float),
                    "error": np.zeros(len(ordered_codes), dtype=float),
                    "code_type": code_column,
                }
            )
        )

    replicate_counts = []
    for table in zscore_replicates.values():
        replicate_counts.append(
            table.head(top_n)[code_column].value_counts().reindex(ordered_codes, fill_value=0).to_numpy(dtype=float)
        )

    replicate_matrix = np.vstack(replicate_counts) if replicate_counts else np.zeros((0, len(ordered_codes)))
    zscore_mean = replicate_matrix.mean(axis=0) if replicate_counts else np.zeros(len(ordered_codes), dtype=float)
    zscore_std = replicate_matrix.std(axis=0, ddof=0) if replicate_counts else np.zeros(len(ordered_codes), dtype=float)
    frames.append(
        pd.DataFrame(
            {
                "code_value": ordered_codes,
                "method": "z_score",
                "frequency": zscore_mean,
                "error": zscore_std,
                "code_type": code_column,
            }
        )
    )

    freq_df = pd.concat(frames, ignore_index=True)
    top_codes = (
        freq_df.groupby("code_value", as_index=False)["frequency"]
        .sum()
        .sort_values(["frequency", "code_value"], ascending=[False, True])
        .head(MAX_CODES_PER_PLOT)["code_value"]
        .tolist()
    )
    return freq_df.loc[freq_df["code_value"].isin(sorted(top_codes))].copy()


def plot_frequency_table(freq_df: pd.DataFrame, top_n: int, output_path: Path) -> None:
    if freq_df.empty:
        return

    code_type = str(freq_df["code_type"].iat[0])
    code_values = sorted(freq_df["code_value"].unique().tolist())
    x = np.arange(len(code_values), dtype=float)
    width = 0.24

    fig_width = max(8.0, 0.6 * len(code_values) * len(METHOD_ORDER))
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    for idx, method in enumerate(METHOD_ORDER):
        subset = freq_df.loc[freq_df["method"] == method].set_index("code_value").reindex(code_values)
        positions = x + (idx - (len(METHOD_ORDER) - 1) / 2) * width
        ax.bar(
            positions,
            subset["frequency"].to_numpy(dtype=float),
            width=width,
            yerr=subset["error"].to_numpy(dtype=float),
            capsize=4 if method == "z_score" else 0,
            color=PALETTE[method],
            label=method,
        )

    ax.set_title(f"{code_type} frequency among top {top_n} compounds")
    ax.set_xlabel(code_type)
    ax.set_ylabel("Frequency")
    ax.set_xticks(x)
    ax.set_xticklabels([str(code_value) for code_value in code_values])
    ax.legend(title="Method")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    fig.savefig(output_path.with_suffix(".png"), dpi=200)
    plt.close(fig)


def print_method_table(method: str, table: pd.DataFrame, code_columns: list[str]) -> None:
    rich_table = Table(title=f"Top {min(TOP_HITS_EXPORT_N, len(table))} compounds by {method}")
    rich_table.add_column("Rank", justify="right")
    for code_column in code_columns:
        rich_table.add_column(code_column, justify="right")
    rich_table.add_column("Score", justify="right")

    extra_columns: list[str] = []
    if method == "z_score":
        extra_columns = ["mean_count", "sd_norm_z_score"]
        rich_table.add_column("Mean Count", justify="right")
        rich_table.add_column("Score SD", justify="right")

    for _, row in table.head(TOP_HITS_EXPORT_N).iterrows():
        values = [
            str(int(row["rank"])),
            *[str(int(row[code_column])) for code_column in code_columns],
            f"{row['score']:.6g}",
        ]
        if extra_columns:
            values.extend([f"{row['mean_count']:.3f}", f"{row['sd_norm_z_score']:.6g}"])
        rich_table.add_row(*values)

    console.print(rich_table)


def print_query_table(ranked_tables: dict[str, pd.DataFrame], code_columns: list[str]) -> None:
    rich_table = Table(title="Queried compound enrichment by method")
    rich_table.add_column("Compound")
    rich_table.add_column("Method")
    rich_table.add_column("Rank", justify="right")
    rich_table.add_column("Score", justify="right")

    key_columns = code_columns[:2]
    for code_0, code_1 in QUERY_COMPOUNDS:
        compound = f"{code_0};{code_1}"
        for method in METHOD_ORDER:
            table = ranked_tables[method]
            match = table.loc[(table[key_columns[0]] == code_0) & (table[key_columns[1]] == code_1)]
            if match.empty:
                rich_table.add_row(compound, method, "-", "not found")
                continue
            row = match.iloc[0]
            rich_table.add_row(compound, method, str(int(row["rank"])), f"{row['score']:.6g}")

    console.print(rich_table)


def build_controls_table(ranked_tables: dict[str, pd.DataFrame], code_columns: list[str]) -> pd.DataFrame:
    rows = []
    key_columns = code_columns[:2]
    for code_0, code_1 in QUERY_COMPOUNDS:
        compound = f"{code_0};{code_1}"
        for method in METHOD_ORDER:
            table = ranked_tables[method]
            match = table.loc[(table[key_columns[0]] == code_0) & (table[key_columns[1]] == code_1)]
            if match.empty:
                rows.append({"compound": compound, "method": method, "rank": pd.NA, "score": pd.NA})
                continue
            row = match.iloc[0]
            rows.append({"compound": compound, "method": method, "rank": int(row["rank"]), "score": float(row["score"])})
    return pd.DataFrame(rows)


def export_hits_csv(ranked_tables: dict[str, pd.DataFrame], code_columns: list[str], output_path: Path) -> None:
    frames = []
    base_columns = ["method", "rank", *code_columns, "score"]
    for method in METHOD_ORDER:
        table = ranked_tables[method].head(TOP_HITS_EXPORT_N).copy()
        if method == "z_score":
            export_columns = base_columns + ["mean_count", "sd_norm_z_score"]
        else:
            export_columns = base_columns
        frames.append(table.loc[:, export_columns])

    pd.concat(frames, ignore_index=True).to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()
    validate_top_n(args.top_n)

    analysis_root = args.analysis_root.expanduser().resolve()
    analysis_config = args.analysis_config.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else Path(__file__).resolve().parent / "enrichment" / args.name
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    zscore_selections = args.zscore_selections or infer_zscore_selections(analysis_config=analysis_config, name=args.name)

    counts_path = analysis_root / "counts" / args.name / "stats.csv"
    edger_path = analysis_root / "edgeR" / args.name / "stats.csv"
    zscore_paths = {
        selection: analysis_root / "z_score" / selection / "stats.csv"
        for selection in zscore_selections
    }

    counts_ranked = load_counts_scores(counts_path)
    edger_ranked = load_edger_scores(edger_path)
    code_columns, zscore_replicates = load_zscore_replicates(zscore_paths)
    zscore_ranked = build_zscore_aggregate(zscore_replicates, code_columns)

    ranked_tables: dict[str, pd.DataFrame] = {
        "counts": counts_ranked,
        "edgeR": edger_ranked,
        "z_score": zscore_ranked,
        "z_score_replicates": zscore_replicates,  # type: ignore[dict-item]
    }

    for top_n in get_plot_top_ns(args.top_n):
        for method in METHOD_ORDER:
            print(f"{method}: selected {len(ranked_tables[method].head(top_n))} compounds for top {top_n}")

        for code_column in code_columns:
            freq_df = compute_frequency_table(ranked_tables=ranked_tables, top_n=top_n, code_column=code_column)
            if freq_df.empty:
                continue
            output_path = output_dir / f"top_{top_n}_{code_column}_frequency.pdf"
            plot_frequency_table(freq_df=freq_df, top_n=top_n, output_path=output_path)
            print(f"Wrote {output_path}")

    print_query_table(ranked_tables, code_columns)
    for method in METHOD_ORDER:
        print_method_table(method, ranked_tables[method], code_columns)

    hits_path = output_dir / "hits.csv"
    export_hits_csv(ranked_tables=ranked_tables, code_columns=code_columns, output_path=hits_path)
    print(f"Wrote {hits_path}")

    controls_path = output_dir / "controls.csv"
    build_controls_table(ranked_tables=ranked_tables, code_columns=code_columns).to_csv(controls_path, index=False)
    print(f"Wrote {controls_path}")


if __name__ == "__main__":
    main()
