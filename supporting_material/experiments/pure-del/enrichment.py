from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rich.console import Console
from rich.table import Table

PSEUDOCOUNT = 1e-6
DEFAULT_TOP_N = 10
DEFAULT_PLOT_TOP_NS = [10, 50, 100, 500, 1000]
METHOD_ORDER = ["diff", "ratio", "edgeR"]
MAX_CODES_PER_PLOT = 10
TOP_HITS_EXPORT_N = 10
PALETTE = {
    "diff": "#4C78A8",
    "ratio": "#F58518",
    "edgeR": "#54A24B",
}
console = Console()


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_data_dir = script_dir / "lane-2" / "analysis" / "his_pure_up"
    default_output_dir = script_dir / "enrichment" / "his_pure_up"

    parser = argparse.ArgumentParser(
        description="Plot code frequencies among top-ranked Pure-DEL compounds."
    )
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N, help="Number of top compounds to select per method.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=default_data_dir,
        help="Analysis directory containing counts/stats.csv and edgeR/enrichment_stats.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir,
        help="Directory where plots and CSV summaries will be written.",
    )
    return parser.parse_args()


def validate_top_n(top_n: int) -> None:
    if top_n <= 0:
        raise ValueError("--top-n must be a positive integer.")


def get_code_columns(df: pd.DataFrame) -> list[str]:
    return [column for column in df.columns if column.startswith("code_")]


def get_plot_top_ns(top_n: int) -> list[int]:
    return sorted({*DEFAULT_PLOT_TOP_NS, top_n})


def build_ranked_tables(counts_df: pd.DataFrame, edger_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    counts_ranked = counts_df.assign(
        diff=counts_df["protein"] - counts_df["no_protein"],
        ratio=(counts_df["protein"] + PSEUDOCOUNT) / (counts_df["no_protein"] + PSEUDOCOUNT),
    )

    ranked_tables = {
        "diff": counts_ranked.sort_values(["diff", "protein", "no_protein"], ascending=[False, False, True]).copy(),
        "ratio": counts_ranked.sort_values(["ratio", "protein", "no_protein"], ascending=[False, False, True]).copy(),
    }

    edger_hits = edger_df.loc[(edger_df["FDR"] < 0.05) & (edger_df["logFC"] > 0)].copy()
    ranked_tables["edgeR"] = edger_hits.sort_values(["logFC", "FDR"], ascending=[False, True]).copy()

    score_columns = {"diff": "diff", "ratio": "ratio", "edgeR": "logFC"}
    for method, table in ranked_tables.items():
        ranked_tables[method] = table.reset_index(drop=True).copy()
        ranked_tables[method]["rank"] = ranked_tables[method].index + 1
        ranked_tables[method]["method"] = method
        ranked_tables[method]["score"] = ranked_tables[method][score_columns[method]]

    return ranked_tables


def compute_frequency_table(top_tables: dict[str, pd.DataFrame], code_column: str) -> pd.DataFrame:
    frames = []
    all_codes: set[int] = set()

    for method in METHOD_ORDER:
        table = top_tables[method]
        if table.empty:
            continue
        all_codes.update(table[code_column].astype(int).tolist())

    if not all_codes:
        return pd.DataFrame(columns=["code_value", "frequency", "method", "code_type"])

    ordered_codes = sorted(all_codes)

    for method in METHOD_ORDER:
        table = top_tables[method]
        if table.empty:
            counts = pd.Series(0, index=ordered_codes, dtype=int)
        else:
            counts = table[code_column].value_counts().reindex(ordered_codes, fill_value=0)
        frames.append(
            pd.DataFrame(
                {
                    "code_value": ordered_codes,
                    "frequency": counts.to_numpy(),
                    "method": method,
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
    ordered_top_codes = sorted(top_codes)
    return freq_df.loc[freq_df["code_value"].isin(ordered_top_codes)].copy()


def plot_frequency_table(freq_df: pd.DataFrame, top_n: int, output_path: Path) -> None:
    if freq_df.empty:
        return

    code_type = freq_df["code_type"].iat[0]
    fig_width = max(8, 0.45 * freq_df["code_value"].nunique() * len(METHOD_ORDER))
    plt.figure(figsize=(fig_width, 6))
    ax = sns.barplot(
        data=freq_df,
        x="code_value",
        y="frequency",
        hue="method",
        hue_order=METHOD_ORDER,
        palette=PALETTE,
    )
    ax.set_title(f"{code_type} frequency among top {top_n} compounds")
    ax.set_xlabel(code_type)
    ax.set_ylabel("Frequency")
    ax.legend(title="Method")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.savefig(output_path.with_suffix(".png"), dpi=200)
    plt.close()


def print_method_table(method: str, table: pd.DataFrame, code_columns: list[str]) -> None:
    rich_table = Table(title=f"Top {min(TOP_HITS_EXPORT_N, len(table))} enriched compounds by {method}")
    rich_table.add_column("Rank", justify="right")
    for code_column in code_columns:
        rich_table.add_column(code_column, justify="right")
    rich_table.add_column("Score", justify="right")

    for _, row in table.head(TOP_HITS_EXPORT_N).iterrows():
        rich_table.add_row(
            str(int(row["rank"])),
            *[str(int(row[code_column])) for code_column in code_columns],
            f"{row['score']:.6g}",
        )

    console.print(rich_table)


def export_hits_csv(ranked_tables: dict[str, pd.DataFrame], code_columns: list[str], output_path: Path) -> None:
    frames = []
    columns = ["method", "rank", *code_columns, "score"]
    for method in METHOD_ORDER:
        table = ranked_tables[method].head(TOP_HITS_EXPORT_N).copy()
        frames.append(table.loc[:, columns])

    pd.concat(frames, ignore_index=True).to_csv(output_path, index=False)


def export_summary_csv(ranked_tables: dict[str, pd.DataFrame], code_columns: list[str], output_path: Path) -> None:
    rows = []
    for method in METHOD_ORDER:
        table = ranked_tables[method]
        top_row = table.iloc[0] if not table.empty else None
        row = {
            "method": method,
            "selected_compounds": int(len(table)),
            "top_score": float(top_row["score"]) if top_row is not None else pd.NA,
            "top_rank": int(top_row["rank"]) if top_row is not None else pd.NA,
        }
        for code_column in code_columns:
            row[f"top_{code_column}"] = int(top_row[code_column]) if top_row is not None else pd.NA
        rows.append(row)

    pd.DataFrame(rows).to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()
    validate_top_n(args.top_n)

    data_dir = args.data_dir.expanduser().resolve()
    counts_path = data_dir / "counts" / "stats.csv"
    edger_path = data_dir / "edgeR" / "enrichment_stats.csv"
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    counts_df = pd.read_csv(counts_path)
    edger_df = pd.read_csv(edger_path)
    code_columns = get_code_columns(counts_df)

    ranked_tables = build_ranked_tables(counts_df=counts_df, edger_df=edger_df)
    for top_n in get_plot_top_ns(args.top_n):
        top_tables = {method: table.head(top_n).copy() for method, table in ranked_tables.items()}

        for method in METHOD_ORDER:
            print(f"{method}: selected {len(top_tables[method])} compounds for top {top_n}")

        for code_column in code_columns:
            freq_df = compute_frequency_table(top_tables=top_tables, code_column=code_column)
            if freq_df.empty:
                print(f"Skipped {code_column} for top {top_n}: no compounds available to plot")
                continue
            output_path = output_dir / f"top_{top_n}_{code_column}_frequency.pdf"
            plot_frequency_table(freq_df=freq_df, top_n=top_n, output_path=output_path)
            print(f"Wrote {output_path}")

    for method in METHOD_ORDER:
        print_method_table(method, ranked_tables[method], code_columns)

    hits_path = output_dir / "hits.csv"
    export_hits_csv(ranked_tables=ranked_tables, code_columns=code_columns, output_path=hits_path)
    print(f"Wrote {hits_path}")

    summary_path = output_dir / "summary.csv"
    export_summary_csv(ranked_tables=ranked_tables, code_columns=code_columns, output_path=summary_path)
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
