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
METHOD_ORDER = ["diff", "ratio", "edgeR"]
MAX_CODES_PER_PLOT = 10
PALETTE = {
    "diff": "#4C78A8",
    "ratio": "#F58518",
    "edgeR": "#54A24B",
}
QUERY_COMPOUNDS = [(159, 474), (357, 474), (68, 474)]
TOP_HITS_EXPORT_N = 10
console = Console()


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_data_dir = script_dir / "lane-1-fasta" / "analysis" / "protein_vs_no_protein"

    parser = argparse.ArgumentParser(
        description="Plot code_0 and code_1 frequencies among top-ranked compounds."
    )
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N, help="Number of top compounds to select per method.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=default_data_dir,
        help="Analysis directory containing counts/stats.csv and edgeR/enrichment_stats.csv, for example lane-1/analysis/protein_vs_no_protein.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where plots and CSV summaries will be written. Defaults to data-dir.",
    )
    return parser.parse_args()


def validate_top_n(top_n: int) -> None:
    if top_n <= 0:
        raise ValueError("--top-n must be a positive integer.")


def build_top_tables(counts_df: pd.DataFrame, edger_df: pd.DataFrame, top_n: int) -> dict[str, pd.DataFrame]:
    ranked_tables = build_ranked_tables(counts_df=counts_df, edger_df=edger_df)
    return {method: table.head(top_n).copy() for method, table in ranked_tables.items()}


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
        all_codes.update(table[code_column].astype(int).tolist())

    ordered_codes = sorted(all_codes)

    for method in METHOD_ORDER:
        table = top_tables[method]
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
    plt.close()


def print_method_table(method: str, table: pd.DataFrame) -> None:
    rich_table = Table(title=f"Top {min(TOP_HITS_EXPORT_N, len(table))} enriched compounds by {method}")
    rich_table.add_column("Rank", justify="right")
    rich_table.add_column("code_0", justify="right")
    rich_table.add_column("code_1", justify="right")
    rich_table.add_column("Score", justify="right")

    for _, row in table.head(TOP_HITS_EXPORT_N).iterrows():
        rich_table.add_row(
            str(int(row["rank"])),
            str(int(row["code_0"])),
            str(int(row["code_1"])),
            f"{row['score']:.6g}",
        )

    console.print(rich_table)


def print_query_table(ranked_tables: dict[str, pd.DataFrame]) -> None:
    rich_table = Table(title="Queried compound enrichment by method")
    rich_table.add_column("Compound")
    rich_table.add_column("Method")
    rich_table.add_column("Rank", justify="right")
    rich_table.add_column("Score", justify="right")

    for code_0, code_1 in QUERY_COMPOUNDS:
        compound = f"{code_0};{code_1}"
        for method in METHOD_ORDER:
            table = ranked_tables[method]
            match = table.loc[(table["code_0"] == code_0) & (table["code_1"] == code_1)]
            if match.empty:
                rich_table.add_row(compound, method, "-", "not found")
                continue
            row = match.iloc[0]
            rich_table.add_row(compound, method, str(int(row["rank"])), f"{row['score']:.6g}")

    console.print(rich_table)


def build_controls_table(ranked_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for code_0, code_1 in QUERY_COMPOUNDS:
        compound = f"{code_0};{code_1}"
        for method in METHOD_ORDER:
            table = ranked_tables[method]
            match = table.loc[(table["code_0"] == code_0) & (table["code_1"] == code_1)]
            if match.empty:
                rows.append(
                    {"compound": compound, "method": method, "rank": pd.NA, "score": pd.NA}
                )
                continue
            row = match.iloc[0]
            rows.append(
                {
                    "compound": compound,
                    "method": method,
                    "rank": int(row["rank"]),
                    "score": float(row["score"]),
                }
            )
    return pd.DataFrame(rows)


def export_hits_csv(ranked_tables: dict[str, pd.DataFrame], output_path: Path) -> None:
    frames = []
    for method in METHOD_ORDER:
        table = ranked_tables[method].head(TOP_HITS_EXPORT_N).copy()
        frames.append(
            table.loc[:, ["method", "rank", "code_0", "code_1", "score"]]
        )

    pd.concat(frames, ignore_index=True).to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()
    validate_top_n(args.top_n)

    data_dir = args.data_dir.expanduser().resolve()
    counts_path = data_dir / "counts" / "stats.csv"
    edger_path = data_dir / "edgeR" / "enrichment_stats.csv"
    output_dir = data_dir if args.output_dir is None else args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    counts_df = pd.read_csv(counts_path)
    edger_df = pd.read_csv(edger_path)

    ranked_tables = build_ranked_tables(counts_df=counts_df, edger_df=edger_df)
    top_tables = {method: table.head(args.top_n).copy() for method, table in ranked_tables.items()}

    for method in METHOD_ORDER:
        print(f"{method}: selected {len(top_tables[method])} compounds")

    for code_column in ["code_0", "code_1"]:
        freq_df = compute_frequency_table(top_tables=top_tables, code_column=code_column)
        output_path = output_dir / f"top_{args.top_n}_{code_column}_frequency.pdf"
        plot_frequency_table(freq_df=freq_df, top_n=args.top_n, output_path=output_path)
        print(f"Wrote {output_path}")

    print_query_table(ranked_tables)
    for method in METHOD_ORDER:
        print_method_table(method, ranked_tables[method])

    hits_path = output_dir / "hits.csv"
    export_hits_csv(ranked_tables=ranked_tables, output_path=hits_path)
    print(f"Wrote {hits_path}")
    controls_path = output_dir / "controls.csv"
    build_controls_table(ranked_tables=ranked_tables).to_csv(controls_path, index=False)
    print(f"Wrote {controls_path}")


if __name__ == "__main__":
    main()
