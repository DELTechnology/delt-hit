import json
import math
from pathlib import Path

import numpy as np
from rich.console import Console
from rich.table import Table


def print_report(output_dir: Path, save_path: Path) -> None:
    """Pretty, aligned cutadapt pipeline report using Rich."""
    report_files = sorted(output_dir.glob("*.cutadapt.json"))
    reports = {
        p.stem.replace(".cutadapt", ""): json.load(p.open("r"))
        for p in report_files
    }

    if not reports:
        console = Console()
        console.print("[red]No '*.cutadapt.json' files found.[/red]")
        save_path.write_text("No '*.cutadapt.json' files found.")
        return

    # Build row dicts
    rows = []
    for region_id, rep in sorted(reports.items(), key=lambda kv: kv[0]):
        c_input = int(rep["read_counts"]["input"])
        c_output = int(rep["read_counts"]["output"])
        prop_out = np.nan if c_input == 0 else c_output / c_input
        rows.append(
            dict(
                region_id=region_id,
                total_in=c_input,
                with_adapters=c_output,
                discarded=c_input - c_output,
                p_with=prop_out,
                p_discarded=(np.nan if math.isnan(prop_out) else 1 - prop_out),
            )
        )

    first, last = rows[0], rows[-1]
    overall_in = first["total_in"]
    overall_out = last["with_adapters"]
    overall_p_with = (overall_out / overall_in) if overall_in else float("nan")
    overall_discarded = overall_in - overall_out
    overall_p_discarded = (1 - overall_p_with) if not math.isnan(overall_p_with) else float("nan")

    def fmt_int(x: int) -> str:
        return f"{x:,}"

    def fmt_pct(x: float) -> str:
        return "N/A" if (x is None or math.isnan(x)) else f"{x:.2%}"

    console = Console(record=True)  # record=True lets us export to text later

    table = Table(title="Cutadapt Pipeline Report")
    table.add_column("Region", justify="left", style="cyan", no_wrap=True)
    table.add_column("Input", justify="right")
    table.add_column("With adapters", justify="right")
    table.add_column("Discarded", justify="right")
    table.add_column("% with", justify="right")
    table.add_column("% discarded", justify="right")

    for r in rows:
        p_dis = r["p_discarded"]
        pct_dis = fmt_pct(p_dis)
        pct_dis = f"[red]{pct_dis}[/red]" if not math.isnan(p_dis) and p_dis > 0.10 else pct_dis
        table.add_row(
            str(r["region_id"]),
            fmt_int(r["total_in"]),
            fmt_int(r["with_adapters"]),
            fmt_int(r["discarded"]),
            fmt_pct(r["p_with"]),
            pct_dis,
        )

    console.print(table)

    # Overall summary
    console.print("\n[bold]Overall[/bold]:")
    console.print(
        f"  With adapters : {fmt_int(overall_out)} ({fmt_pct(overall_p_with)})"
    )
    console.print(
        f"  Discarded     : {fmt_int(overall_discarded)} ({fmt_pct(overall_p_discarded)})"
    )

    # Save report as plain text (without ANSI codes)
    save_path.write_text(console.export_text())
