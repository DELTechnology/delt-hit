from pathlib import Path
import sys

import pandas as pd
from loguru import logger


base_dir = Path(__file__).resolve().parent
lane = sys.argv[1] if len(sys.argv) > 1 else "lane-1"
published_dir = base_dir / "published"
selections_dir = base_dir / lane / "selections"

published_paths = sorted(published_dir.glob("selection_*_.txt"))
selection_dirs = [path for path in selections_dir.glob("*") if path.is_dir()]
if not selection_dirs:
    raise SystemExit(
        f"No selection directories found under {selections_dir}. "
        "Run the DELT-Hit workflow for this lane before generating comparisons."
    )

comparison_dir = base_dir / "comparison" / lane
comparison_dir.mkdir(parents=True, exist_ok=True)
report_rows: list[dict[str, object]] = []

for published_path in published_paths:
    selection_name = published_path.stem.removeprefix("selection_").removesuffix("_")
    delt_path = selections_dir / selection_name / "counts.txt"

    if not delt_path.exists():
        logger.warning(f"Selection {selection_name} not found in DELT-Hit output, skipping")
        continue

    delt = pd.read_csv(delt_path, sep="\t")
    delt = delt.rename(columns={"count": "delt"})
    delt.loc[:, ["code_0", "code_1", "code_2"]] = delt.loc[:, ["code_0", "code_1", "code_2"]] + 1

    legacy = pd.read_csv(published_path, sep="\t").rename(
        columns={"Count": "legacy", "Code1": "code_0", "Code2": "code_1", "Code3": "code_2"}
    )

    counts = pd.concat(
        [
            legacy.set_index(["code_0", "code_1", "code_2"]),
            delt.set_index(["code_0", "code_1", "code_2"]),
        ],
        axis=1,
    )
    counts = counts[counts.legacy != 0]
    counts = counts.convert_dtypes()
    counts["identical"] = counts.legacy == counts.delt

    counts.to_csv(comparison_dir / f"{selection_name}.csv", index=True)

    mismatch_count = int((~counts["identical"].fillna(False)).sum())
    report_rows.append(
        {
            "selection": selection_name,
            "legacy_observed_compounds": int(counts["legacy"].notna().sum()),
            "delt_observed_compounds": int(counts["delt"].notna().sum()),
            "legacy_total_counts": int(counts["legacy"].fillna(0).sum()),
            "delt_total_counts": int(counts["delt"].fillna(0).sum()),
            "mismatch_count": mismatch_count,
            "all_identical": mismatch_count == 0,
        }
    )

    if counts.identical.all():
        logger.info(f"All counts identical for selection {selection_name}")
    else:
        logger.warning(f"Not all counts identical for selection {selection_name}")

report = pd.DataFrame(
    report_rows,
    columns=[
        "selection",
        "legacy_observed_compounds",
        "delt_observed_compounds",
        "legacy_total_counts",
        "delt_total_counts",
        "mismatch_count",
        "all_identical",
    ],
)
report.to_csv(comparison_dir / "report.csv", index=False)
