from pathlib import Path
import sys

import pandas as pd
from loguru import logger


base_dir = Path(__file__).resolve().parent
lane = sys.argv[1] if len(sys.argv) > 1 else "lane-1"
published_dir = base_dir / "published"
selections_dir = base_dir / lane / "selections"

published_paths = sorted(published_dir.glob("selection_*_.txt"))

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

    save_path = base_dir / "comparison" / lane / f"{selection_name}.csv"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(save_path, index=True)

    if counts.identical.all():
        logger.info(f"All counts identical for selection {selection_name}")
    else:
        logger.warning(f"Not all counts identical for selection {selection_name}")
