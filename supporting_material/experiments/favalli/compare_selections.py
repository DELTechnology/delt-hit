from pathlib import Path
import sys

import pandas as pd
from loguru import logger


base_dir = Path(__file__).resolve().parent
legacy_path = base_dir / "published" / "1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt"
lane = sys.argv[1] if len(sys.argv) > 1 else "lane-1"
selection_names = [path.stem for path in (base_dir / lane / "selections").glob("*") if path.is_dir()]
selection_names = sorted(selection_names)

for selection_name in sorted(selection_names):
    delt = pd.read_csv(base_dir / lane / "selections" / selection_name / "counts.txt", sep="\t")
    delt = delt.rename(columns={"count": "delt"})
    delt.loc[:, ["code_0", "code_1"]] = delt.loc[:, ["code_0", "code_1"]] + 1

    legacy = pd.read_csv(legacy_path)
    if selection_name not in legacy.columns:
        logger.warning(f"Selection {selection_name} not found in legacy data, skipping")
        continue

    legacy = legacy[["CodeA", "CodeB", selection_name]].rename(
        columns={"CodeA": "code_0", "CodeB": "code_1", selection_name: "legacy"}
    )

    counts = pd.concat(
        [
            legacy.set_index(["code_0", "code_1"]),
            delt.set_index(["code_0", "code_1"]),
        ],
        axis=1,
        join="outer",
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
