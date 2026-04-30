from pathlib import Path

import pandas as pd
from loguru import logger


base_dir = Path(__file__).resolve().parent
legacy_path = base_dir / "published" / "1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt"
selection_names = [path.stem for path in (base_dir / "lane-1" / "selections").glob("*_1") if path.is_dir()]

for selection_name in sorted(selection_names):
    delt = pd.read_csv(base_dir / "lane-1" / "selections" / selection_name / "counts.txt", sep="\t")
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
    )
    counts = counts[counts.legacy != 0]
    counts = counts.convert_dtypes()
    counts["identical"] = counts.legacy == counts.delt

    save_path = base_dir / "comparison" / f"{selection_name}.csv"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(save_path, index=True)

    if counts.identical.all():
        logger.info(f"All counts identical for selection {selection_name}")
    else:
        logger.warning(f"Not all counts identical for selection {selection_name}")
