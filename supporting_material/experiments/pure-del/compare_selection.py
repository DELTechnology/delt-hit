from pathlib import Path

import pandas as pd
from loguru import logger


base_dir = Path(__file__).resolve().parent

for lane_name in ["lane-1.bk", "lane-2.bk"]:
    selections_dir = base_dir / lane_name / "selections"
    lane_index = lane_name.split(".")[0].split("-")[1]
    selection_names = [path.name for path in selections_dir.glob("MK_*") if path.is_dir()]

    for selection_name in sorted(selection_names, key=lambda name: int(name.removeprefix("MK_"))):
        published_name = f"selection_{selection_name.removeprefix('MK_')}_{lane_index}_.txt"
        published_path = base_dir / "published" / published_name
        if not published_path.exists():
            logger.warning(f"Published file not found for selection {selection_name}, skipping")
            continue

        delt = pd.read_csv(selections_dir / selection_name / "counts.txt", sep="\t")
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

        save_path = base_dir / "comparison" / f"{selection_name}.csv"
        save_path.parent.mkdir(parents=True, exist_ok=True)
        counts.to_csv(save_path, index=True)

        if counts.identical.all():
            logger.info(f"All counts identical for selection {selection_name}")
        else:
            logger.warning(f"Not all counts identical for selection {selection_name}")
