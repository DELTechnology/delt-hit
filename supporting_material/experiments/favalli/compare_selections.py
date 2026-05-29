import pandas as pd
from loguru import logger
from pathlib import Path

base_dir = Path('/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/favalli')
legacy_path = base_dir / "published" / "1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt"
lane = "lane-1-fasta"
selections_dir = base_dir / lane / "selections"
selection_names = sorted(path.stem for path in selections_dir.glob("*") if path.is_dir())
if not selection_names:
    raise SystemExit(
        f"No selection directories found under {selections_dir}. "
        "Run the DELT-Hit workflow for this prefix before generating comparisons."
    )

legacy = pd.read_csv(legacy_path)
comparison_dir = base_dir / "comparison" / lane
comparison_dir.mkdir(parents=True, exist_ok=True)
report_rows: list[dict[str, object]] = []

for selection_name in selection_names:
    delt = pd.read_csv(selections_dir / selection_name / "counts.txt", sep="\t")
    delt = delt.rename(columns={"count": "delt"})
    delt.loc[:, ["code_0", "code_1"]] = delt.loc[:, ["code_0", "code_1"]] + 1

    if selection_name not in legacy.columns:
        logger.warning(f"Selection {selection_name} not found in legacy data, skipping")
        continue

    legacy_selection = legacy[["CodeA", "CodeB", selection_name]].rename(
        columns={"CodeA": "code_0", "CodeB": "code_1", selection_name: "legacy"}
    )

    counts = pd.concat(
        [
            legacy_selection.set_index(["code_0", "code_1"]),
            delt.set_index(["code_0", "code_1"]),
        ],
        axis=1,
        join="outer",
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
