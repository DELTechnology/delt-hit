from __future__ import annotations

import argparse
from pathlib import Path
import re

import pandas as pd
from loguru import logger


def _selection_sort_key(selection_name: str) -> tuple[int, int, str]:
    match = re.search(r"(\d+)(?:_(\d+))?$", selection_name)
    if not match:
        return (10**9, 10**9, selection_name)
    first = int(match.group(1))
    second = int(match.group(2) or 0)
    return (first, second, selection_name)


def _detect_mode(base_dir: Path) -> str:
    favalli_file = base_dir / "published" / "1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt"
    if favalli_file.exists():
        return "favalli"

    pure_del_files = list((base_dir / "published").glob("selection_*_*.txt"))
    if pure_del_files:
        return "pure-del"

    raise FileNotFoundError(f"Could not determine comparison mode from {base_dir}")


def _candidate_selection_dirs(base_dir: Path, mode: str) -> list[Path]:
    selection_dirs: list[Path] = []
    if mode == "favalli":
        lane_dirs = [base_dir / "lane-1"]
    else:
        lane_dirs = sorted(base_dir.glob("lane-*.bk")) + sorted(base_dir.glob("lane-*"))

    for lane_dir in lane_dirs:
        selections_dir = lane_dir / "selections"
        if selections_dir.is_dir():
            selection_dirs.extend(sorted(path for path in selections_dir.iterdir() if path.is_dir()))
    return selection_dirs


def _load_delt_counts(counts_path: Path) -> pd.DataFrame:
    delt = pd.read_csv(counts_path, sep="\t")
    delt = delt.rename(columns={"count": "delt"})
    code_cols = [column for column in delt.columns if column.startswith("code_")]
    delt.loc[:, code_cols] = delt.loc[:, code_cols] + 1
    return delt


def _compare_favalli(base_dir: Path, selection_dir: Path) -> pd.DataFrame | None:
    selection_name = selection_dir.name
    legacy_path = base_dir / "published" / "1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt"
    legacy = pd.read_csv(legacy_path)
    if selection_name not in legacy.columns:
        logger.warning(f"Selection {selection_name} not found in legacy data, skipping")
        return None

    delt = _load_delt_counts(selection_dir / "counts.txt")
    legacy = legacy[["CodeA", "CodeB", selection_name]].rename(
        columns={"CodeA": "code_0", "CodeB": "code_1", selection_name: "legacy"}
    )
    return _finalize_counts(delt=delt, legacy=legacy, index_cols=["code_0", "code_1"])


def _compare_pure_del(base_dir: Path, selection_dir: Path) -> pd.DataFrame | None:
    selection_name = selection_dir.name
    selection_suffix = selection_name.removeprefix("MK_")
    lane_match = re.search(r"lane-(\d+)", str(selection_dir))
    if lane_match is None:
        raise ValueError(f"Could not determine lane index from {selection_dir}")
    lane_index = lane_match.group(1)

    legacy_path = base_dir / "published" / f"selection_{selection_suffix}_{lane_index}_.txt"
    if not legacy_path.exists():
        logger.warning(f"Published file not found for selection {selection_name}, skipping")
        return None

    delt = _load_delt_counts(selection_dir / "counts.txt")
    legacy = pd.read_csv(legacy_path, sep="\t").rename(
        columns={"Count": "legacy", "Code1": "code_0", "Code2": "code_1", "Code3": "code_2"}
    )
    return _finalize_counts(delt=delt, legacy=legacy, index_cols=["code_0", "code_1", "code_2"])


def _finalize_counts(delt: pd.DataFrame, legacy: pd.DataFrame, index_cols: list[str]) -> pd.DataFrame:
    delt = delt.set_index(index_cols)
    legacy = legacy.set_index(index_cols)
    counts = pd.concat([legacy, delt], axis=1)
    counts = counts[counts.legacy != 0]
    counts = counts.convert_dtypes()
    counts["identical"] = counts.legacy == counts.delt
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=Path, nargs="?", default=Path(__file__).resolve().parent / "favalli")
    args = parser.parse_args()

    base_dir = args.base_dir.resolve()
    mode = _detect_mode(base_dir)
    compare_dir = base_dir / "comparison"
    compare_dir.mkdir(parents=True, exist_ok=True)

    for selection_dir in sorted(_candidate_selection_dirs(base_dir, mode), key=lambda path: _selection_sort_key(path.name)):
        if mode == "favalli":
            counts = _compare_favalli(base_dir, selection_dir)
            if counts is None:
                continue
        else:
            counts = _compare_pure_del(base_dir, selection_dir)
            if counts is None:
                continue

        save_path = compare_dir / f"{selection_dir.name}.csv"
        counts.to_csv(save_path, index=True)
        if counts.identical.all():
            logger.info(f"All counts identical for selection {selection_dir.name}")
        else:
            logger.warning(f"Not all counts identical for selection {selection_dir.name}")


if __name__ == "__main__":
    main()
