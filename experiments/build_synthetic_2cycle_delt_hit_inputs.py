#!/usr/bin/env python3
"""Build DELT-Hit inputs for the synthetic two-cycle FASTQ fixture."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "experiments"
EXPERIMENT_NAME = "synthetic_2cycle"
LIBRARY_TAG = "ACGTACGTAC"
UMI_PLUS_CLOSING = "NNNNNNNNTTGGAACC"
BB1_COUNT = 100
BB2_COUNT = 10
TAG_LENGTH = 8
QUALITY_CHAR = "I"
DEFAULT_NUM_READS_PER_COMPOUND = 1

AMIDE_SMIRKS = (
    "[CX3:1](=[O:2])[OX2;H1]."
    "[N;$([N;H1,H2]),$([N;H3]);"
    "!$([N;H1,H2][S]);"
    "!$([NX3][CX3](=[OX1])[#6]);"
    "!$([NX3][CX3](=[OX1])[OX2]):4]"
    ">>[CX3:1](=[O:2])[N:4]"
)


def int_to_dna(value: int, length: int = TAG_LENGTH) -> str:
    alphabet = "ACGT"
    bases: list[str] = []
    for _ in range(length):
        bases.append(alphabet[value % 4])
        value //= 4
    return "".join(reversed(bases))


def make_tags(count: int, offset: int) -> list[str]:
    return [int_to_dna(offset + idx) for idx in range(count)]


def alkyl_amine(index: int) -> str:
    """Return a simple primary alkyl amine SMILES."""
    motifs = [
        "NC",
        "NCC",
        "NCCC",
        "NCCCC",
        "NCCO",
        "NCCN",
        "NC1CC1",
        "NCCF",
        "NCCCl",
        "NCCBr",
    ]
    return motifs[index % len(motifs)]


def build_sheets(*, output_dir: Path) -> dict[str, pd.DataFrame]:
    bb1_tags = make_tags(BB1_COUNT, offset=0)
    bb2_tags = make_tags(BB2_COUNT, offset=10_000)
    fastq_gz = output_dir / f"{EXPERIMENT_NAME}.fastq.gz"

    sheets: dict[str, pd.DataFrame] = {
        "experiment": pd.DataFrame(
            [
                ("name", EXPERIMENT_NAME),
                ("fastq_path", str(fastq_gz)),
                ("save_dir", str(output_dir)),
                ("num_cores", 1),
            ],
            columns=["variable", "value"],
        ),
        "selection": pd.DataFrame(
            [
                {
                    "name": "synthetic_selection",
                    "operator": "Codex",
                    "date": "2026-04-21",
                    "group": "synthetic",
                    "target": "synthetic_target",
                    "S0": LIBRARY_TAG,
                    "S1": UMI_PLUS_CLOSING,
                }
            ]
        ),
        "structure": pd.DataFrame(
            [
                {"name": "S0", "type": "selection", "max_error_rate": 0, "indels": False},
                {"name": "B0", "type": "building_block", "max_error_rate": 0, "indels": False},
                {"name": "B1", "type": "building_block", "max_error_rate": 0, "indels": False},
                {"name": "S1", "type": "selection", "max_error_rate": 0, "indels": False},
            ]
        ),
        "constant": pd.DataFrame([{"name": "C0", "codon": "TTGGAACC"}]),
        "B0": pd.DataFrame(
            [
                {
                    "codon": tag,
                    "smiles": alkyl_amine(i),
                    "educt": "scaffold_1",
                    "reaction": "AMIDE1",
                    "product": "product_1",
                }
                for i, tag in enumerate(bb1_tags)
            ]
        ),
        "B1": pd.DataFrame(
            [
                {
                    "codon": tag,
                    "smiles": alkyl_amine(i + 100),
                    "educt": "product_1",
                    "reaction": "AMIDE2",
                    "product": "product_2",
                }
                for i, tag in enumerate(bb2_tags)
            ]
        ),
        "compounds": pd.DataFrame(
            [{"name": "scaffold_1", "smiles": "O=C(O)CCC(=O)O"}]
        ),
        "reactions": pd.DataFrame(
            [
                {"name": "AMIDE1", "smirks": AMIDE_SMIRKS},
                {"name": "AMIDE2", "smirks": AMIDE_SMIRKS},
            ]
        ),
        "reaction_graph": pd.DataFrame(columns=["educt_1", "educt_2", "reaction", "product"]),
    }
    return sheets


def write_long_csv(sheets: dict[str, pd.DataFrame], *, path: Path) -> None:
    rows = []
    for sheet_name, frame in sheets.items():
        for row_index, row in frame.reset_index(drop=True).iterrows():
            for column_name, value in row.items():
                if pd.notna(value):
                    rows.append(
                        {
                            "sheet": sheet_name,
                            "row": row_index,
                            "column": column_name,
                            "value": value,
                        }
                    )
    pd.DataFrame(rows).to_csv(path, index=False)


def write_workbook(sheets: dict[str, pd.DataFrame], *, path: Path) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False)


def write_gzipped_fastq(
    *,
    fastq_path: Path,
    fastq_gz_path: Path,
    num_reads_per_compound: int,
) -> None:
    bb1_tags = make_tags(BB1_COUNT, offset=0)
    bb2_tags = make_tags(BB2_COUNT, offset=10_000)
    with fastq_path.open("w") as handle:
        read_idx = 1
        for bb1_idx, bb1_tag in enumerate(bb1_tags, start=1):
            for bb2_idx, bb2_tag in enumerate(bb2_tags, start=1):
                for replicate_idx in range(1, num_reads_per_compound + 1):
                    umi = int_to_dna(20_000 + read_idx, TAG_LENGTH)
                    sequence = f"{LIBRARY_TAG}{bb1_tag}{bb2_tag}{umi}TTGGAACC"
                    quality = QUALITY_CHAR * len(sequence)
                    handle.write(
                        f"@synthetic_combo_{read_idx:04d}_read_{replicate_idx:03d}_"
                        f"BB1_{bb1_idx:03d}_BB2_{bb2_idx:03d}\n"
                        f"{sequence}\n"
                        "+\n"
                        f"{quality}\n"
                    )
                    read_idx += 1

    with fastq_path.open("rb") as src, gzip.open(fastq_gz_path, "wb") as dst:
        shutil.copyfileobj(src, dst)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the workbook, FASTQ, and long-form CSV are written.",
    )
    parser.add_argument(
        "--num-reads-per-compound",
        type=int,
        default=DEFAULT_NUM_READS_PER_COMPOUND,
        help="FASTQ reads to emit for each BB1/BB2 combination.",
    )
    args = parser.parse_args()
    if args.num_reads_per_compound < 1:
        parser.error("num_reads_per_compound must be at least 1")
    return args


def main(*, output_dir: Path = DEFAULT_OUTPUT_DIR, num_reads_per_compound: int = DEFAULT_NUM_READS_PER_COMPOUND) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    long_csv = output_dir / f"{EXPERIMENT_NAME}.csv"
    workbook = output_dir / f"{EXPERIMENT_NAME}.xlsx"
    fastq_path = output_dir / f"{EXPERIMENT_NAME}.fastq"
    fastq_gz_path = output_dir / f"{EXPERIMENT_NAME}.fastq.gz"

    sheets = build_sheets(output_dir=output_dir)
    write_long_csv(sheets, path=long_csv)
    write_workbook(sheets, path=workbook)
    write_gzipped_fastq(
        fastq_path=fastq_path,
        fastq_gz_path=fastq_gz_path,
        num_reads_per_compound=num_reads_per_compound,
    )
    print(f"Wrote long-form CSV source to {long_csv}")
    print(f"Wrote DELT-Hit workbook to {workbook}")
    print(f"Wrote FASTQ to {fastq_path}")
    print(f"Wrote gzipped FASTQ copy to {fastq_gz_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(
        output_dir=cli_args.output_dir,
        num_reads_per_compound=cli_args.num_reads_per_compound,
    )
