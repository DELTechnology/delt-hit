#!/usr/bin/env python3
"""Build DELT-Hit inputs for the synthetic two-cycle FASTQ fixture."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_DIR = ROOT / "experiments"
LONG_CSV = EXPERIMENT_DIR / "synthetic_2cycle.csv"
WORKBOOK = EXPERIMENT_DIR / "synthetic_2cycle.xlsx"
FASTQ_IN = EXPERIMENT_DIR / "synthetic_2cycle.fastq"
FASTQ_GZ = EXPERIMENT_DIR / "synthetic_2cycle.fastq.gz"

LIBRARY_TAG = "ACGTACGTAC"
UMI_PLUS_CLOSING = "NNNNNNNNTTGGAACC"
BB1_COUNT = 100
BB2_COUNT = 10
TAG_LENGTH = 8
QUALITY_CHAR = "I"

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


def build_sheets() -> dict[str, pd.DataFrame]:
    bb1_tags = make_tags(BB1_COUNT, offset=0)
    bb2_tags = make_tags(BB2_COUNT, offset=10_000)

    sheets: dict[str, pd.DataFrame] = {
        "experiment": pd.DataFrame(
            [
                ("name", "synthetic_2cycle"),
                ("fastq_path", str(FASTQ_GZ)),
                ("save_dir", str(EXPERIMENT_DIR)),
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


def write_long_csv(sheets: dict[str, pd.DataFrame]) -> None:
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
    pd.DataFrame(rows).to_csv(LONG_CSV, index=False)


def write_workbook(sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(WORKBOOK, engine="openpyxl") as writer:
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False)


def write_gzipped_fastq() -> None:
    bb1_tags = make_tags(BB1_COUNT, offset=0)
    bb2_tags = make_tags(BB2_COUNT, offset=10_000)
    with FASTQ_IN.open("w") as handle:
        read_idx = 1
        for bb1_idx, bb1_tag in enumerate(bb1_tags, start=1):
            for bb2_idx, bb2_tag in enumerate(bb2_tags, start=1):
                umi = int_to_dna(20_000 + read_idx, TAG_LENGTH)
                sequence = f"{LIBRARY_TAG}{bb1_tag}{bb2_tag}{umi}TTGGAACC"
                quality = QUALITY_CHAR * len(sequence)
                handle.write(
                    f"@synthetic_combo_{read_idx:04d}_BB1_{bb1_idx:03d}_BB2_{bb2_idx:03d}\n"
                    f"{sequence}\n"
                    "+\n"
                    f"{quality}\n"
                )
                read_idx += 1

    with FASTQ_IN.open("rb") as src, gzip.open(FASTQ_GZ, "wb") as dst:
        shutil.copyfileobj(src, dst)


def main() -> None:
    EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)
    sheets = build_sheets()
    write_long_csv(sheets)
    write_workbook(sheets)
    write_gzipped_fastq()
    print(f"Wrote long-form CSV source to {LONG_CSV}")
    print(f"Wrote DELT-Hit workbook to {WORKBOOK}")
    print(f"Wrote FASTQ to {FASTQ_IN}")
    print(f"Wrote gzipped FASTQ copy to {FASTQ_GZ}")


if __name__ == "__main__":
    main()
