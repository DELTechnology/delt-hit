from collections import defaultdict
import gzip

import pandas as pd

from delt_hit.demultiplex.postprocess import (
    _count_bytes_chunk,
    extract_ids,
    get_counts,
    save_counts,
)
from delt_hit.demultiplex.parser import whitelists_from_excel


def test_building_block_whitelist_indices_are_zero_based():
    whitelists = whitelists_from_excel("templates/library.xlsx")

    first_bb_name = sorted(name for name in whitelists if name.startswith("B"))[0]
    first_indices = [item["index"] for item in whitelists[first_bb_name][:3]]

    assert first_indices == [0, 1, 2]


def test_extract_ids_keeps_building_block_codes_zero_based():
    ids = extract_ids("@read?0-S0.0?1-B0.0?2-B1.3?3-S1.0")

    assert ids["selection_ids"] == (0, 0)
    assert ids["barcodes"] == (0, 3)


def test_bytes_chunk_parser_matches_zero_based_extract_ids():
    counts = _count_bytes_chunk(b"@read?0-S0.0?1-B0.0?2-B1.3?3-S1.0\n")

    assert counts == {(0, 0): {(0, 3): 1}}


def test_save_counts_writes_zero_based_code_columns(tmp_path):
    counts = defaultdict(lambda: defaultdict(int))
    counts[(0, 0)][(0, 3)] = 2

    save_counts(counts, output_dir=tmp_path, as_files=True, sort_by_counts=False)

    table = pd.read_csv(tmp_path / "0-0_counts.txt", sep="\t")
    assert table.to_dict("records") == [{"code_0": 0, "code_1": 3, "count": 2, "id": "0_3"}]


def test_get_counts_and_save_counts_preserve_expected_codon_counts(tmp_path):
    input_path = tmp_path / "reads_with_adapters.gz"
    lines = [
        "@read_0001?0-S0.0?1-B0.0?2-B1.3?3-S1.0\n",
        "@read_0002?0-S0.0?1-B0.0?2-B1.3?3-S1.0\n",
        "@read_0003?0-S0.0?1-B0.2?2-B1.1?3-S1.0\n",
        "@read_0004?0-S0.1?1-B0.4?2-B1.0?3-S1.0\n",
        "@read_0005?0-S0.1?1-B0.4?2-B1.0?3-S1.0\n",
        "@read_0006?0-S0.1?1-B0.4?2-B1.0?3-S1.0\n",
    ]

    with gzip.open(input_path, "wt") as handle:
        handle.writelines(lines)

    counts = get_counts(input_path=input_path, num_reads=len(lines), num_workers=1)

    assert {selection_ids: dict(barcode_counts) for selection_ids, barcode_counts in counts.items()} == {
        (0, 0): {(0, 3): 2, (2, 1): 1},
        (1, 0): {(4, 0): 3},
    }

    output_dir = tmp_path / "selections"
    ids_to_name = {(0, 0): "selection_a", (1, 0): "selection_b"}
    save_counts(counts, output_dir=output_dir, ids_to_name=ids_to_name, as_files=True, sort_by_counts=True)

    selection_a = pd.read_csv(output_dir / "selection_a_counts.txt", sep="\t")
    selection_b = pd.read_csv(output_dir / "selection_b_counts.txt", sep="\t")

    assert selection_a.to_dict("records") == [
        {"code_0": 0, "code_1": 3, "count": 2, "id": "0_3"},
        {"code_0": 2, "code_1": 1, "count": 1, "id": "2_1"},
    ]
    assert selection_b.to_dict("records") == [
        {"code_0": 4, "code_1": 0, "count": 3, "id": "4_0"},
    ]
