from pathlib import Path

import pytest

from delt_hit.cli.dashboard.api import resolve_selection_name


def test_resolve_selection_name_from_counts_subdirectory() -> None:
    selections = {"AG24_4": {}, "AG24_5": {}}

    assert (
        resolve_selection_name(
            counts_path=Path("campaign/selections/AG24_4/counts.txt"),
            selections=selections,
        )
        == "AG24_4"
    )


def test_resolve_selection_name_from_flat_counts_filename() -> None:
    selections = {"AG24_4": {}, "AG24_5": {}}

    assert (
        resolve_selection_name(
            counts_path=Path("campaign/selections/AG24_4_counts.txt"),
            selections=selections,
        )
        == "AG24_4"
    )


def test_resolve_selection_name_prefers_explicit_override() -> None:
    selections = {"AG24_4": {}, "AG24_5": {}}

    assert (
        resolve_selection_name(
            counts_path=Path("campaign/selections/AG24_4_counts.txt"),
            selections=selections,
            selection_name="AG24_5",
        )
        == "AG24_5"
    )


def test_resolve_selection_name_raises_for_unknown_selection() -> None:
    selections = {"AG24_4": {}, "AG24_5": {}}

    with pytest.raises(ValueError, match="available selections"):
        resolve_selection_name(
            counts_path=Path("campaign/selections/unknown_counts.txt"),
            selections=selections,
        )
