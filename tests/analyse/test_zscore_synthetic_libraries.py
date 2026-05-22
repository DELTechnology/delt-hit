"""Synthetic library tests verifying correct diversity handling in z-score.

These tests establish that:
  1. diversity = theoretical library size (product of per-cycle BB counts), NOT
     the number of observed compounds in the count file.
  2. count_observed_compounds() and infer_diversity_from_cycle_bbs() return
     different values under partial sequencing coverage.
  3. z-scores are numerically correct against the Faver formula using the
     theoretical diversity.
  4. The library_size fallback (when not configured) uses the per-cycle BB
     product, not the observed compound count.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from delt_hit.analyse.enrichment import calculate_zscore
from delt_hit.cli.analyse.api import Analyse
from delt_hit.cli.analyse.zscore import (
    count_observed_compounds,
    infer_diversity_from_cycle_bbs,
    infer_diversity_from_library_config,
    resolve_library_size_from_experiment,
    zscore_analysis,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_full_2bb_data(n_bb1: int, n_bb2: int, count: float = 10.0) -> pd.DataFrame:
    """All n_bb1 × n_bb2 compounds with uniform count, single selection."""
    rows = [
        {"code_1": f"A{i:03d}", "code_2": f"B{j:03d}", "count": count, "name": "sel_1"}
        for i in range(n_bb1)
        for j in range(n_bb2)
    ]
    return pd.DataFrame(rows)


def _make_partial_2bb_data(
    n_bb1: int, n_bb2: int, n_drop: int, count: float = 10.0
) -> pd.DataFrame:
    """n_bb1 × n_bb2 compounds with n_drop removed, all BB identities still represented.

    Anchors are kept: (A*000, B000) for all BB1 values, and (A000, B*) for all BB2
    values. Drops come only from the non-anchor set, so every BB identity remains
    in the retained rows and infer_diversity_from_cycle_bbs returns n_bb1 × n_bb2.
    """
    all_cmpds = [(f"A{i:03d}", f"B{j:03d}") for i in range(n_bb1) for j in range(n_bb2)]
    anchors = (
        {(f"A{i:03d}", "B000") for i in range(n_bb1)}
        | {("A000", f"B{j:03d}") for j in range(n_bb2)}
    )
    droppable = [c for c in all_cmpds if c not in anchors]
    assert n_drop <= len(droppable), (
        f"Cannot drop {n_drop} compounds from {n_bb1}×{n_bb2} library "
        f"while keeping all BB identities (max droppable: {len(droppable)})"
    )
    dropped = set(droppable[:n_drop])
    kept = [c for c in all_cmpds if c not in dropped]
    return pd.DataFrame({
        "code_1": [c[0] for c in kept],
        "code_2": [c[1] for c in kept],
        "count": count,
        "name": "sel_1",
    })


def _make_full_3bb_data(n_bb1: int, n_bb2: int, n_bb3: int, count: float = 5.0) -> pd.DataFrame:
    rows = [
        {"code_1": f"A{i}", "code_2": f"B{j}", "code_3": f"C{k}", "count": count, "name": "sel_1"}
        for i in range(n_bb1)
        for j in range(n_bb2)
        for k in range(n_bb3)
    ]
    return pd.DataFrame(rows)


def _make_delthit_library_config(cycle_sizes: list[int]) -> dict:
    """Create a minimal DELT-Hit demultiplex config.yaml-style dictionary."""
    building_blocks = [f"B{i}" for i in range(len(cycle_sizes))]
    whitelists = {
        name: [
            {
                "index": idx,
                "codon": f"{name}_{idx}",
                "smiles": "C",
                "educt": "scaffold_1",
                "reaction": "R1",
                "product": f"{name}_product_{idx}",
            }
            for idx in range(size)
        ]
        for name, size in zip(building_blocks, cycle_sizes)
    }
    return {
        "experiment": {"name": "synthetic", "save_dir": "campaign"},
        "selections": {},
        "library": {
            "building_blocks": building_blocks,
            "products": [],
            "educts": ["scaffold_1"],
            "reactions": ["R1"],
            "bb_edges": [],
            "other_edges": [],
        },
        "catalog": {
            "compounds": {"scaffold_1": {"smiles": "C"}},
            "reactions": {"R1": {"smirks": "[C:1]>>[C:1]"}},
        },
        "structure": [
            {"name": name, "type": "building_block"}
            for name in building_blocks
        ],
        "whitelists": whitelists,
    }


def _write_yaml(path: Path, data: dict) -> None:
    import yaml

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False))


def _write_delthit_counts(
    path: Path,
    cycle_sizes: list[int],
    rows_to_keep: int,
) -> None:
    """Write a sparse DELT-Hit counts.txt with 0-indexed code columns."""
    import itertools

    code_columns = [f"code_{idx}" for idx in range(len(cycle_sizes))]
    combinations = list(itertools.product(*[range(size) for size in cycle_sizes]))
    kept = combinations[:rows_to_keep]
    frame = pd.DataFrame(kept, columns=code_columns)
    frame["count"] = [100] + [10] * (len(frame) - 1)
    frame["id"] = [
        "_".join(str(value) for value in row)
        for row in kept
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def _write_delthit_zscore_fixture(
    tmp_path: Path,
    cycle_sizes: list[int],
) -> tuple[Path, Path]:
    """Write analysis.yaml + campaign/config.yaml + sparse selection counts."""
    campaign = tmp_path / "campaign"
    config_path = campaign / "config.yaml"
    _write_yaml(config_path, _make_delthit_library_config(cycle_sizes))

    # Keep fewer rows than the full product so observed compound count cannot
    # accidentally equal the designed diversity.
    diversity = int(np.prod(cycle_sizes))
    rows_to_keep = max(2, diversity - 1)
    protein_counts = campaign / "selections" / "protein_1" / "counts.txt"
    control_counts = campaign / "selections" / "control_1" / "counts.txt"
    _write_delthit_counts(protein_counts, cycle_sizes, rows_to_keep)
    _write_delthit_counts(control_counts, cycle_sizes, rows_to_keep)

    analysis_path = tmp_path / "analysis.yaml"
    _write_yaml(
        analysis_path,
        {
            "experiments": [
                {
                    "name": "synthetic",
                    "save_dir": str(tmp_path / "analysis_output"),
                    "selections": [
                        {
                            "name": "protein_1",
                            "counts_path": str(protein_counts),
                            "group": "protein",
                        },
                        {
                            "name": "control_1",
                            "counts_path": str(control_counts),
                            "group": "no_protein",
                        },
                    ],
                }
            ]
        },
    )
    return analysis_path, config_path


# ---------------------------------------------------------------------------
# infer_diversity_from_cycle_bbs — unit tests
# ---------------------------------------------------------------------------

class TestInferDiversity:
    def test_2bb_full_coverage(self):
        data = _make_full_2bb_data(n_bb1=10, n_bb2=20)
        assert infer_diversity_from_cycle_bbs(data, ["code_1", "code_2"]) == 200

    def test_2bb_partial_coverage_all_bbs_present(self):
        # 5×8 = 40 theoretical; drop 10 compounds but keep all BB identities
        data = _make_partial_2bb_data(n_bb1=5, n_bb2=8, n_drop=10)
        assert infer_diversity_from_cycle_bbs(data, ["code_1", "code_2"]) == 40

    def test_3bb_full_coverage(self):
        data = _make_full_3bb_data(n_bb1=5, n_bb2=6, n_bb3=4)
        assert infer_diversity_from_cycle_bbs(data, ["code_1", "code_2", "code_3"]) == 120

    def test_single_cycle_equals_unique_bb_count(self):
        data = pd.DataFrame({"code_1": ["A", "B", "C", "D"], "count": [1.0]*4, "name": "sel_1"})
        assert infer_diversity_from_cycle_bbs(data, ["code_1"]) == 4

    def test_returns_int(self):
        data = _make_full_2bb_data(n_bb1=3, n_bb2=4)
        result = infer_diversity_from_cycle_bbs(data, ["code_1", "code_2"])
        assert isinstance(result, int)


# ---------------------------------------------------------------------------
# count_observed_compounds — unit tests
# ---------------------------------------------------------------------------

class TestCountObservedCompounds:
    def test_full_coverage_equals_theoretical(self):
        data = _make_full_2bb_data(n_bb1=5, n_bb2=6)
        assert count_observed_compounds(data, ["code_1", "code_2"]) == 30

    def test_partial_coverage_less_than_theoretical(self):
        # 5×8 = 40 theoretical, drop 10 → 30 observed
        data = _make_partial_2bb_data(n_bb1=5, n_bb2=8, n_drop=10)
        assert count_observed_compounds(data, ["code_1", "code_2"]) == 30

    def test_returns_int(self):
        data = _make_full_2bb_data(n_bb1=2, n_bb2=3)
        assert isinstance(count_observed_compounds(data, ["code_1", "code_2"]), int)


# ---------------------------------------------------------------------------
# Diversity vs observed count — the core conceptual distinction
# ---------------------------------------------------------------------------

class TestDiversityVsObservedCount:
    def test_inferred_diversity_exceeds_observed_count_under_partial_coverage(self):
        """The key invariant: theoretical diversity > observed compound count when
        sequencing coverage is incomplete."""
        data = _make_partial_2bb_data(n_bb1=10, n_bb2=20, n_drop=50)  # 150 of 200
        code_cols = ["code_1", "code_2"]

        inferred = infer_diversity_from_cycle_bbs(data, code_cols)
        observed = count_observed_compounds(data, code_cols)

        assert inferred == 200       # all 10×20 BBs still appear
        assert observed == 150       # only 150 compounds in the file
        assert inferred > observed

    def test_equal_when_coverage_is_complete(self):
        data = _make_full_2bb_data(n_bb1=6, n_bb2=7)
        code_cols = ["code_1", "code_2"]
        assert infer_diversity_from_cycle_bbs(data, code_cols) == count_observed_compounds(data, code_cols)

    def test_zscore_expected_fraction_differs_by_denominator(self, tmp_path):
        """Using library_size=theoretical (12) vs library_size=observed (10) gives
        different expected fractions and therefore different z-scores."""
        bb1 = ["A", "B", "C"]
        bb2 = ["X", "Y", "Z", "W"]
        all_compounds = [(b1, b2) for b1 in bb1 for b2 in bb2]  # 12 total
        kept = all_compounds[:10]  # 10 observed

        data = pd.DataFrame({
            "code_1": [c[0] for c in kept],
            "code_2": [c[1] for c in kept],
            "count": [20.0] * 10,
            "name": "sel_1",
        })
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=12,
                        save_dir=tmp_path / "theoretical")
        zscore_analysis(data=data, samples=samples, library_size=10,
                        save_dir=tmp_path / "observed")

        theoretical = pd.read_csv(tmp_path / "theoretical" / "sel_1.csv")
        observed_div = pd.read_csv(tmp_path / "observed" / "sel_1.csv")

        assert pytest.approx(theoretical["expected_fraction"].iloc[0]) == 1 / 12
        assert pytest.approx(observed_div["expected_fraction"].iloc[0]) == 1 / 10

        # z-scores must differ because p_i differs
        assert not np.allclose(
            theoretical["zscore"].values,
            observed_div["zscore"].values,
        )

    def test_wrong_diversity_inflates_expected_fraction(self, tmp_path):
        """Using too-small diversity (n_observed=80 instead of theoretical=100) as
        denominator inflates p_i and lowers z-scores for enriched compounds."""
        # 10×10 = 100 theoretical; drop 20 non-anchor compounds → 80 observed
        data = _make_partial_2bb_data(n_bb1=10, n_bb2=10, n_drop=20)
        # Make the anchor compound A000-B000 a strong hit
        data.loc[(data["code_1"] == "A000") & (data["code_2"] == "B000"), "count"] = 500.0
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        # Correct: theoretical library_size = 100 (>= 80 observed → passes validation)
        zscore_analysis(data=data, samples=samples, library_size=100,
                        save_dir=tmp_path / "correct")
        # Wrong: library_size = 80 = n_observed; just barely passes validation (80 >= 80)
        zscore_analysis(data=data, samples=samples, library_size=80,
                        save_dir=tmp_path / "wrong")

        correct = pd.read_csv(tmp_path / "correct" / "sel_1.csv")
        wrong = pd.read_csv(tmp_path / "wrong" / "sel_1.csv")

        hit_correct = correct[correct["zscore"] == correct["zscore"].max()].iloc[0]
        hit_wrong = wrong[wrong["zscore"] == wrong["zscore"].max()].iloc[0]

        # With a smaller denominator (wrong), p_i is larger → hit z-score is lower
        assert hit_correct["zscore"] > hit_wrong["zscore"]


# ---------------------------------------------------------------------------
# Numerical correctness against the Faver formula
# ---------------------------------------------------------------------------

class TestNumericalCorrectness:
    """Verify z-score values exactly against the closed-form Faver formula."""

    def test_2bb_uniform_counts_all_zero_zscore(self, tmp_path):
        """When every compound has the expected count, all z-scores = 0."""
        n_bb1, n_bb2 = 5, 4          # diversity = 20
        total_reads = 200.0          # expected_count = 200/20 = 10
        data = _make_full_2bb_data(n_bb1=n_bb1, n_bb2=n_bb2, count=10.0)
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=20,
                        save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")

        assert pytest.approx(result["zscore"].values, abs=1e-10) == np.zeros(len(result))

    def test_2bb_hit_compound_exact_zscore(self, tmp_path):
        """One enriched compound — verify z-score against hand-computed formula."""
        # Library: 5 BB1 × 5 BB2 = 25 compounds
        # Hit: count=100, all others: count=5
        n_bb1, n_bb2 = 5, 5
        diversity = n_bb1 * n_bb2  # 25
        hit_count = 100.0
        background_count = 5.0
        n_non_hit = diversity - 1   # 24

        total_reads = hit_count + n_non_hit * background_count  # 220

        rows = []
        for i in range(n_bb1):
            for j in range(n_bb2):
                c = hit_count if (i == 0 and j == 0) else background_count
                rows.append({"code_1": f"A{i}", "code_2": f"B{j}", "count": c, "name": "sel_1"})
        data = pd.DataFrame(rows)
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=diversity, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")
        hit_row = result[(result["code_1"] == "A0") & (result["code_2"] == "B0")].iloc[0]

        p_i = 1.0 / diversity
        p_obs = hit_count / total_reads
        expected_z = (p_obs - p_i) / np.sqrt(p_i * (1 - p_i))

        assert pytest.approx(hit_row["zscore"], abs=1e-8) == expected_z
        assert pytest.approx(hit_row["expected_fraction"], abs=1e-12) == p_i
        assert pytest.approx(hit_row["expected_count"], abs=1e-10) == total_reads / diversity
        assert pytest.approx(hit_row["fold_enrichment"], abs=1e-10) == hit_count / (total_reads / diversity)

    def test_3bb_hit_compound_exact_zscore(self, tmp_path):
        """Three-cycle library — verify z-score for a hit compound."""
        n_bb1, n_bb2, n_bb3 = 4, 5, 3   # diversity = 60
        diversity = n_bb1 * n_bb2 * n_bb3
        hit_count = 200.0
        background_count = 2.0
        n_non_hit = diversity - 1

        total_reads = hit_count + n_non_hit * background_count

        rows = []
        for i in range(n_bb1):
            for j in range(n_bb2):
                for k in range(n_bb3):
                    c = hit_count if (i == 0 and j == 0 and k == 0) else background_count
                    rows.append({"code_1": f"A{i}", "code_2": f"B{j}", "code_3": f"C{k}",
                                 "count": c, "name": "sel_1"})
        data = pd.DataFrame(rows)
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=diversity, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")
        hit_row = result[
            (result["code_1"] == "A0") & (result["code_2"] == "B0") & (result["code_3"] == "C0")
        ].iloc[0]

        p_i = 1.0 / diversity
        p_obs = hit_count / total_reads
        expected_z = (p_obs - p_i) / np.sqrt(p_i * (1 - p_i))

        assert pytest.approx(hit_row["zscore"], abs=1e-8) == expected_z

    def test_depleted_compound_has_negative_zscore(self, tmp_path):
        data = _make_full_2bb_data(n_bb1=5, n_bb2=4, count=10.0)
        # Replace one compound with a very low count (depleted)
        data.loc[(data["code_1"] == "A000") & (data["code_2"] == "B000"), "count"] = 1.0
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=20, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")
        depleted = result[(result["code_1"] == "A000") & (result["code_2"] == "B000")].iloc[0]

        assert depleted["zscore"] < 0.0

    def test_zero_count_gives_large_negative_zscore(self, tmp_path):
        """A zero-filled compound has the most negative z-score."""
        # protein_1 sees A and B; control_1 sees A, B, C — so C is zero-filled in protein_1
        data = pd.DataFrame([
            {"code_1": "A", "name": "protein_1", "count": 50.0},
            {"code_1": "B", "name": "protein_1", "count": 50.0},
            {"code_1": "A", "name": "control_1", "count": 30.0},
            {"code_1": "B", "name": "control_1", "count": 30.0},
            {"code_1": "C", "name": "control_1", "count": 40.0},
        ])
        samples = pd.DataFrame([
            {"name": "protein_1", "group": "protein"},
            {"name": "control_1", "group": "no_protein"},
        ])

        zscore_analysis(data=data, samples=samples, library_size=3, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "protein_1.csv")

        zero_row = result[result["code_1"] == "C"].iloc[0]
        assert zero_row["count"] == 0.0

        p_i = 1.0 / 3
        p_obs = 0.0 / 100.0
        expected_z = (p_obs - p_i) / np.sqrt(p_i * (1 - p_i))
        assert pytest.approx(zero_row["zscore"], abs=1e-8) == expected_z
        assert zero_row["zscore"] < 0.0
        assert zero_row["zscore"] == result["zscore"].min()


# ---------------------------------------------------------------------------
# library_size fallback behaviour
# ---------------------------------------------------------------------------

class TestLibrarySizeFallback:
    def test_fallback_uses_bb_product_not_observed_count(self, tmp_path):
        """When library_size is None, the fallback is the per-cycle BB product."""
        bb1 = ["A", "B", "C"]       # 3 BBs
        bb2 = ["X", "Y", "Z", "W"]  # 4 BBs → theoretical = 12
        all_cmpds = [(b1, b2) for b1 in bb1 for b2 in bb2]

        data = pd.DataFrame({
            "code_1": [c[0] for c in all_cmpds],
            "code_2": [c[1] for c in all_cmpds],
            "count": [10.0] * 12,
            "name": "sel_1",
        })
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=None, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")

        assert pytest.approx(result["expected_fraction"].iloc[0]) == 1 / 12

    def test_fallback_partial_coverage_still_uses_bb_product(self, tmp_path):
        """Partial coverage: fallback is 3×4=12 (from BBs), not 8 (from observed)."""
        bb1 = ["A", "B", "C"]
        bb2 = ["X", "Y", "Z", "W"]
        all_cmpds = [(b1, b2) for b1 in bb1 for b2 in bb2]

        # Drop 4 compounds but keep all BB identities (A,B,C and X,Y,Z,W all appear)
        to_drop = {("A", "X"), ("B", "Y"), ("C", "Z"), ("A", "W")}
        kept = [c for c in all_cmpds if c not in to_drop]
        assert len(kept) == 8

        data = pd.DataFrame({
            "code_1": [c[0] for c in kept],
            "code_2": [c[1] for c in kept],
            "count": [10.0] * 8,
            "name": "sel_1",
        })
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        zscore_analysis(data=data, samples=samples, library_size=None, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")

        # All 3+4 BBs are present → inferred diversity = 3×4 = 12
        assert pytest.approx(result["expected_fraction"].iloc[0]) == 1 / 12

    def test_configured_library_size_beats_fallback(self, tmp_path):
        """Explicitly configured library_size takes precedence over inference."""
        data = _make_full_2bb_data(n_bb1=3, n_bb2=4)  # 12 observed = 12 inferred
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        # Configure a larger theoretical library (some BBs not in this run)
        zscore_analysis(data=data, samples=samples, library_size=50, save_dir=tmp_path)
        result = pd.read_csv(tmp_path / "sel_1.csv")

        assert pytest.approx(result["expected_fraction"].iloc[0]) == 1 / 50

    def test_library_size_smaller_than_observed_raises(self, tmp_path):
        data = _make_full_2bb_data(n_bb1=5, n_bb2=4)  # 20 observed
        samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

        with pytest.raises(ValueError, match="smaller than the number of observed"):
            zscore_analysis(data=data, samples=samples, library_size=10, save_dir=tmp_path)


class TestDelthitConfigDiversity:
    @pytest.mark.parametrize(
        ("cycle_sizes", "expected"),
        [
            ([5], 5),
            ([3, 7], 21),
            ([4, 5, 6], 120),
            ([2, 5, 3, 4], 120),
        ],
    )
    def test_infers_exact_diversity_from_config_whitelists(
        self,
        cycle_sizes,
        expected,
    ):
        config = _make_delthit_library_config(cycle_sizes)
        code_columns = [f"code_{idx}" for idx in range(len(cycle_sizes))]

        assert infer_diversity_from_library_config(config, code_columns) == expected

    def test_config_diversity_rejects_code_column_mismatch(self):
        config = _make_delthit_library_config([2, 3, 4])

        with pytest.raises(ValueError, match="does not match count-file code columns"):
            infer_diversity_from_library_config(config, ["code_0", "code_1"])

    @pytest.mark.parametrize("cycle_sizes", [[2, 5], [3, 4, 2], [2, 5, 3, 4]])
    def test_analysis_api_auto_discovers_delthit_config_yaml(
        self,
        tmp_path,
        cycle_sizes,
    ):
        analysis_path, _ = _write_delthit_zscore_fixture(tmp_path, cycle_sizes)
        expected_diversity = int(np.prod(cycle_sizes))

        Analyse().enrichment(
            config_path=analysis_path,
            name="synthetic",
            method="zscore",
        )

        output = pd.read_csv(
            tmp_path
            / "analysis_output"
            / "synthetic"
            / "zscore"
            / "protein_1.csv"
        )
        assert pytest.approx(output["expected_fraction"].iloc[0]) == (
            1 / expected_diversity
        )

    def test_explicit_library_config_path_is_resolved_relative_to_analysis_yaml(
        self,
        tmp_path,
    ):
        analysis_path, _ = _write_delthit_zscore_fixture(tmp_path, [3, 4, 5])
        exp = {
            "library_config_path": "campaign/config.yaml",
            "selections": [],
        }

        resolved = resolve_library_size_from_experiment(
            exp=exp,
            analysis_config_path=analysis_path,
            code_columns=["code_0", "code_1", "code_2"],
        )

        assert resolved == 60


# ---------------------------------------------------------------------------
# Multi-replicate group statistics
# ---------------------------------------------------------------------------

class TestGroupStats:
    def test_replicate_zscores_averaged_correctly(self, tmp_path):
        """Group z-score in stats.csv is the mean of per-replicate z-scores."""
        # Two protein replicates with the same compound, different counts
        data = pd.DataFrame([
            {"code_1": "A", "code_2": "B", "count": 80.0, "name": "protein_1"},
            {"code_1": "A", "code_2": "C", "count": 20.0, "name": "protein_1"},
            {"code_1": "A", "code_2": "B", "count": 60.0, "name": "protein_2"},
            {"code_1": "A", "code_2": "C", "count": 40.0, "name": "protein_2"},
        ])
        samples = pd.DataFrame([
            {"name": "protein_1", "group": "protein"},
            {"name": "protein_2", "group": "protein"},
        ])

        zscore_analysis(data=data, samples=samples, library_size=4, save_dir=tmp_path)

        r1 = pd.read_csv(tmp_path / "protein_1.csv")
        r2 = pd.read_csv(tmp_path / "protein_2.csv")
        stats = pd.read_csv(tmp_path / "stats.csv")

        hit_z_r1 = r1[r1["code_2"] == "B"]["zscore"].iloc[0]
        hit_z_r2 = r2[r2["code_2"] == "B"]["zscore"].iloc[0]
        expected_mean = (hit_z_r1 + hit_z_r2) / 2.0

        stats_hit = stats[stats["code_2"] == "B"].iloc[0]
        assert pytest.approx(stats_hit["protein"], abs=1e-8) == expected_mean

    def test_enrichment_column_is_protein_minus_no_protein(self, tmp_path):
        data = pd.DataFrame([
            {"code_1": "A", "count": 80.0, "name": "protein_1"},
            {"code_1": "B", "count": 20.0, "name": "protein_1"},
            {"code_1": "A", "count": 20.0, "name": "control_1"},
            {"code_1": "B", "count": 80.0, "name": "control_1"},
        ])
        samples = pd.DataFrame([
            {"name": "protein_1", "group": "protein"},
            {"name": "control_1", "group": "no_protein"},
        ])

        zscore_analysis(data=data, samples=samples, library_size=2, save_dir=tmp_path)
        stats = pd.read_csv(tmp_path / "stats.csv")

        for _, row in stats.iterrows():
            assert pytest.approx(row["enrichment"], abs=1e-10) == row["protein"] - row["no_protein"]
