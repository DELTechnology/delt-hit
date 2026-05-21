"""Tests for normalized z-score enrichment."""

import numpy as np
import pytest

from delt_hit.analyse.enrichment import calculate_fold_enrichment, calculate_zscore


def _zscore_expected(count: float, total_reads: int, diversity: int) -> float:
    p_i = 1.0 / diversity
    p_obs = count / total_reads
    return (p_obs - p_i) / np.sqrt(p_i * (1.0 - p_i))


def test_zscore_canonical_value():
    z, _, _ = calculate_zscore(np.array([100.0]), total_reads=1000, diversity=100)
    assert pytest.approx(z[0], abs=1e-3) == _zscore_expected(100, 1000, 100)
    assert pytest.approx(z[0], abs=1e-3) == 0.9045


def test_zscore_zero_when_count_equals_expected():
    z, _, _ = calculate_zscore(np.array([10.0]), total_reads=1000, diversity=100)
    assert pytest.approx(z[0], abs=1e-10) == 0.0


def test_zscore_negative_when_below_expected():
    z, _, _ = calculate_zscore(np.array([5.0]), total_reads=1000, diversity=100)
    assert z[0] < 0.0


def test_zscore_positive_when_above_expected():
    z, _, _ = calculate_zscore(np.array([100.0]), total_reads=1000, diversity=100)
    assert z[0] > 0.0


def test_ci_ordering_positive_counts():
    counts = np.array([1.0, 5.0, 10.0, 50.0, 100.0, 500.0])
    z, lo, hi = calculate_zscore(counts, total_reads=1000, diversity=100)
    assert np.all(lo < z)
    assert np.all(z < hi)


def test_ci_ordering_depleted_counts():
    counts = np.array([0.0, 1.0, 5.0])
    z, lo, hi = calculate_zscore(counts, total_reads=1000, diversity=100)
    assert np.all(lo <= z + 1e-12)
    assert np.all(z <= hi + 1e-12)


def test_ci_width_shrinks_with_more_reads():
    _, lo_lo, hi_lo = calculate_zscore(
        np.array([10.0]), total_reads=1000, diversity=100
    )
    _, lo_hi, hi_hi = calculate_zscore(
        np.array([100.0]), total_reads=10000, diversity=100
    )
    assert (hi_hi[0] - lo_hi[0]) < (hi_lo[0] - lo_lo[0])


def test_zscore_count_zero():
    z, _, _ = calculate_zscore(np.array([0.0]), total_reads=1000, diversity=100)
    assert pytest.approx(z[0], abs=1e-6) == _zscore_expected(0, 1000, 100)
    assert z[0] < 0.0


def test_zscore_count_equals_total_reads():
    z, lo, _ = calculate_zscore(np.array([1000.0]), total_reads=1000, diversity=100)
    assert pytest.approx(z[0], abs=1e-6) == _zscore_expected(1000, 1000, 100)
    assert lo[0] < z[0]


def test_zscore_diversity_one_returns_nan():
    z, lo, hi = calculate_zscore(np.array([100.0]), total_reads=1000, diversity=1)
    assert np.isnan(z[0])
    assert np.isnan(lo[0])
    assert np.isnan(hi[0])


def test_zscore_array_shape_preserved():
    counts = np.array([0.0, 5.0, 10.0, 50.0, 100.0])
    z, lo, hi = calculate_zscore(counts, total_reads=1000, diversity=100)
    assert z.shape == counts.shape
    assert lo.shape == counts.shape
    assert hi.shape == counts.shape


def test_zscore_matches_scalar_formula_elementwise():
    counts = np.array([0.0, 5.0, 10.0, 20.0, 100.0])
    z, _, _ = calculate_zscore(counts, total_reads=1000, diversity=100)
    for i, count in enumerate(counts):
        assert pytest.approx(z[i], abs=1e-8) == _zscore_expected(count, 1000, 100)


def test_zscore_raises_on_zero_total_reads():
    with pytest.raises(ValueError, match="total_reads"):
        calculate_zscore(np.array([10.0]), total_reads=0, diversity=100)


def test_zscore_raises_on_zero_diversity():
    with pytest.raises(ValueError, match="diversity"):
        calculate_zscore(np.array([10.0]), total_reads=1000, diversity=0)


def test_fold_enrichment_canonical():
    fold = calculate_fold_enrichment(np.array([100.0]), total_reads=1000, diversity=100)
    assert pytest.approx(fold[0], abs=1e-9) == 10.0


def test_fold_enrichment_at_expected_is_one():
    fold = calculate_fold_enrichment(np.array([10.0]), total_reads=1000, diversity=100)
    assert pytest.approx(fold[0], abs=1e-9) == 1.0


def test_fold_enrichment_zero_count():
    fold = calculate_fold_enrichment(np.array([0.0]), total_reads=1000, diversity=100)
    assert fold[0] == 0.0
