"""Enrichment scoring methods for DEL count data."""

from __future__ import annotations

import numpy as np


def calculate_zscore(
    counts: np.ndarray,
    total_reads: int,
    diversity: int,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the normalized z-score enrichment (Faver et al. 2019).

    The normalized z-score is:

    ``z_n = (p_observed - p_i) / sqrt(p_i * (1 - p_i))``

    where ``p_observed = count / total_reads`` and ``p_i = 1 / diversity``.
    Confidence intervals use an Agresti-Coull interval on ``p_observed`` and
    are transformed back into z-score units.

    Args:
        counts: Observed read counts, one value per compound.
        total_reads: Total decoded reads for the selection. Must be positive.
        diversity: Compound-level library diversity. Must be at least one.
        alpha: Significance level for confidence intervals.

    Returns:
        Tuple of ``(zscore, ci_lower, ci_upper)`` arrays.
    """
    if total_reads <= 0:
        raise ValueError(f"total_reads must be > 0, got {total_reads}")
    if diversity < 1:
        raise ValueError(f"diversity must be >= 1, got {diversity}")

    counts = np.asarray(counts, dtype=float)

    p_i = 1.0 / diversity
    denominator = np.sqrt(p_i * (1.0 - p_i))
    if np.isclose(denominator, 0.0):
        nan_array = np.full_like(counts, np.nan)
        return nan_array, nan_array, nan_array

    p_observed = counts / total_reads
    zscore = (p_observed - p_i) / denominator

    z_alpha = _z_alpha_from_significance(alpha)
    n_prime = total_reads + z_alpha**2
    p_adj = (counts + z_alpha**2 / 2.0) / n_prime
    se_adj = np.sqrt(p_adj * (1.0 - p_adj) / n_prime)

    ci_p_lower = p_adj - z_alpha * se_adj
    ci_p_upper = p_adj + z_alpha * se_adj
    ci_lower = (ci_p_lower - p_i) / denominator
    ci_upper = (ci_p_upper - p_i) / denominator

    return zscore, ci_lower, ci_upper


def calculate_fold_enrichment(
    counts: np.ndarray,
    total_reads: int,
    diversity: int,
) -> np.ndarray:
    """Compute fold enrichment relative to uniform library expectation."""
    if total_reads <= 0:
        raise ValueError(f"total_reads must be > 0, got {total_reads}")
    if diversity < 1:
        raise ValueError(f"diversity must be >= 1, got {diversity}")

    counts = np.asarray(counts, dtype=float)
    expected_count = total_reads / diversity
    return counts / expected_count


def _z_alpha_from_significance(alpha: float) -> float:
    """Return the two-tailed normal critical value for ``alpha``."""
    from scipy.stats import norm

    return float(norm.ppf(1.0 - alpha / 2.0))
