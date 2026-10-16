"""
Unit tests for robustness module.
"""

import numpy as np
import pytest

from fingerprint.robustness import (
    compute_analytical_power,
    evaluate_panel,
    simulate_assay_dropout,
    simulate_genotyping_error,
    simulate_missing_data,
    simulate_seed_mixture,
)


class TestSimulateGenotypingError:
    """Tests for genotyping error simulation."""

    def test_no_error(self):
        """Test with zero error rate."""
        gt = np.array([[0, 2, 1], [0, 2, 1], [2, 0, 1]], dtype=np.int8).T
        rng = np.random.default_rng(42)

        gt_err = simulate_genotyping_error(gt, 0.0, rng)
        assert np.array_equal(gt_err, gt)

    def test_with_error(self):
        """Test with non-zero error rate."""
        gt = np.zeros((1000, 10), dtype=np.int8)
        gt[:, :5] = 0
        gt[:, 5:] = 2
        rng = np.random.default_rng(42)

        gt_err = simulate_genotyping_error(gt, 0.1, rng)

        # About 10% of called genotypes should be flipped
        called = gt >= 0
        n_called = called.sum()
        n_changed = (gt_err[called] != gt[called]).sum()

        # Should be roughly 10% (binomial variation)
        assert 0.05 < n_changed / n_called < 0.15

    def test_error_preserves_missing(self):
        """Test that missing values stay missing."""
        gt = np.array([[0, -1, 2], [-1, 2, 0]], dtype=np.int8).T
        rng = np.random.default_rng(42)

        gt_err = simulate_genotyping_error(gt, 1.0, rng)  # 100% error

        # Missing should stay missing
        assert gt_err[1, 0] == -1
        assert gt_err[0, 1] == -1


class TestSimulateMissingData:
    """Tests for missing data simulation."""

    def test_no_missing(self):
        """Test with zero missing rate."""
        gt = np.array([[0, 2, 1], [0, 2, 1]], dtype=np.int8).T
        rng = np.random.default_rng(42)

        gt_miss = simulate_missing_data(gt, 0.0, rng)
        assert np.array_equal(gt_miss, gt)

    def test_with_missing(self):
        """Test with non-zero missing rate."""
        gt = np.ones((1000, 10), dtype=np.int8) * 2
        rng = np.random.default_rng(42)

        gt_miss = simulate_missing_data(gt, 0.2, rng)

        # About 20% should be missing
        n_missing = (gt_miss == -1).sum()
        n_total = gt_miss.size
        assert 0.15 < n_missing / n_total < 0.25


class TestSimulateSeedMixture:
    """Tests for seed mixture simulation."""

    def test_no_mixture(self):
        """Test with zero mixture proportion."""
        gt = np.array([[0, 2], [0, 2], [1, 1]], dtype=np.int8).T
        rng = np.random.default_rng(42)

        gt_mix = simulate_seed_mixture(gt, 0.0, rng)
        assert np.array_equal(gt_mix, gt)

    def test_with_mixture(self):
        """Test with non-zero mixture."""
        # Two distinct varieties
        gt = np.zeros((100, 2), dtype=np.int8)
        gt[:, 0] = 0
        gt[:, 1] = 2
        rng = np.random.default_rng(42)

        gt_mix = simulate_seed_mixture(gt, 0.5, rng)

        # About 50% of samples should be mixed
        # Mixed samples show het at differentiating loci
        n_het = (gt_mix == 1).sum()
        n_total = (gt_mix >= 0).sum()

        # In mixture, differentiating loci become het
        # Only sample 0 gets mixed (sample 1 is the contaminant source)
        # At 50% mixture rate, ~50 samples mixed
        # Each mixed sample has 100 het loci
        assert n_het > 0


class TestSimulateAssayDropout:
    """Tests for assay dropout simulation."""

    def test_no_dropout(self):
        """Test with zero dropout."""
        panel = np.array([0, 1, 2, 3, 4])
        backup = np.array([5, 6, 7])
        rng = np.random.default_rng(42)

        new_panel, mask = simulate_assay_dropout(panel, 0.0, rng, backup)
        assert np.array_equal(new_panel, panel)
        assert mask.all()

    def test_with_dropout(self):
        """Test with dropout."""
        panel = np.arange(10)
        backup = np.arange(10, 15)
        rng = np.random.default_rng(42)

        new_panel, mask = simulate_assay_dropout(panel, 0.3, rng, backup)

        # 30% dropout = 3 dropped
        assert mask.sum() == 7  # 7 kept
        assert len(new_panel) == 10  # 7 kept + 3 backup

    def test_dropout_no_backup(self):
        """Test dropout without backup."""
        panel = np.arange(10)
        rng = np.random.default_rng(42)

        new_panel, mask = simulate_assay_dropout(panel, 0.3, rng, None)

        assert mask.sum() == 7
        assert len(new_panel) == 7


class TestEvaluatePanel:
    """Tests for panel evaluation."""

    def test_perfect_separation(self):
        """Test perfectly separated varieties."""
        # 3 varieties, 10 SNPs, all different
        gt = np.array([
            [0, 0, 2, 2, 1, 1, 0, 2, 1, 0],
            [2, 2, 0, 0, 1, 1, 2, 0, 1, 2],
            [1, 1, 1, 1, 0, 2, 1, 1, 0, 1],
        ], dtype=np.int8).T  # (10, 3)

        panel = np.arange(10)
        n_corr, n_tot, acc, mean_sep, min_sep = evaluate_panel(gt, panel, k=1)

        assert n_tot == 3
        assert n_corr == 3
        assert acc == 1.0

    def test_identical_varieties(self):
        """Test identical varieties (clones)."""
        gt = np.array([
            [0, 2, 1],
            [0, 2, 1],
            [1, 0, 2],
        ], dtype=np.int8).T  # (3, 3) - first two identical

        panel = np.arange(3)
        n_corr, n_tot, acc, mean_sep, min_sep = evaluate_panel(gt, panel, k=1)

        # First two varieties cannot be distinguished
        assert n_corr == 1  # Only third variety correctly identified
        assert acc == 1/3


class TestAnalyticalPower:
    """Tests for analytical power computation."""

    def test_basic(self):
        """Test basic power calculation."""
        result = compute_analytical_power(
            panel_size=50,
            k=3,
            error_rate=0.01,
            missing_rate=0.05,
        )

        assert 'prob_k_or_more' in result
        assert 0 <= result['prob_k_or_more'] <= 1
        assert result['expected_called'] == 50 * 0.95
        assert result['expected_correct'] == 50 * 0.95 * 0.99

    def test_high_error(self):
        """Test with high error rate and small panel."""
        # Use small panel size and high k so high error actually reduces power
        result = compute_analytical_power(
            panel_size=10,
            k=9,
            error_rate=0.10,
            missing_rate=0.10,
        )

        # High error should reduce power with small panel and high k
        assert result['prob_k_or_more'] < 0.5

    def test_large_panel(self):
        """Test with large panel."""
        result = compute_analytical_power(
            panel_size=384,
            k=3,
            error_rate=0.02,
            missing_rate=0.10,
        )

        # Large panel should have high power even with errors
        assert result['prob_k_or_more'] > 0.9


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
