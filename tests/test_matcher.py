"""
Unit tests for matcher module.
"""

import numpy as np
import pytest

from fingerprint.matcher import (
    Thresholds,
    calibrate_thresholds,
    classify_query,
    pair_distance,
    seed_lot_purity_test,
    verify_hybrid,
)


class TestPairDistance:
    """Tests for pair_distance function."""

    def test_identical_homozygotes(self):
        """Test identical homozygous calls."""
        q = np.array([0, 0, 2, 2], dtype=np.int8)
        r = np.array([0, 0, 2, 2], dtype=np.int8)
        dist, n = pair_distance(q, r)
        assert dist == 0.0
        assert n == 4

    def test_opposite_homozygotes(self):
        """Test opposite homozygous calls."""
        q = np.array([0, 0, 2, 2], dtype=np.int8)
        r = np.array([2, 2, 0, 0], dtype=np.int8)
        dist, n = pair_distance(q, r)
        assert dist == 1.0
        assert n == 4

    def test_het_vs_homozygote(self):
        """Test heterozygote vs homozygote."""
        q = np.array([1, 1, 1, 1], dtype=np.int8)
        r = np.array([0, 0, 2, 2], dtype=np.int8)
        dist, n = pair_distance(q, r)
        # Each het-hom comparison = 0.5
        assert dist == 0.5
        assert n == 4

    def test_with_missing(self):
        """Test with missing data."""
        q = np.array([0, 2, -1, 1], dtype=np.int8)
        r = np.array([0, 2, 0, 1], dtype=np.int8)
        dist, n = pair_distance(q, r)
        # 3 loci compared: 0-0=0, 2-2=0, 1-1=0
        assert dist == 0.0
        assert n == 3

    def test_all_missing(self):
        """Test all missing."""
        q = np.array([-1, -1], dtype=np.int8)
        r = np.array([-1, -1], dtype=np.int8)
        dist, n = pair_distance(q, r)
        assert np.isnan(dist)
        assert n == 0


class TestCalibrateThresholds:
    """Tests for calibrate_thresholds function."""

    def test_basic_calibration(self):
        """Test threshold calibration with replicate groups."""
        # 3 varieties, 2 reps each, 10 SNPs
        # Var 0: reps 0,1 - identical
        # Var 1: reps 2,3 - identical
        # Var 2: reps 4,5 - identical
        gt = np.zeros((10, 6), dtype=np.int8)
        gt[:, 0:2] = 0  # Var 0
        gt[:, 2:4] = 2  # Var 1
        gt[:, 4:6] = 1  # Var 2

        ref_ids = np.array(['V0', 'V0', 'V1', 'V1', 'V2', 'V2'])
        replicate_groups = {'V0': [0, 1], 'V1': [2, 3], 'V2': [4, 5]}

        thresholds = calibrate_thresholds(
            ref_gt=gt,
            ref_ids=ref_ids,
            replicate_groups=replicate_groups,
            match_percentile=99.9,
        )

        # Replicate distances should be 0
        assert thresholds.replicate_mean == 0.0
        assert thresholds.replicate_std == 0.0

        # Nearest neighbor distances: each sample's nearest is 0.5 (to V2 or V0/V1)
        # 6 samples, all have nearest distance 0.5
        assert np.isclose(thresholds.nearest_mean, 0.5)

        # Match threshold at 99.9th percentile of replicate (0)
        assert thresholds.match_threshold == 0.0

        # Related threshold = min(3*match, nearest/2)
        assert thresholds.related_threshold == 0.0

    def test_with_noise(self):
        """Test with some genotyping errors."""
        gt = np.zeros((100, 4), dtype=np.int8)
        gt[:, 0:2] = 0  # Var 0 reps
        gt[:, 2:4] = 2  # Var 1 reps

        # Add some errors to replicates
        gt[0, 1] = 1  # Error in rep 1
        gt[1, 1] = 1  # Error in rep 1

        ref_ids = np.array(['V0', 'V0', 'V1', 'V1'])
        replicate_groups = {'V0': [0, 1], 'V1': [2, 3]}

        thresholds = calibrate_thresholds(
            ref_gt=gt,
            ref_ids=ref_ids,
            replicate_groups=replicate_groups,
        )

        # Replicate distance > 0 due to errors
        assert thresholds.replicate_mean > 0
        assert thresholds.match_threshold > 0


class TestClassifyQuery:
    """Tests for classify_query function."""

    def test_perfect_match(self):
        """Test perfect match case."""
        gt = np.zeros((10, 3), dtype=np.int8)
        gt[:, 0] = 0
        gt[:, 1] = 2
        gt[:, 2] = 1

        ref_ids = np.array(['V0', 'V1', 'V2'])
        query = gt[:, 0].copy()  # Exact match to V0

        thresholds = Thresholds(
            match_threshold=0.02,
            related_threshold=0.15,
            margin_to_second=0.05,
            min_loci_fraction=0.9,
            replicate_mean=0.0,
            replicate_std=0.0,
            nearest_mean=0.5,
            nearest_std=0.1,
            auc=1.0,
            eer=0.0,
        )

        result = classify_query(query, gt, ref_ids, thresholds)

        assert result.decision == "MATCH"
        assert result.best_match == "V0"
        assert result.best_distance == 0.0
        assert result.n_loci_used == 10

    def test_no_match(self):
        """Test no match case."""
        gt = np.zeros((10, 3), dtype=np.int8)
        gt[:, 0] = 0
        gt[:, 1] = 2
        gt[:, 2] = 1

        ref_ids = np.array(['V0', 'V1', 'V2'])
        # Query differs from all - use a pattern not matching any variety
        # V0=all 0, V1=all 2, V2=all 1. Use alternating pattern.
        query = np.array([0, 2, 0, 2, 0, 2, 0, 2, 0, 2], dtype=np.int8)

        thresholds = Thresholds(
            match_threshold=0.02,
            related_threshold=0.15,
            margin_to_second=0.05,
            min_loci_fraction=0.9,
            replicate_mean=0.0,
            replicate_std=0.0,
            nearest_mean=0.5,
            nearest_std=0.1,
            auc=1.0,
            eer=0.0,
        )

        result = classify_query(query, gt, ref_ids, thresholds)

        assert result.decision == "NO MATCH"

    def test_insufficient_data(self):
        """Test insufficient data case."""
        gt = np.zeros((100, 2), dtype=np.int8)
        gt[:, 0] = 0
        gt[:, 1] = 2

        ref_ids = np.array(['V0', 'V1'])
        # Query only has 50 calls out of 100
        query = np.full(100, -1, dtype=np.int8)
        query[:50] = 0

        thresholds = Thresholds(
            match_threshold=0.02,
            related_threshold=0.15,
            margin_to_second=0.05,
            min_loci_fraction=0.9,  # Need 90%
            replicate_mean=0.0,
            replicate_std=0.0,
            nearest_mean=0.5,
            nearest_std=0.1,
            auc=1.0,
            eer=0.0,
        )

        result = classify_query(query, gt, ref_ids, thresholds)

        assert result.decision == "INSUFFICIENT DATA"


class TestVerifyHybrid:
    """Tests for verify_hybrid function."""

    def test_true_f1(self):
        """Test true F1 hybrid."""
        # Parents: P1 = all 0, P2 = all 2
        # Expected F1: all 1 (het)
        p1 = np.array([0] * 10, dtype=np.int8)
        p2 = np.array([2] * 10, dtype=np.int8)
        h = np.array([1] * 10, dtype=np.int8)

        result = verify_hybrid(h, p1, p2)

        assert result['concordance'] == 1.0
        assert result['selfed_like'] == 0.0
        assert result['loci_compared'] == 10

    def test_selfed_hybrid(self):
        """Test selfed (homozygous) hybrid."""
        p1 = np.array([0] * 10, dtype=np.int8)
        p2 = np.array([2] * 10, dtype=np.int8)
        # Selfed: some loci become homozygous
        h = np.array([0, 0, 1, 1, 2, 2, 1, 1, 0, 2], dtype=np.int8)

        result = verify_hybrid(h, p1, p2)

        # Concordance: 4/10 match expected het
        assert np.isclose(result['concordance'], 0.4)
        # Selfed-like: 6/10 are homozygous where het expected
        assert np.isclose(result['selfed_like'], 0.6)

    def test_missing_data(self):
        """Test with missing data."""
        p1 = np.array([0, 0, 2, 2], dtype=np.int8)
        p2 = np.array([2, 2, 0, 0], dtype=np.int8)
        h = np.array([1, -1, 1, -1], dtype=np.int8)  # 2 missing

        result = verify_hybrid(h, p1, p2)

        assert result['loci_compared'] == 2
        assert result['concordance'] == 1.0


class TestSeedLotPurity:
    """Tests for seed_lot_purity_test function."""

    def test_zero_offtypes(self):
        """Test zero off-types detected."""
        result = seed_lot_purity_test(n_seeds=100, n_offtypes=0, confidence=0.95)

        assert result['observed_rate'] == 0.0
        assert result['ci_lower'] == 0.0
        assert result['ci_upper'] > 0.0
        assert result['n_for_5pct'] == 59  # ln(0.05)/ln(0.95) ≈ 59

    def test_all_offtypes(self):
        """Test all off-types."""
        result = seed_lot_purity_test(n_seeds=50, n_offtypes=50, confidence=0.95)

        assert result['observed_rate'] == 1.0
        assert result['ci_lower'] < 1.0
        assert result['ci_upper'] == 1.0

    def test_intermediate(self):
        """Test intermediate case."""
        result = seed_lot_purity_test(n_seeds=100, n_offtypes=5, confidence=0.95)

        assert result['observed_rate'] == 0.05
        assert result['ci_lower'] < 0.05
        assert result['ci_upper'] > 0.05


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
