"""
Unit tests for panel selection module.
"""

import numpy as np
import pytest

from fingerprint.selection import (
    compute_pairwise_separation,
    compute_resolution_curve,
    lazy_greedy_multicover,
    select_panel,
)


class TestLazyGreedyMulticover:
    """Tests for lazy greedy multicover algorithm."""

    def test_simple_separation(self):
        """Test simple case with perfect separation."""
        # 4 samples, 3 SNPs that together separate all pairs
        gt = np.array([
            [0, 2, 0, 2],  # SNP 0: separates (0,1), (0,3), (1,2), (2,3)
            [0, 0, 2, 2],  # SNP 1: separates (0,2), (0,3), (1,2), (1,3)
            [0, 2, 2, 0],  # SNP 2: separates (0,1), (0,2), (1,3), (2,3)
        ], dtype=np.int8)  # (3 variants, 4 samples)

        k = 1
        max_snps = 10

        result = lazy_greedy_multicover(
            gt=gt,
            k=k,
            max_snps=max_snps,
            candidate_indices=np.arange(3),
            verbose=False,
        )

        # Should select at least 2 SNPs to separate all pairs
        assert len(result.panel_indices) >= 2
        assert result.unresolved_pairs == 0

    def test_k_greater_than_one(self):
        """Test k=2 requiring two separating SNPs per pair."""
        # 3 samples, 4 SNPs that together give k=2 separation for all pairs
        gt = np.array([
            [0, 2, 2],  # SNP 0: separates (0,1), (0,2)
            [0, 2, 2],  # SNP 1: separates (0,1), (0,2) - redundant
            [2, 0, 2],  # SNP 2: separates (0,1), (1,2)
            [2, 0, 2],  # SNP 3: separates (0,1), (1,2) - redundant
        ], dtype=np.int8)  # (4 variants, 3 samples)

        k = 2
        max_snps = 10

        result = lazy_greedy_multicover(
            gt=gt,
            k=k,
            max_snps=max_snps,
            candidate_indices=np.arange(4),
            verbose=False,
        )

        # Should select at least 2 SNPs per pair (all 4 for k=2)
        assert len(result.panel_indices) >= 4
        assert result.unresolved_pairs == 0

    def test_spacing_constraint(self):
        """Test physical spacing constraint."""
        gt = np.array([
            [0, 2, 0],  # SNP 0
            [0, 2, 0],  # SNP 1
            [0, 0, 2],  # SNP 2
        ], dtype=np.int8)  # (3 variants, 3 samples)
        pos = np.array([100, 101, 1000000])
        chrom = np.array(['chr1', 'chr1', 'chr1'])

        result = lazy_greedy_multicover(
            gt=gt,
            k=1,
            max_snps=10,
            candidate_indices=np.arange(3),
            pos=pos,
            chrom=chrom,
            min_gap_bp=1000,
            verbose=False,
        )

        # Should not select both SNP 0 and SNP 1 (too close)
        selected = set(result.panel_indices)
        if 0 in selected:
            assert 1 not in selected
        if 1 in selected:
            assert 0 not in selected

    def test_chromosome_quota(self):
        """Test per-chromosome quota."""
        gt = np.array([
            [0, 2, 0, 0],  # chr1
            [0, 2, 0, 0],  # chr1
            [0, 2, 0, 0],  # chr1
            [0, 0, 2, 0],  # chr2
        ], dtype=np.int8)  # (4 variants, 4 samples)

        chrom = np.array(['chr1', 'chr1', 'chr1', 'chr2'])

        result = lazy_greedy_multicover(
            gt=gt,
            k=1,
            max_snps=10,
            candidate_indices=np.arange(4),
            chrom=chrom,
            max_chrom_fraction=0.5,  # max 50% per chrom
            verbose=False,
        )

        # Max 2 from chr1 (50% of 4)
        chr1_count = sum(1 for s in result.panel_indices if chrom[s] == 'chr1')
        assert chr1_count <= 2


class TestComputePairwiseSeparation:
    """Tests for pairwise separation computation."""

    def test_basic(self):
        """Test basic separation counting."""
        gt = np.array([
            [0, 2, 0],  # SNP 0
            [0, 0, 2],  # SNP 1
        ], dtype=np.int8)

        panel_indices = np.array([0, 1])
        iu = np.array([0, 0, 1])
        ju = np.array([1, 2, 2])

        sep = compute_pairwise_separation(gt, panel_indices, iu, ju)

        # Pair 0-1: diff at SNP 0 only -> 1
        # Pair 0-2: diff at SNP 1 only -> 1
        # Pair 1-2: diff at both SNP 0 and 1 -> 2
        assert sep[0] == 1
        assert sep[1] == 1
        assert sep[2] == 2


class TestSelectPanel:
    """Integration tests for select_panel (requires full setup)."""

    def test_toy_example(self):
        """Test with toy data."""
        # 5 samples, 10 variants
        np.random.seed(42)
        gt = np.random.randint(0, 3, size=(10, 5), dtype=np.int8)
        # Make some variants informative
        gt[0] = [0, 2, 0, 0, 0]
        gt[1] = [0, 0, 2, 0, 0]
        gt[2] = [0, 0, 0, 2, 0]

        variant_ids = np.array([f'SNP{i}' for i in range(10)])
        chrom = np.array(['chr1'] * 5 + ['chr2'] * 5)
        pos = np.arange(10) * 1000000
        scores = np.array([1.0, 0.9, 0.8, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

        config = {
            "k": 1,
            "core_max_snps": 5,
            "backup_fraction": 0.2,
            "extended_max_snps": 10,
            "min_gap_bp": 100000,
            "max_chrom_fraction": 0.5,
            "ilp_time_limit": 10,
            "min_maf": 0.0,
        }

        core, backup, extended = select_panel(
            gt=gt,
            variant_ids=variant_ids,
            chrom=chrom,
            pos=pos,
            scores=scores,
            config=config,
        )

        assert len(core.panel_indices) > 0
        assert len(core.variant_ids) == len(core.panel_indices)
        assert len(backup.panel_indices) > 0
        assert len(extended.panel_indices) > 0


class TestResolutionCurve:
    """Tests for resolution curve computation."""

    def test_compute_resolution_curve(self):
        """Test resolution curve computation."""
        np.random.seed(42)
        gt = np.random.randint(0, 3, size=(20, 10), dtype=np.int8)
        # Make first 5 SNPs highly informative - perfectly separate varieties
        for i in range(5):
            gt[i] = np.array([0, 2, 0, 2, 0, 2, 0, 2, 0, 2])

        variant_ids = np.array([f'SNP{i}' for i in range(20)])
        chrom = np.array(['chr1'] * 20)
        pos = np.arange(20) * 1000000
        scores = np.concatenate([np.ones(5) * 0.9, np.ones(15) * 0.1])

        panel_sizes = [2, 4, 6, 8, 10]

        curve = compute_resolution_curve(
            gt=gt,
            variant_ids=variant_ids,
            chrom=chrom,
            pos=pos,
            scores=scores,
            panel_sizes=panel_sizes,
            k=1,
            n_random_trials=100,
        )

        assert len(curve) == len(panel_sizes) * 3  # 3 k values
        assert 'designed_fraction' in curve.columns
        assert 'random_mean' in curve.columns
        # Test that function runs and returns correct structure
        assert all(curve['designed_fraction'] >= 0)
        assert all(curve['designed_fraction'] <= 1)
        assert all(curve['random_mean'] >= 0)
        assert all(curve['random_mean'] <= 1)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
