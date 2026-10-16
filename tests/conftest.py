"""
Pytest configuration and shared fixtures.
"""

import numpy as np
import pytest


@pytest.fixture
def sample_genotype_data():
    """Create sample genotype data for testing."""
    np.random.seed(42)
    # 20 variants, 10 samples
    gt = np.random.randint(0, 3, size=(20, 10), dtype=np.int8)
    # Add some missing
    gt[gt == 3] = -1  # Not possible but ensures some -1

    variant_ids = np.array([f'SNP{i:03d}' for i in range(20)])
    sample_ids = np.array([f'VAR{i//2}_REP{i%2}' for i in range(10)])
    chrom = np.array(['chr1'] * 10 + ['chr2'] * 10)
    pos = np.arange(20) * 1000000
    ref = np.random.choice(['A', 'T', 'C', 'G'], size=20)
    alt = np.random.choice(['A', 'T', 'C', 'G'], size=20)

    return {
        'gt': gt,
        'variant_ids': variant_ids,
        'sample_ids': sample_ids,
        'chrom': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
    }


@pytest.fixture
def perfect_genotype_data():
    """Create perfectly separated genotype data."""
    # 3 varieties, 2 reps each, 10 SNPs
    gt = np.zeros((10, 6), dtype=np.int8)
    gt[:, 0:2] = 0   # Var 0
    gt[:, 2:4] = 2   # Var 1
    gt[:, 4:6] = 1   # Var 2 (all het)

    variant_ids = np.array([f'SNP{i:03d}' for i in range(10)])
    sample_ids = np.array(['V0_R1', 'V0_R2', 'V1_R1', 'V1_R2', 'V2_R1', 'V2_R2'])
    chrom = np.array(['chr1'] * 10)
    pos = np.arange(10) * 1000000
    ref = np.array(['A'] * 10)
    alt = np.array(['T'] * 10)

    return {
        'gt': gt,
        'variant_ids': variant_ids,
        'sample_ids': sample_ids,
        'chrom': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
    }


@pytest.fixture
def replicate_groups():
    """Standard replicate groups for testing."""
    return {
        'V0': [0, 1],
        'V1': [2, 3],
        'V2': [4, 5],
    }


@pytest.fixture(scope="session")
def test_config():
    """Test configuration."""
    return {
        "qc": {
            "sample_missing_threshold": 0.10,
            "snp_missing_threshold": 0.05,
            "maf_threshold_fingerprint": 0.20,
            "maf_threshold_close_relatives": 0.05,
            "het_threshold_inbreds": 0.05,
            "call_rate_threshold": 0.95,
        },
        "scoring": {
            "pic_weight": 0.4,
            "he_weight": 0.3,
            "flank_uniqueness_weight": 0.2,
            "spacing_weight": 0.1,
        },
        "selection": {
            "k": 3,
            "core_max_snps": 96,
            "backup_fraction": 0.20,
            "extended_max_snps": 384,
            "min_gap_bp": 2_000_000,
            "max_chrom_fraction": 0.15,
            "ilp_time_limit": 300,
            "min_maf": 0.20,
        },
        "robustness": {
            "error_rates": [0.0, 0.01, 0.02, 0.05],
            "missing_rates": [0.0, 0.05, 0.10, 0.15],
            "mixture_props": [0.0, 0.05, 0.10, 0.20],
            "dropout_rates": [0.0, 0.05, 0.10, 0.15],
            "n_simulations": 10,  # Small for testing
        },
        "matching": {
            "min_loci_called_fraction": 0.90,
            "match_threshold_percentile": 99.9,
            "related_threshold": 0.15,
            "margin_to_second_best": 0.05,
        },
        "assay": {
            "amplicon_size_range": [50, 100],
            "primer_tm_range": [57, 63],
            "primer_gc_range": [40, 60],
            "flank_length": 60,
            "universal_tails": {
                "fam": "GAAGGTGACCAAGTTCATGCT",
                "hex": "GAAGGTCGGAGTCAACGGATT",
            },
        },
        "random_seed": 42,
    }


def pytest_configure(config):
    """Configure pytest."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
