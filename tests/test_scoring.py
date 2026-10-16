"""
Unit tests for SNP scoring module.
"""

import numpy as np
import pytest

from fingerprint.scoring import calculate_informativeness, score_snps


class TestCalculateInformativeness:
    """Tests for calculate_informativeness function."""

    def test_basic_calculation(self):
        """Test basic informativeness calculations."""
        # Create simple genotype matrix
        # 3 variants, 4 samples
        gt = np.array([
            [0, 0, 2, 2],  # MAF = 0.5
            [0, 0, 0, 2],  # MAF = 0.125
            [1, 1, 1, 1],  # All het, MAF = 0.5 but high het
        ], dtype=np.int8)

        variant_ids = np.array(['SNP1', 'SNP2', 'SNP3'])

        result = calculate_informativeness(gt, variant_ids, min_called=1)

        # Check shapes
        assert result['maf'].shape == (3,)
        assert result['pic'].shape == (3,)
        assert result['he'].shape == (3,)
        assert result['pi_locus'].shape == (3,)

        # SNP1: MAF=0.5, PIC=0.375, HE=0.5, PI=0.5
        assert np.isclose(result['maf'][0], 0.5)
        assert np.isclose(result['pic'][0], 0.375)
        assert np.isclose(result['he'][0], 0.5)
        assert np.isclose(result['pi_locus'][0], 0.5)

        # SNP2: [0, 0, 0, 2] -> alt_count=2, alt_freq=2/8=0.25, maf=0.25
        assert np.isclose(result['maf'][1], 0.25)

        # SNP3: All het -> alt_freq=0.5, but het_rate=1.0
        assert np.isclose(result['maf'][2], 0.5)
        assert np.isclose(result['het_rate'][2], 1.0)

    def test_missing_data_handling(self):
        """Test handling of missing genotypes."""
        gt = np.array([
            [0, 0, 2, -1],  # One missing
            [-1, -1, -1, -1],  # All missing
        ], dtype=np.int8)

        variant_ids = np.array(['SNP1', 'SNP2'])
        result = calculate_informativeness(gt, variant_ids, min_called=1)

        # SNP1: 3 called, alt_freq = 2/(2*3) = 0.333
        assert np.isclose(result['maf'][0], 1/3)
        assert result['call_rate'][0] == 0.75

        # SNP2: 0 called -> should be masked
        assert result['call_rate'][1] == 0.0

    def test_monomorphic_snp(self):
        """Test monomorphic SNP (MAF=0)."""
        gt = np.array([
            [0, 0, 0, 0],  # All same
        ], dtype=np.int8)

        variant_ids = np.array(['SNP1'])
        result = calculate_informativeness(gt, variant_ids)

        assert result['maf'][0] == 0.0
        assert result['pic'][0] == 0.0
        assert result['he'][0] == 0.0
        assert result['pi_locus'][0] == 1.0


class TestScoreSNPs:
    """Tests for score_snps function (requires FASTA)."""

    @pytest.fixture
    def mock_fasta(self, tmp_path):
        """Create a mock FASTA file for testing."""
        fasta_path = tmp_path / "test.fa"
        fasta_path.write_text(">chr1\n" + "A" * 1000 + "\n")
        # Create index
        import subprocess
        subprocess.run(["samtools", "faidx", str(fasta_path)], check=True, capture_output=True)
        return str(fasta_path)

    def test_score_snps_basic(self, mock_fasta):
        """Test basic scoring with mock FASTA."""
        gt = np.array([
            [0, 2, 0, 2],  # Good SNP
            [0, 0, 0, 2],  # Low MAF
        ], dtype=np.int8)

        variant_ids = np.array(['SNP1', 'SNP2'])
        chrom = np.array(['chr1', 'chr1'])
        pos = np.array([100, 200])
        ref = np.array(['A', 'T'])
        alt = np.array(['G', 'C'])

        scores = score_snps(
            gt=gt,
            variant_ids=variant_ids,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            fasta_path=mock_fasta,
            config={
                "min_maf": 0.20,
                "max_het": 0.05,
                "min_call_rate": 0.95,
                "flank_length": 50,
            }
        )

        assert len(scores.variant_ids) == 2
        assert scores.composite_score[0] > 0  # Good SNP should have positive score
        assert scores.composite_score[1] == 0  # Low MAF should be zeroed


class TestEdgeCases:
    """Edge case tests."""

    def test_single_sample(self):
        """Test with single sample."""
        gt = np.array([[0, 2, 1]], dtype=np.int8).T  # (3 variants, 1 sample)
        variant_ids = np.array(['SNP1', 'SNP2', 'SNP3'])

        result = calculate_informativeness(gt, variant_ids)
        assert result['maf'].shape == (3,)

    def test_all_missing(self):
        """Test all missing genotypes."""
        gt = np.full((2, 4), -1, dtype=np.int8)
        variant_ids = np.array(['SNP1', 'SNP2'])

        result = calculate_informativeness(gt, variant_ids, min_called=1)
        assert np.all(result['maf'] == 0)
        assert np.all(result['call_rate'] == 0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
