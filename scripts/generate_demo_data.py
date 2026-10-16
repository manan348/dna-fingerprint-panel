"""
Generate demo data for the Streamlit app so users can immediately
see the fingerprint matching workflow in action.

Creates:
- data/demo/reference_panel.parquet   - Small reference panel (10 varieties, 20 SNPs)
- data/demo/query_sample.csv          - Query sample matching variety V01
- data/demo/hybrid_genotype.csv       - F1 hybrid genotype
- data/demo/parent1_genotype.csv       - Parent 1 genotype
- data/demo/parent2_genotype.csv       - Parent 2 genotype
- data/demo/calibrated_thresholds.json - Calibrated decision thresholds
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from fingerprint.io import write_genotype_parquet, GenotypeData


def main():
    rng = np.random.default_rng(42)
    demo_dir = Path("data/demo")
    demo_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------
    # 1. Reference panel: 10 varieties, 20 SNPs
    # ---------------------------------------------------------
    n_varieties = 10
    n_snps = 20
    n_reps = 2

    # Variety IDs
    variety_ids = [f"V{i:02d}" for i in range(1, n_varieties + 1)]

    # SNP IDs
    snp_ids = [f"SNP{j:03d}" for j in range(1, n_snps + 1)]

    # Chromosomes and positions
    chrom = np.array(["chr1"] * 10 + ["chr2"] * 10)
    pos = np.arange(1, n_snps + 1) * 1_000_000
    ref_allele = np.array(["A"] * n_snps)
    alt_allele = np.array(["T"] * n_snps)

    # Generate genotype matrix: (n_snps, n_samples)
    # Each variety has 2 replicates with some noise
    gt = np.zeros((n_snps, n_varieties * n_reps), dtype=np.int8)

    for v_idx, var_id in enumerate(variety_ids):
        # Each variety has a distinct genotype pattern
        base_pattern = rng.integers(0, 3, size=n_snps)
        for r in range(n_reps):
            sample_idx = v_idx * n_reps + r
            # Add some genotyping noise (5% error rate)
            noise = rng.random(n_snps) < 0.05
            gt[:, sample_idx] = base_pattern
            gt[noise, sample_idx] = rng.integers(0, 3, size=noise.sum())

    # Sample IDs
    sample_ids = []
    for var_id in variety_ids:
        for r in range(n_reps):
            sample_ids.append(f"{var_id}_R{r + 1}")
    sample_ids = np.array(sample_ids)

    # Create GenotypeData
    ref_data = GenotypeData(
        gt=gt,
        variant_ids=np.array(snp_ids),
        sample_ids=sample_ids,
        chrom=chrom,
        pos=pos,
        ref=ref_allele,
        alt=alt_allele,
        metadata={"demo": True, "n_varieties": n_varieties, "n_snps": n_snps},
    )

    # Write reference panel
    ref_path = demo_dir / "reference_panel"
    write_genotype_parquet(ref_data, str(ref_path))
    print(f"Reference panel written to {ref_path}.*")

    # ---------------------------------------------------------
    # 2. Query sample: matches V01
    # ---------------------------------------------------------
    query_rows = []
    for j, snp_id in enumerate(snp_ids):
        query_rows.append({"snp_id": snp_id, "genotype": int(gt[j, 0])})

    query_df = pd.DataFrame(query_rows)
    query_path = demo_dir / "query_sample.csv"
    query_df.to_csv(query_path, index=False)
    print(f"Query sample written to {query_path}")

    # ---------------------------------------------------------
    # 3. Hybrid verification demo
    #    Parent 1: all 0s, Parent 2: all 2s, Hybrid: all 1s (true F1)
    # ---------------------------------------------------------
    p1_rows = [{"snp_id": snp_id, "genotype": 0} for snp_id in snp_ids]
    p2_rows = [{"snp_id": snp_id, "genotype": 2} for snp_id in snp_ids]
    h_rows = [{"snp_id": snp_id, "genotype": 1} for snp_id in snp_ids]

    pd.DataFrame(p1_rows).to_csv(demo_dir / "parent1_genotype.csv", index=False)
    pd.DataFrame(p2_rows).to_csv(demo_dir / "parent2_genotype.csv", index=False)
    pd.DataFrame(h_rows).to_csv(demo_dir / "hybrid_genotype.csv", index=False)
    print("Hybrid verification demo files written.")

    # ---------------------------------------------------------
    # 4. Calibrated thresholds
    # ---------------------------------------------------------
    thresholds = {
        "match_threshold": 0.02,
        "related_threshold": 0.10,
        "margin_to_second": 0.05,
        "min_loci_fraction": 0.90,
        "replicate_mean": 0.005,
        "replicate_std": 0.010,
        "nearest_mean": 0.45,
        "nearest_std": 0.080,
        "auc": 0.995,
        "eer": 0.003,
    }
    thresh_path = demo_dir / "calibrated_thresholds.json"
    with open(thresh_path, "w") as f:
        json.dump(thresholds, f, indent=2)
    print(f"Thresholds written to {thresh_path}")

    print("\nDemo data generation complete!")
    print(f"\nDemo files in: {demo_dir.absolute()}")
    print("  - reference_panel.parquet (variants + genotypes + meta)")
    print("  - query_sample.csv")
    print("  - hybrid_genotype.csv")
    print("  - parent1_genotype.csv")
    print("  - parent2_genotype.csv")
    print("  - calibrated_thresholds.json")


if __name__ == "__main__":
    main()