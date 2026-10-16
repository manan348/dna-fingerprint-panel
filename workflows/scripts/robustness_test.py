"""
Robustness testing (workflow script).
"""

import argparse
import numpy as np
import pandas as pd
import yaml
from fingerprint.io import read_vcf_to_genotype
from fingerprint.robustness import run_robustness_simulation, summarize_robustness, recommend_parameters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", required=True, help="Core panel SNP list")
    parser.add_argument("--backup", help="Backup panel SNP list")
    parser.add_argument("--vcf", required=True, help="Genotype VCF")
    parser.add_argument("--out", required=True, help="Output results TSV")
    parser.add_argument("--heatmap", required=True, help="Output heatmap PNG")
    parser.add_argument("--config", required=True, help="Config YAML")
    parser.add_argument("--n-sims", type=int, default=200, help="Simulations per condition")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Read core panel
    core_snps = np.loadtxt(args.core, dtype=str)

    backup_snps = None
    if args.backup:
        backup_snps = np.loadtxt(args.backup, dtype=str)

    # Read genotypes
    data = read_vcf_to_genotype(args.vcf)

    # Map SNP IDs to indices
    core_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in core_snps if s in data.variant_ids])
    backup_idx = None
    if backup_snps is not None:
        backup_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in backup_snps if s in data.variant_ids])

    print(f"Running robustness simulation with {len(core_idx)} core SNPs...")
    results = run_robustness_simulation(
        gt=data.gt,
        core_panel=core_idx,
        backup_panel=backup_idx,
        config=config["robustness"],
        n_simulations=args.n_sims,
        seed=config["random_seed"],
    )

    results.to_csv(args.out, sep="\t", index=False)
    print(f"Results written to {args.out}")

    # Summary
    summary = summarize_robustness(results)
    summary_out = args.out.replace(".tsv", "_summary.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"Summary written to {summary_out}")

    # Recommendations
    rec = recommend_parameters(summary)
    print("\nRecommendations:")
    for k, v in rec.items():
        print(f"  {k}: {v}")

    # Heatmap
    try:
        from fingerprint.robustness import plot_robustness_heatmap
        plot_robustness_heatmap(results, args.heatmap)
        print(f"Heatmap written to {args.heatmap}")
    except ImportError:
        print("Matplotlib not available, skipping heatmap")


if __name__ == "__main__":
    main()