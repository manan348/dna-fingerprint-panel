"""
Panel selection using greedy + ILP (workflow script).
"""

import argparse
import numpy as np
import pandas as pd
import yaml
from fingerprint.io import read_vcf_to_genotype
from fingerprint.selection import select_panel, write_panel_files
from fingerprint.scoring import read_scores_table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, help="Candidate SNP scores TSV")
    parser.add_argument("--vcf", required=True, help="Filtered VCF")
    parser.add_argument("--out-dir", required=True, help="Output directory for panels")
    parser.add_argument("--resolution-out", required=True, help="Resolution curve output")
    parser.add_argument("--config", required=True, help="Config YAML")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Read scores
    score_data = read_scores_table(args.scores)
    score_array = score_data.composite_score

    # Read genotypes
    data = read_vcf_to_genotype(args.vcf)

    # Align scores with genotype variants
    # Create mapping from variant_id to index in genotype data
    var_to_idx = {v: i for i, v in enumerate(data.variant_ids)}
    aligned_scores = np.zeros(len(data.variant_ids))
    for i, v in enumerate(score_data.variant_ids):
        if v in var_to_idx:
            aligned_scores[var_to_idx[v]] = score_data.composite_score[i]

    # Select panel
    core, backup, extended = select_panel(
        gt=data.gt,
        variant_ids=data.variant_ids,
        chrom=data.chrom,
        pos=data.pos,
        scores=aligned_scores,
        config=config["selection"],
    )

    # Write panels
    write_panel_files(core, backup, extended, args.out_dir)

    # Resolution curve
    from fingerprint.selection import compute_resolution_curve
    panel_sizes = list(range(12, config["selection"]["core_max_snps"] + 1, 4))
    curve = compute_resolution_curve(
        gt=data.gt,
        variant_ids=data.variant_ids,
        chrom=data.chrom,
        pos=data.pos,
        scores=aligned_scores,
        panel_sizes=panel_sizes,
        k=config["selection"]["k"],
    )
    curve.to_csv(args.resolution_out, sep="\t", index=False)
    print(f"Resolution curve written to {args.resolution_out}")


if __name__ == "__main__":
    main()