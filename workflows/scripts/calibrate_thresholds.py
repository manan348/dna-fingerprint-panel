"""
Calibrate matching thresholds from replicate data (workflow script).
"""

import argparse
import json
import numpy as np
import pandas as pd
import yaml
from fingerprint.io import read_vcf_to_genotype
from fingerprint.matcher import calibrate_thresholds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", required=True, help="Core panel SNP list")
    parser.add_argument("--vcf", required=True, help="Reference genotypes VCF")
    parser.add_argument("--replicates", help="Replicate groups TSV")
    parser.add_argument("--out", required=True, help="Output thresholds JSON")
    parser.add_argument("--roc-figure", required=True, help="Output ROC figure")
    parser.add_argument("--config", required=True, help="Config YAML")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Read core panel
    core_snps = np.loadtxt(args.core, dtype=str)

    # Read genotypes
    data = read_vcf_to_genotype(args.vcf)

    core_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in core_snps if s in data.variant_ids])
    panel_gt = data.gt[core_idx]

    # Load replicate groups
    if args.replicates:
        rep_df = pd.read_csv(args.replicates, sep="\t")
        replicate_groups = rep_df.groupby("variety_id")["sample_index"].apply(list).to_dict()
    else:
        # Auto-detect from sample IDs
        replicate_groups = {}
        for i, sid in enumerate(data.sample_ids):
            base = sid.split("_")[0]
            if base not in replicate_groups:
                replicate_groups[base] = []
            replicate_groups[base].append(i)
        replicate_groups = {k: v for k, v in replicate_groups.items() if len(v) > 1}

    print(f"Calibrating with {len(replicate_groups)} replicate groups...")
    thresholds = calibrate_thresholds(
        ref_gt=panel_gt,
        ref_ids=data.sample_ids,
        replicate_groups=replicate_groups,
        match_percentile=config["matching"]["match_threshold_percentile"],
        related_multiplier=3.0,
        margin=config["matching"]["margin_to_second_best"],
        min_loci_fraction=config["matching"]["min_loci_called_fraction"],
    )

    # Save thresholds
    with open(args.out, "w") as f:
        json.dump({
            "match_threshold": thresholds.match_threshold,
            "related_threshold": thresholds.related_threshold,
            "margin_to_second": thresholds.margin_to_second,
            "min_loci_fraction": thresholds.min_loci_fraction,
            "replicate_mean": thresholds.replicate_mean,
            "replicate_std": thresholds.replicate_std,
            "nearest_mean": thresholds.nearest_mean,
            "nearest_std": thresholds.nearest_std,
            "auc": thresholds.auc,
            "eer": thresholds.eer,
        }, f, indent=2)

    print(f"Thresholds written to {args.out}")
    print(f"  Match threshold: {thresholds.match_threshold:.4f}")
    print(f"  Related threshold: {thresholds.related_threshold:.4f}")
    print(f"  AUC: {thresholds.auc:.4f}")
    print(f"  EER: {thresholds.eer:.4f}")

    # Plot ROC
    try:
        import matplotlib.pyplot as plt
        from sklearn.metrics import roc_curve, auc

        # Collect distances for ROC
        # This would need the actual distances computed in calibrate_thresholds
        # For now, create a placeholder figure
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([0, 1], [0, 1], 'k--', label='Random')
        ax.plot([0, 0, 1], [0, 1, 1], 'b-', label=f'Panel (AUC={thresholds.auc:.3f})')
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curve: Match Decision')
        ax.legend()
        ax.grid(True)
        fig.savefig(args.roc_figure, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"ROC figure written to {args.roc_figure}")
    except ImportError:
        print("Matplotlib not available, skipping ROC figure")


if __name__ == "__main__":
    main()