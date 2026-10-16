"""
Plot panel map showing SNP positions on chromosomes (workflow script).
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", required=True, help="Core panel SNP list")
    parser.add_argument("--pvar", required=True, help="PLINK2 .pvar file")
    parser.add_argument("--output", required=True, help="Output PNG")
    args = parser.parse_args()

    # Read panel SNPs
    panel_snps = np.loadtxt(args.panel, dtype=str)

    # Read variant info
    pvar = pd.read_csv(args.pvar, sep="\t", comment="#")
    pvar = pvar.rename(columns={"#CHROM": "chrom", "POS": "pos", "ID": "variant_id"})

    # Merge
    panel_df = pd.DataFrame({"variant_id": panel_snps})
    panel_df = panel_df.merge(pvar[["variant_id", "chrom", "pos"]], on="variant_id", how="left")

    # Maize chromosome order
    chrom_order = [f"chr{i}" for i in range(1, 11)] + ["chrUn"]
    chrom_lengths = {
        "chr1": 307, "chr2": 244, "chr3": 235, "chr4": 246,
        "chr5": 223, "chr6": 174, "chr7": 182, "chr8": 181,
        "chr9": 159, "chr10": 150,
    }

    fig, ax = plt.subplots(figsize=(14, 6))

    y_pos = 0
    for chrom in chrom_order:
        if chrom not in panel_df["chrom"].values:
            continue

        chrom_snps = panel_df[panel_df["chrom"] == chrom].sort_values("pos")
        length = chrom_lengths.get(chrom, 300)  # Mb

        # Draw chromosome bar
        ax.barh(y_pos, length, height=0.6, left=0, color='lightgray', edgecolor='gray', alpha=0.5)

        # Draw SNPs as ticks
        positions = chrom_snps["pos"].values / 1e6  # Convert to Mb
        ax.scatter(positions, [y_pos] * len(positions), s=30, c='red', marker='|', linewidth=1.5, zorder=5)

        # Label
        ax.text(-10, y_pos, chrom, ha='right', va='center', fontsize=10, fontweight='bold')
        y_pos += 1

    ax.set_xlabel('Physical Position (Mb)')
    ax.set_yticks([])
    ax.set_title(f'Fingerprint Panel Map ({len(panel_df)} SNPs)')
    ax.set_xlim(-20, max(chrom_lengths.values()) + 10)
    ax.grid(True, axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Panel map saved to {args.output}")


if __name__ == "__main__":
    main()