"""
Plot resolution curve (workflow script).
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Resolution curve TSV")
    parser.add_argument("--output", required=True, help="Output PNG")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, k in enumerate([1, 2, 3]):
        ax = axes[idx]
        subset = df[df["k"] == k]

        ax.plot(subset["panel_size"], subset["designed_fraction"], 'b-o', label='Designed', linewidth=2)
        ax.plot(subset["panel_size"], subset["random_mean"], 'r--', label='Random (mean)', linewidth=1)
        ax.fill_between(
            subset["panel_size"],
            subset["random_min"],
            subset["random_max"],
            alpha=0.2, color='red', label='Random (range)'
        )

        ax.set_xlabel('Panel Size')
        ax.set_ylabel(f'Fraction of Pairs Separated (k≥{k})')
        ax.set_title(f'Separation at k={k}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1.05)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Resolution curve plot saved to {args.output}")


if __name__ == "__main__":
    main()