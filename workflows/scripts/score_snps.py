"""
Score SNPs for informativeness and assay-friendliness (workflow script).
"""

import argparse
import numpy as np
import pandas as pd
from fingerprint.io import read_vcf_to_genotype, GenotypeData
from fingerprint.scoring import score_snps, write_scores_table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="Filtered VCF")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--out", required=True, help="Output scores TSV")
    parser.add_argument("--config", required=True, help="Config YAML")
    args = parser.parse_args()

    import yaml
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Read genotypes
    data = read_vcf_to_genotype(args.vcf)

    # Score SNPs
    scores = score_snps(
        gt=data.gt,
        variant_ids=data.variant_ids,
        chrom=data.chrom,
        pos=data.pos,
        ref=data.ref,
        alt=data.alt,
        fasta_path=args.fasta,
        config=config["scoring"],
    )

    # Write output
    write_scores_table(scores, args.out)
    print(f"Scored {len(scores.variant_ids)} SNPs, written to {args.out}")


if __name__ == "__main__":
    main()