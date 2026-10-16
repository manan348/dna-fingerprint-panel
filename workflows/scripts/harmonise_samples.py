"""
Harmonise sample IDs across datasets using alias table.
"""

import sys
import pandas as pd
import numpy as np
from fingerprint.io import read_vcf_to_genotype, write_genotype_parquet, GenotypeData, harmonise_sample_ids, build_alias_table


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--alias-table", required=True)
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--out-metadata", required=True)
    args = parser.parse_args()

    # Read VCF
    data = read_vcf_to_genotype(args.vcf)

    # Read alias table
    alias_df = pd.read_csv(args.alias_table, sep="\t")

    # Harmonise
    harmonised_ids, mapping = harmonise_sample_ids(
        data.sample_ids, alias_df, "alias", "canonical_id"
    )

    # Update metadata
    metadata = pd.DataFrame({
        "sample_id": harmonised_ids,
        "original_id": data.sample_ids,
    })
    metadata.to_csv(args.out_metadata, sep="\t", index=False)

    # Update GenotypeData
    data.sample_ids = harmonised_ids

    # Write back (simplified - just write parquet)
    write_genotype_parquet(data, args.out_vcf.replace(".vcf.gz", ""))

    print(f"Harmonised {len(data.sample_ids)} samples")
    print(f"Metadata written to {args.out_metadata}")


if __name__ == "__main__":
    main()