"""
Clone detection from PLINK2 files.
"""

import argparse
import numpy as np
import pandas as pd
import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pfile", required=True, help="PLINK2 file prefix")
    parser.add_argument("--out", required=True, help="Output clone groups TSV")
    parser.add_argument("--threshold", type=float, default=0.99, help="IBS threshold")
    parser.add_argument("--min-shared", type=int, default=10000, help="Minimum shared SNPs")
    args = parser.parse_args()

    # Read PLINK2 files
    pvar = pd.read_csv(f"{args.pfile}.pvar", sep="\t", comment="#")
    psam = pd.read_csv(f"{args.pfile}.psam", sep="\t", comment="#")

    # Read genotype data using pysam (pgen)
    # For efficiency, we'll sample SNPs
    n_variants = len(pvar)
    n_samples = len(psam)

    # Sample SNPs for IBS calculation
    if n_variants > args.min_shared:
        sample_idx = np.random.choice(n_variants, args.min_shared, replace=False)
        sample_idx = np.sort(sample_idx)
    else:
        sample_idx = np.arange(n_variants)

    # Read sampled genotypes
    # Note: This is simplified; in practice use plink2 --export or a proper pgen reader
    gt = np.zeros((len(sample_idx), n_samples), dtype=np.int8)

    # Use plink2 to export to VCF for the sampled SNPs
    # This is a placeholder - real implementation would read pgen directly
    import subprocess
    snp_ids = pvar.iloc[sample_idx]["ID"].tolist()
    with open("temp_snps.txt", "w") as f:
        f.write("\n".join(snp_ids))

    subprocess.run([
        "plink2", "--pfile", args.pfile,
        "--extract", "temp_snps.txt",
        "--export", "vcf",
        "--out", "temp_export"
    ], check=True, capture_output=True)

    # Read VCF
    import allel
    callset = allel.read_vcf("temp_export.vcf.gz", fields=["calldata/GT"])
    gt_array = allel.GenotypeArray(callset["calldata/GT"])
    gt = gt_array.to_n_alt(fill=-1).astype(np.int8).T  # (n_variants, n_samples)

    # Cleanup
    import os
    for f in ["temp_snps.txt", "temp_export.vcf.gz", "temp_export.vcf.gz.tbi"]:
        if os.path.exists(f):
            os.remove(f)

    # Compute pairwise IBS
    sample_ids = psam["#IID"].values if "#IID" in psam.columns else psam.iloc[:, 0].values
    clone_pairs = []

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            both_called = (gt[:, i] >= 0) & (gt[:, j] >= 0)
            n_shared = both_called.sum()
            if n_shared < args.min_shared:
                continue
            matches = (gt[both_called, i] == gt[both_called, j]).sum()
            ibs = matches / n_shared
            if ibs >= args.threshold:
                clone_pairs.append({
                    "sample1": sample_ids[i],
                    "sample2": sample_ids[j],
                    "ibs": ibs,
                    "n_shared": n_shared,
                })

    if clone_pairs:
        pairs_df = pd.DataFrame(clone_pairs)
        groups = assign_clone_groups(pairs_df)
        groups.to_csv(args.out, sep="\t", index=False)
        print(f"Found {len(groups)} clone group assignments")
    else:
        pd.DataFrame(columns=["group_id", "sample", "representative"]).to_csv(args.out, sep="\t", index=False)
        print("No clone pairs found")


def assign_clone_groups(pairs_df):
    from collections import defaultdict
    adj = defaultdict(set)
    for _, row in pairs_df.iterrows():
        adj[row["sample1"]].add(row["sample2"])
        adj[row["sample2"]].add(row["sample1"])

    visited = set()
    groups = []
    group_id = 0

    for sample in adj:
        if sample in visited:
            continue
        component = set()
        stack = [sample]
        while stack:
            s = stack.pop()
            if s in visited:
                continue
            visited.add(s)
            component.add(s)
            stack.extend(adj[s] - visited)

        rep = sorted(component)[0]
        for s in component:
            groups.append({
                "group_id": group_id,
                "sample": s,
                "representative": rep,
            })
        group_id += 1

    return pd.DataFrame(groups)


if __name__ == "__main__":
    main()