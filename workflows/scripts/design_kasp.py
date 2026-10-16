"""
Design KASP assays for panel SNPs (workflow script).
"""

import argparse
import pandas as pd
import yaml
from fingerprint.assay import design_panel_assays, write_kasp_design_sheet, generate_lab_sop


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", required=True, help="Core panel SNP list")
    parser.add_argument("--backup", help="Backup panel SNP list")
    parser.add_argument("--pvar", required=True, help="PLINK2 .pvar file with variant info")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--out-xlsx", required=True, help="Output Excel design sheet")
    parser.add_argument("--out-sop", required=True, help="Output lab SOP markdown")
    parser.add_argument("--config", required=True, help="Config YAML")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Read panel SNPs
    core_snps = pd.read_csv(args.core, header=None, names=["variant_id"])

    # Read variant info from .pvar
    pvar = pd.read_csv(args.pvar, sep="\t", comment="#")
    # Columns: CHROM, POS, ID, REF, ALT
    pvar = pvar.rename(columns={"#CHROM": "chrom", "POS": "pos", "ID": "variant_id", "REF": "ref", "ALT": "alt"})

    # Merge to get coordinates
    core_df = core_snps.merge(pvar[["variant_id", "chrom", "pos", "ref", "alt"]], on="variant_id", how="left")

    if core_df["chrom"].isna().any():
        missing = core_df[core_df["chrom"].isna()]["variant_id"].tolist()
        print(f"Warning: {len(missing)} SNPs not found in .pvar: {missing[:10]}")

    core_df = core_df.dropna(subset=["chrom"]).reset_index(drop=True)

    # Add backup SNPs if provided
    if args.backup:
        backup_snps = pd.read_csv(args.backup, header=None, names=["variant_id"])
        backup_df = backup_snps.merge(pvar[["variant_id", "chrom", "pos", "ref", "alt"]], on="variant_id", how="left")
        backup_df = backup_df.dropna(subset=["chrom"])
        panel_df = pd.concat([core_df, backup_df], ignore_index=True)
    else:
        panel_df = core_df

    print(f"Designing assays for {len(panel_df)} SNPs...")
    assays = design_panel_assays(
        panel_df=panel_df,
        fasta_path=args.fasta,
        flank_length=config["assay"]["flank_length"],
        fam_tail=config["assay"]["universal_tails"]["fam"],
        hex_tail=config["assay"]["universal_tails"]["hex"],
    )

    # Write design sheet
    write_kasp_design_sheet(assays, args.out_xlsx)

    # Write lab SOP
    generate_lab_sop(assays, args.out_sop, panel_name="Core_Panel")

    passed = sum(1 for a in assays if a.pass_all)
    print(f"Designed {len(assays)} assays, {passed} passed all checks")


if __name__ == "__main__":
    main()