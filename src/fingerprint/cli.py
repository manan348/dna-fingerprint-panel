"""
Command-line interface for the DNA Fingerprint Panel toolkit.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd
import yaml

from .assay import design_panel_assays, generate_lab_sop, write_kasp_design_sheet
from .io import GenotypeData, load_reference_panel, write_genotype_parquet
from .matcher import (
    Thresholds,
    calibrate_thresholds,
    classify_query,
    seed_lot_purity_test,
    verify_hybrid,
)
from .robustness import recommend_parameters, run_robustness_simulation, summarize_robustness
from .scoring import rank_candidates, read_scores_table, score_snps, write_scores_table
from .selection import compute_resolution_curve, select_panel, write_panel_files


@click.group()
@click.version_option(version="1.0.0")
@click.option("--config", "-c", default="config.yaml", help="Path to config.yaml")
@click.pass_context
def cli(ctx: click.Context, config: str):
    """DNA Fingerprint Panel for Varietal Identification"""
    ctx.ensure_object(dict)
    with open(config) as f:
        ctx.obj["config"] = yaml.safe_load(f)
    ctx.obj["config_path"] = config


@cli.command()
@click.option("--vcf", required=True, help="Input VCF file")
@click.option("--fasta", required=True, help="Reference genome FASTA")
@click.option("--vcf-out", help="Output VCF for candidates")
@click.option("--scores-out", default="results/tables/candidate_snp_scores.tsv", help="Output scores table")
@click.pass_context
def score(ctx: click.Context, vcf: str, fasta: str, vcf_out: str | None, scores_out: str):
    """Score SNPs for informativeness and assay-friendliness (Step 4)."""
    from .io import read_vcf_to_genotype

    config = ctx.obj["config"]
    click.echo(f"Reading VCF: {vcf}")
    data = read_vcf_to_genotype(vcf, max_missing=config["qc"]["snp_missing_threshold"])

    click.echo(f"Scoring {data.n_variants} SNPs across {data.n_samples} samples")
    scores = score_snps(
        gt=data.gt,
        variant_ids=data.variant_ids,
        chrom=data.chrom,
        pos=data.pos,
        ref=data.ref,
        alt=data.alt,
        fasta_path=fasta,
        config=config["scoring"],
    )

    Path(scores_out).parent.mkdir(parents=True, exist_ok=True)
    write_scores_table(scores, scores_out)
    click.echo(f"Scores written to {scores_out}")

    # Rank and show top candidates
    ranked = rank_candidates(scores, top_n=20)
    click.echo("\nTop 20 candidates:")
    click.echo(ranked[["rank", "variant_id", "chrom", "pos", "maf", "pic", "composite_score"]].to_string(index=False))


@cli.command()
@click.option("--scores", required=True, help="Candidate SNP scores table")
@click.option("--vcf", required=True, help="Filtered VCF/genotype file")
@click.option("--out-dir", default="results/panels", help="Output directory for panels")
@click.option("--resolution-out", default="results/tables/panel_resolution_curve.tsv", help="Resolution curve output")
@click.pass_context
def select(ctx: click.Context, scores: str, vcf: str, out_dir: str, resolution_out: str):
    """Select fingerprint panel using greedy + ILP (Step 5)."""
    from .io import read_vcf_to_genotype

    config = ctx.obj["config"]
    sel_config = config["selection"]

    click.echo(f"Reading scores: {scores}")
    score_data = read_scores_table(scores)

    click.echo(f"Reading genotypes: {vcf}")
    data = read_vcf_to_genotype(vcf)

    # Align scores with genotype data
    # This is simplified - in practice need proper alignment
    score_array = score_data.composite_score

    click.echo("Running panel selection...")
    core, backup, extended = select_panel(
        gt=data.gt,
        variant_ids=data.variant_ids,
        chrom=data.chrom,
        pos=data.pos,
        scores=score_array,
        config=sel_config,
    )

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    write_panel_files(core, backup, extended, out_dir)
    click.echo(f"Panels written to {out_dir}")

    # Resolution curve
    panel_sizes = list(range(12, sel_config["core_max_snps"] + 1, 4))
    click.echo("Computing resolution curve...")
    curve = compute_resolution_curve(
        gt=data.gt,
        variant_ids=data.variant_ids,
        chrom=data.chrom,
        pos=data.pos,
        scores=score_array,
        panel_sizes=panel_sizes,
        k=sel_config["k"],
    )
    Path(resolution_out).parent.mkdir(parents=True, exist_ok=True)
    curve.to_csv(resolution_out, sep="\t", index=False)
    click.echo(f"Resolution curve written to {resolution_out}")


@cli.command()
@click.option("--core-panel", required=True, help="Core panel SNP list")
@click.option("--backup-panel", help="Backup panel SNP list")
@click.option("--vcf", required=True, help="Genotype VCF for simulation")
@click.option("--out", default="results/tables/robustness_results.tsv", help="Output table")
@click.option("--n-sims", default=200, help="Simulations per condition")
@click.pass_context
def robustness(ctx: click.Context, core_panel: str, backup_panel: str | None, vcf: str, out: str, n_sims: int):
    """Run robustness simulations (Step 7)."""
    from .io import read_vcf_to_genotype

    config = ctx.obj["config"]
    rob_config = config["robustness"]

    click.echo(f"Reading core panel: {core_panel}")
    core_snps = np.loadtxt(core_panel, dtype=str)

    backup_snps = None
    if backup_panel:
        click.echo(f"Reading backup panel: {backup_panel}")
        backup_snps = np.loadtxt(backup_panel, dtype=str)

    click.echo(f"Reading genotypes: {vcf}")
    data = read_vcf_to_genotype(vcf)

    # Map SNP IDs to indices
    core_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in core_snps if s in data.variant_ids])
    backup_idx = None
    if backup_snps is not None:
        backup_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in backup_snps if s in data.variant_ids])

    click.echo(f"Running {n_sims} simulations per condition...")
    results = run_robustness_simulation(
        gt=data.gt,
        core_panel=core_idx,
        backup_panel=backup_idx,
        config=rob_config,
        n_simulations=n_sims,
        seed=config["random_seed"],
    )

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(out, sep="\t", index=False)
    click.echo(f"Results written to {out}")

    # Summary
    summary = summarize_robustness(results)
    summary_out = out.replace(".tsv", "_summary.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)
    click.echo(f"Summary written to {summary_out}")

    # Recommendations
    rec = recommend_parameters(summary)
    click.echo("\nRecommendations:")
    for k, v in rec.items():
        click.echo(f"  {k}: {v}")


@cli.command()
@click.option("--core-panel", required=True, help="Core panel SNP list")
@click.option("--vcf", required=True, help="Reference genotypes VCF")
@click.option("--replicates", help="TSV mapping variety_id to sample indices for replicates")
@click.option("--out", default="results/tables/calibrated_thresholds.json", help="Output thresholds JSON")
@click.pass_context
def calibrate(ctx: click.Context, core_panel: str, vcf: str, replicates: str | None, out: str):
    """Calibrate matching thresholds from replicate data (Step 8)."""
    from .io import read_vcf_to_genotype

    config = ctx.obj["config"]
    match_config = config["matching"]

    click.echo(f"Reading core panel: {core_panel}")
    core_snps = np.loadtxt(core_panel, dtype=str)

    click.echo(f"Reading genotypes: {vcf}")
    data = read_vcf_to_genotype(vcf)

    core_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in core_snps if s in data.variant_ids])
    panel_gt = data.gt[core_idx]

    # Load replicate groups
    if replicates:
        rep_df = pd.read_csv(replicates, sep="\t")
        replicate_groups = rep_df.groupby("variety_id")["sample_index"].apply(list).to_dict()
    else:
        # Auto-detect from sample IDs (assuming format variety_rep1, variety_rep2)
        replicate_groups = {}
        for i, sid in enumerate(data.sample_ids):
            base = sid.split("_")[0]  # Simple heuristic
            if base not in replicate_groups:
                replicate_groups[base] = []
            replicate_groups[base].append(i)
        # Filter to only those with >1 replicate
        replicate_groups = {k: v for k, v in replicate_groups.items() if len(v) > 1}

    click.echo(f"Calibrating thresholds with {len(replicate_groups)} replicate groups...")
    thresholds = calibrate_thresholds(
        ref_gt=panel_gt,
        ref_ids=data.sample_ids,
        replicate_groups=replicate_groups,
        match_percentile=match_config["match_threshold_percentile"],
        related_multiplier=3.0,
        margin=match_config["margin_to_second_best"],
        min_loci_fraction=match_config["min_loci_called_fraction"],
    )

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
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

    click.echo(f"Thresholds written to {out}")
    click.echo(f"  Match threshold: {thresholds.match_threshold:.4f}")
    click.echo(f"  Related threshold: {thresholds.related_threshold:.4f}")
    click.echo(f"  AUC: {thresholds.auc:.4f}")
    click.echo(f"  EER: {thresholds.eer:.4f}")


@cli.command()
@click.option("--reference", required=True, help="Reference database (Parquet)")
@click.option("--thresholds", required=True, help="Calibrated thresholds JSON")
@click.option("--query", required=True, help="Query genotype CSV (snp_id, genotype)")
@click.option("--output", help="Output JSON result")
@click.pass_context
def match(ctx: click.Context, reference: str, thresholds: str, query: str, output: str | None):
    """Match a query sample against reference panel (Step 8)."""
    with open(thresholds) as f:
        thresh_dict = json.load(f)

    thresholds = Thresholds(**thresh_dict)

    click.echo(f"Loading reference: {reference}")
    ref_gt, ref_var_ids, ref_sample_ids = load_reference_panel(reference)

    click.echo(f"Reading query: {query}")
    qdf = pd.read_csv(query)
    # Build query vector aligned to reference variants
    query_gt = np.full(len(ref_var_ids), -1, dtype=np.int8)
    for _, row in qdf.iterrows():
        idx = np.where(ref_var_ids == row["snp_id"])[0]
        if len(idx) > 0:
            query_gt[idx[0]] = row["genotype"]

    click.echo("Classifying...")
    result = classify_query(query_gt, ref_gt, ref_sample_ids, thresholds)

    output_dict = {
        "decision": result.decision,
        "best_match": result.best_match,
        "best_distance": result.best_distance,
        "n_loci_used": result.n_loci_used,
        "top_hits": result.top_hits,
        "probability_identity": result.probability_identity,
        "warnings": result.warnings,
    }

    if output:
        with open(output, "w") as f:
            json.dump(output_dict, f, indent=2)
        click.echo(f"Result written to {output}")
    else:
        click.echo(json.dumps(output_dict, indent=2))


@cli.command()
@click.option("--panel", required=True, help="Panel SNP list (TSV with variant_id, chrom, pos, ref, alt)")
@click.option("--fasta", required=True, help="Reference genome FASTA")
@click.option("--out-xlsx", default="results/kasp_design_sheet.xlsx", help="Output Excel design sheet")
@click.option("--out-sop", default="docs/lab_sop_kasp.md", help="Output lab SOP")
@click.pass_context
def design_assays(ctx: click.Context, panel: str, fasta: str, out_xlsx: str, out_sop: str):
    """Design KASP assays for a panel (Step 9)."""
    config = ctx.obj["config"]
    assay_config = config["assay"]

    click.echo(f"Reading panel: {panel}")
    panel_df = pd.read_csv(panel, sep="\t")

    click.echo(f"Designing assays for {len(panel_df)} SNPs...")
    assays = design_panel_assays(
        panel_df=panel_df,
        fasta_path=fasta,
        flank_length=assay_config["flank_length"],
        fam_tail=assay_config["universal_tails"]["fam"],
        hex_tail=assay_config["universal_tails"]["hex"],
    )

    Path(out_xlsx).parent.mkdir(parents=True, exist_ok=True)
    write_kasp_design_sheet(assays, out_xlsx)

    Path(out_sop).parent.mkdir(parents=True, exist_ok=True)
    generate_lab_sop(assays, out_sop, panel_name="Core_Panel")

    passed = sum(1 for a in assays if a.pass_all)
    click.echo(f"Designed {len(assays)} assays, {passed} passed all checks")


@cli.command()
@click.option("--core-panel", required=True, help="Core panel SNP list")
@click.option("--vcf", required=True, help="Reference genotypes VCF")
@click.option("--out", default="results/reference_panel.parquet", help="Output reference database")
@click.pass_context
def build_ref(ctx: click.Context, core_panel: str, vcf: str, out: str):
    """Build reference database for matching."""
    from .io import read_vcf_to_genotype

    click.echo(f"Reading core panel: {core_panel}")
    core_snps = np.loadtxt(core_panel, dtype=str)

    click.echo(f"Reading genotypes: {vcf}")
    data = read_vcf_to_genotype(vcf)

    core_idx = np.array([np.where(data.variant_ids == s)[0][0] for s in core_snps if s in data.variant_ids])

    panel_data = GenotypeData(
        gt=data.gt[core_idx],
        variant_ids=data.variant_ids[core_idx],
        sample_ids=data.sample_ids,
        chrom=data.chrom[core_idx],
        pos=data.pos[core_idx],
        ref=data.ref[core_idx],
        alt=data.alt[core_idx],
        metadata={"panel": "core", "n_snps": len(core_idx)},
    )

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    write_genotype_parquet(panel_data, out)
    click.echo(f"Reference database written to {out}")


@cli.command()
@click.option("--hybrid", required=True, help="Hybrid genotype CSV")
@click.option("--parent1", required=True, help="Parent 1 genotype CSV")
@click.option("--parent2", required=True, help="Parent 2 genotype CSV")
@click.pass_context
def verify_hybrid_cmd(ctx: click.Context, hybrid: str, parent1: str, parent2: str):
    """Verify an F1 hybrid against declared parents."""
    h_df = pd.read_csv(hybrid)
    p1_df = pd.read_csv(parent1)
    p2_df = pd.read_csv(parent2)

    # Align by snp_id
    merged = h_df.merge(p1_df, on="snp_id", suffixes=("", "_p1"))
    merged = merged.merge(p2_df, on="snp_id", suffixes=("", "_p2"))

    h_gt = merged["genotype"].values
    p1_gt = merged["genotype_p1"].values
    p2_gt = merged["genotype_p2"].values

    result = verify_hybrid(h_gt, p1_gt, p2_gt)
    click.echo(json.dumps(result, indent=2))


@cli.command()
@click.option("--n-seeds", type=int, required=True, help="Number of seeds tested")
@click.option("--n-offtypes", type=int, required=True, help="Number of off-types detected")
@click.option("--confidence", type=float, default=0.95, help="Confidence level")
@click.pass_context
def seed_purity(ctx: click.Context, n_seeds: int, n_offtypes: int, confidence: float):
    """Calculate seed-lot purity statistics."""
    result = seed_lot_purity_test(n_seeds, n_offtypes, confidence)
    click.echo(json.dumps(result, indent=2))


@cli.command()
@click.pass_context
def run_all(ctx: click.Context):
    """Run the complete pipeline via Snakemake."""
    import subprocess
    result = subprocess.run(["snakemake", "--cores", "4", "--use-conda"], capture_output=True, text=True)
    click.echo(result.stdout)
    if result.stderr:
        click.echo(result.stderr, err=True)
    sys.exit(result.returncode)


def main():
    cli(obj={})


if __name__ == "__main__":
    main()
