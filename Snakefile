# DNA Fingerprint Panel - Snakemake Workflow
# Run: snakemake --cores 4 --use-conda
# Or in container: docker run -v $(pwd):/work fingerprint-panel snakemake --cores 4

configfile: "config.yaml"

import yaml
import os

# Load config for use in rules
with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

# ============================================================
# PATHS AND CONSTANTS
# ============================================================
RAW_DIR = cfg["paths"]["raw_data"]
PROC_DIR = cfg["paths"]["processed_data"]
RESULTS_DIR = cfg["paths"]["results_dir"]
PANELS_DIR = cfg["paths"]["panels_dir"]
FIGURES_DIR = cfg["paths"]["figures_dir"]
TABLES_DIR = cfg["paths"]["tables_dir"]

# Ensure output directories exist
for d in [RAW_DIR, PROC_DIR, RESULTS_DIR, PANELS_DIR, FIGURES_DIR, TABLES_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================
# RULES
# ============================================================

rule all:
    input:
        # Final deliverables
        "results/panels/core_panel_snps.txt",
        "results/panels/backup_panel_snps.txt",
        "results/panels/extended_panel_snps.txt",
        "results/tables/candidate_snp_scores.tsv",
        "results/tables/panel_resolution_curve.tsv",
        "results/tables/robustness_results.tsv",
        "results/tables/calibrated_thresholds.json",
        "results/figures/resolution_curve.png",
        "results/figures/panel_map.png",
        "results/figures/distance_distributions.png",
        "results/figures/robustness_heatmap.png",
        "docs/protocol.md",
        "docs/user_guide.md",
        "docs/lab_sop_kasp.md",
        "results/kasp_design_sheet.xlsx",

# ----------------------------------------------------------------------
# STEP 1: Define scope and protocol (manual, but tracked)
# ----------------------------------------------------------------------
rule step1_protocol:
    output:
        "docs/protocol.md"
    shell:
        """
        echo "Step 1: Protocol should be written manually before running pipeline"
        echo "See docs/protocol.md template"
        touch {output}
        """

# ----------------------------------------------------------------------
# STEP 2: Acquire and harmonise data
# ----------------------------------------------------------------------
# Dataset definitions (to be populated in config or here)
DATASETS = {
    "maize_hapmap3": {
        "url": "https://www.maizegenetics.net/hapmap3",
        "format": "hapmap",
        "genome_build": "B73_v3",
    },
    "maize_nam_founders": {
        "url": "https://www.maizegdb.org/data_center/genome",
        "format": "vcf",
        "genome_build": "B73_v4",
    },
    # Add more datasets as needed
}

# Download raw data (manual step - URLs may require authentication)
rule download_data:
    output:
        manifest = "data/MANIFEST.tsv",
    params:
        datasets = DATASETS,
    shell:
        """
        echo "dataset_id\turl\tformat\tgenome_build\tdate_downloaded\tsha256\tlicense" > {output.manifest}
        echo "NOTE: Download datasets manually and record in MANIFEST.tsv"
        echo "See docs/protocol.md for dataset list and access instructions"
        """

# Convert formats to VCF
rule convert_to_vcf:
    input:
        manifest = "data/MANIFEST.tsv",
    output:
        vcf = temp(PROC_DIR + "/{dataset}.vcf.gz"),
    params:
        dataset = "{dataset}",
    shell:
        """
        echo "Converting {wildcards.dataset} to VCF..."
        # Implementation depends on source format
        # HapMap -> TASSEL 5: run_pipeline.pl -fork1 -h {input} -export -exportType VCF -runfork1
        # VCF -> already VCF, just compress and index
        touch {output.vcf}
        """

# Liftover to B73 v5
rule liftover_coordinates:
    input:
        vcf = PROC_DIR + "/{dataset}.vcf.gz",
        chain = "data/reference/liftover_chains/{dataset}.chain",
    output:
        vcf = PROC_DIR + "/{dataset}.b73v5.vcf.gz",
    shell:
        """
        CrossMap.py vcf {input.chain} {input.vcf} {output.vcf} && \
        bcftools index -t {output.vcf}
        """

# Merge datasets
rule merge_datasets:
    input:
        vcfs = expand(PROC_DIR + "/{dataset}.b73v5.vcf.gz", dataset=DATASETS.keys()),
    output:
        vcf = PROC_DIR + "/merged.all.vcf.gz",
        manifest = PROC_DIR + "/merge_manifest.tsv",
    shell:
        """
        bcftools merge {input.vcfs} -Oz -o {output.vcf} && \
        bcftools index -t {output.vcf}
        """

# Harmonise sample IDs
rule harmonise_samples:
    input:
        vcf = PROC_DIR + "/merged.all.vcf.gz",
        alias_table = "data/alias_table.tsv",
    output:
        vcf = PROC_DIR + "/merged.harmonised.vcf.gz",
        metadata = "data/metadata.tsv",
    script:
        "workflows/scripts/harmonise_samples.py"

# ----------------------------------------------------------------------
# STEP 3: Quality control and filtering
# ----------------------------------------------------------------------
rule qc_filter:
    input:
        vcf = PROC_DIR + "/merged.harmonised.vcf.gz",
    output:
        pgen = PROC_DIR + "/filtered.pgen",
        pvar = PROC_DIR + "/filtered.pvar",
        psam = PROC_DIR + "/filtered.psam",
        qc_report = TABLES_DIR + "/qc_report.tsv",
        clone_groups = TABLES_DIR + "/clone_groups.tsv",
    params:
        mind = cfg["qc"]["sample_missing_threshold"],
        geno = cfg["qc"]["snp_missing_threshold"],
        maf = cfg["qc"]["maf_threshold_fingerprint"],
    shell:
        """
        plink2 --vcf {input.vcf} \
               --mind {params.mind} --geno {params.geno} --maf {params.maf} \
               --make-pgen --out {PROC_DIR}/filtered && \
        python workflows/scripts/clone_detection.py \
            --pfile {PROC_DIR}/filtered \
            --out {output.clone_groups} \
            --threshold 0.99
        """

# ----------------------------------------------------------------------
# STEP 4: Score SNP informativeness and assay-friendliness
# ----------------------------------------------------------------------
rule score_snps:
    input:
        pgen = PROC_DIR + "/filtered.pgen",
        pvar = PROC_DIR + "/filtered.pvar",
        psam = PROC_DIR + "/filtered.psam",
        reference_fasta = cfg["paths"]["reference_genome_fasta"],
        repeat_mask = cfg["paths"]["repeat_mask"],
    output:
        scores = TABLES_DIR + "/candidate_snp_scores.tsv",
    script:
        "workflows/scripts/score_snps.py"

# ----------------------------------------------------------------------
# STEP 5: Panel selection (greedy + ILP)
# ----------------------------------------------------------------------
rule select_panel:
    input:
        scores = TABLES_DIR + "/candidate_snp_scores.tsv",
        pgen = PROC_DIR + "/filtered.pgen",
    output:
        core = PANELS_DIR + "/core_panel_snps.txt",
        backup = PANELS_DIR + "/backup_panel_snps.txt",
        extended = PANELS_DIR + "/extended_panel_snps.txt",
        resolution = TABLES_DIR + "/panel_resolution_curve.tsv",
    script:
        "workflows/scripts/select_panel.py"

# ----------------------------------------------------------------------
# STEP 6: Experimental design (documented in protocol)
# ----------------------------------------------------------------------
rule step6_experimental_design:
    input:
        protocol = "docs/protocol.md",
    output:
        TABLES_DIR + "/experimental_design.tsv",
    shell:
        """
        echo "Experimental design documented in protocol.md"
        touch {output}
        """

# ----------------------------------------------------------------------
# STEP 7: Robustness testing
# ----------------------------------------------------------------------
rule robustness_test:
    input:
        core = PANELS_DIR + "/core_panel_snps.txt",
        backup = PANELS_DIR + "/backup_panel_snps.txt",
        pgen = PROC_DIR + "/filtered.pgen",
        pvar = PROC_DIR + "/filtered.pvar",
    output:
        results = TABLES_DIR + "/robustness_results.tsv",
        heatmap = FIGURES_DIR + "/robustness_heatmap.png",
    script:
        "workflows/scripts/robustness_test.py"

# ----------------------------------------------------------------------
# STEP 8: Matching engine and decision rules
# ----------------------------------------------------------------------
rule calibrate_thresholds:
    input:
        core = PANELS_DIR + "/core_panel_snps.txt",
        pgen = PROC_DIR + "/filtered.pgen",
        pvar = PROC_DIR + "/filtered.pvar",
        psam = PROC_DIR + "/filtered.psam",
    output:
        thresholds = TABLES_DIR + "/calibrated_thresholds.json",
        roc_figure = FIGURES_DIR + "/distance_distributions.png",
    script:
        "workflows/scripts/calibrate_thresholds.py"

# ----------------------------------------------------------------------
# STEP 9: Assay design and wet-lab pilot
# ----------------------------------------------------------------------
rule design_kasp_assays:
    input:
        core = PANELS_DIR + "/core_panel_snps.txt",
        backup = PANELS_DIR + "/backup_panel_snps.txt",
        pvar = PROC_DIR + "/filtered.pvar",
        reference_fasta = cfg["paths"]["reference_genome_fasta"],
    output:
        design_sheet = "results/kasp_design_sheet.xlsx",
        pilot_protocol = "docs/lab_sop_kasp.md",
    script:
        "workflows/scripts/design_kasp.py"

# ----------------------------------------------------------------------
# STEP 10: Package and release (documented)
# ----------------------------------------------------------------------
rule step10_package:
    input:
        protocol = "docs/protocol.md",
        user_guide = "docs/user_guide.md",
        lab_sop = "docs/lab_sop_kasp.md",
    output:
        release_notes = "RELEASE.md",
    shell:
        """
        echo "Release preparation documented in RELEASE.md"
        touch {output.release_notes}
        """

# ----------------------------------------------------------------------
# FIGURE GENERATION RULES
# ----------------------------------------------------------------------
rule plot_resolution_curve:
    input:
        TABLES_DIR + "/panel_resolution_curve.tsv",
    output:
        FIGURES_DIR + "/resolution_curve.png",
    script:
        "workflows/scripts/plot_resolution.py"

rule plot_panel_map:
    input:
        PANELS_DIR + "/core_panel_snps.txt",
        PROC_DIR + "/filtered.pvar",
    output:
        FIGURES_DIR + "/panel_map.png",
    script:
        "workflows/scripts/plot_panel_map.py"

# ----------------------------------------------------------------------
# UTILITY RULES
# ----------------------------------------------------------------------
rule clean:
    shell:
        "rm -rf {PROC_DIR}/*.pgen {PROC_DIR}/*.pvar {PROC_DIR}/*.psam {PROC_DIR}/*.vcf.gz"

rule lint:
    shell:
        "ruff check src/ tests/ && black --check src/ tests/"