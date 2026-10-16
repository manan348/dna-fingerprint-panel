# User Guide: DNA Fingerprint Panel Toolkit

## Table of Contents
1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Data Preparation](#data-preparation)
4. [Running the Pipeline](#running-the-pipeline)
5. [Using the Matcher](#using-the-matcher)
6. [Web Interface](#web-interface)
7. [Interpreting Results](#interpreting-results)
8. [Troubleshooting](#troubleshooting)

---

## Installation

### Option 1: Conda (Recommended)

```bash
# Clone repository
git clone https://github.com/neurocrop/dna-fingerprint-panel.git
cd dna-fingerprint-panel

# Create environment
conda env create -f environment.yml
conda activate fingerprint-panel

# Install package
pip install -e .
```

### Option 2: Docker

```bash
# Build image
docker build -t fingerprint-panel:1.0.0 .

# Run container
docker run -it --rm -v $(pwd):/work fingerprint-panel:1.0.0
# Inside container:
fingerprint --help
```

### Option 3: Minimal pip

```bash
pip install git+https://github.com/neurocrop/dna-fingerprint-panel.git
```

---

## Quick Start

### 1. Prepare Reference Data

```bash
# Download public maize datasets to data/raw/
# See docs/protocol.md for dataset list and URLs

# Record in manifest
# data/MANIFEST.tsv should contain:
# dataset_id	url	format	genome_build	date_downloaded	sha256	license
```

### 2. Run Full Pipeline

```bash
# Via Snakemake (recommended)
snakemake --cores 4 --use-conda

# Or step by step:
fingerprint score --vcf data/processed/merged.vcf.gz --fasta data/reference/B73_v5.fa
fingerprint select --scores results/tables/candidate_snp_scores.tsv --vcf data/processed/filtered.vcf.gz
fingerprint robustness --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz
fingerprint calibrate --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz
fingerprint build-ref --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz
fingerprint design-assays --panel results/panels/core_panel_details.tsv --fasta data/reference/B73_v5.fa
```

### 3. Match a Sample

```bash
# Prepare query CSV with columns: snp_id, genotype
# Genotype coding: 0=hom_ref, 1=het, 2=hom_alt, -1=missing

fingerprint match \
    --reference results/reference_panel.parquet \
    --thresholds results/tables/calibrated_thresholds.json \
    --query my_sample.csv \
    --output result.json
```

---

## Data Preparation

### Required Files

| File | Format | Description |
|------|--------|-------------|
| Reference genome | FASTA (.fa, .fai) | B73 RefGen_v5 for maize |
| Repeat mask | BED | Low-mappability regions |
| Genotype data | VCF/BCF | Public datasets |
| Alias table | TSV | Sample ID harmonisation |

### Manifest File

`data/MANIFEST.tsv` tracks all data sources:

```tsv
dataset_id	url	format	genome_build	date_downloaded	sha256	license
maize_hapmap3	https://...	hapmap	B73_v3	2026-10-01	abc123	CC-BY-4.0
maize_nam	https://...	vcf	B73_v4	2026-10-01	def456	CC-BY-4.0
```

### Metadata File

`data/metadata.tsv` (auto-generated, then curated):

```tsv
variety_id	accepted_name	aliases	type	pedigree	heterotic_group	release_year	country	source_dataset	platform	genome_build	n_replicates	licence
B73	B73	B73,Mo17	inbred	-	SS	1972	USA	maize_hapmap3	GBS	B73_v5	3	CC-BY-4.0
```

### Alias Table

`data/alias_table.tsv` for ID harmonisation:

```tsv
alias	canonical_id
B73_v3	B73
B73_REF	B73
Zm00001d000001	B73
```

---

## Running the Pipeline

### Snakemake Workflow

The `Snakefile` defines all steps with dependencies:

```bash
# Dry run
snakemake -n --cores 4

# Run specific target
snakemake results/panels/core_panel_snps.txt --cores 4

# Run with conda environments
snakemake --cores 4 --use-conda

# Force re-run
snakemake -F --cores 4
```

### Key Targets

| Target | Description |
|--------|-------------|
| `results/panels/core_panel_snps.txt` | Core panel SNP list |
| `results/tables/candidate_snp_scores.tsv` | All SNP scores |
| `results/tables/panel_resolution_curve.tsv` | Resolution curve |
| `results/tables/robustness_results.tsv` | Simulation results |
| `results/tables/calibrated_thresholds.json` | Decision thresholds |
| `results/reference_panel.parquet` | Reference database |
| `results/kasp_design_sheet.xlsx` | KASP primer designs |
| `docs/lab_sop_kasp.md` | Lab protocol |

### Configuration

All parameters in `config.yaml`:

```yaml
panel:
  core_size_range: [48, 96]
  backup_fraction: 0.20

separation:
  k: 3
  min_snp_spacing_bp: 2000000

qc:
  maf_threshold_fingerprint: 0.20
  het_threshold_inbreds: 0.05

matching:
  min_loci_called_fraction: 0.90
  match_threshold_percentile: 99.9
```

---

## Using the Matcher

### Command Line

```bash
# Match single sample
fingerprint match \
    --reference results/reference_panel.parquet \
    --thresholds results/tables/calibrated_thresholds.json \
    --query query.csv

# Verify F1 hybrid
fingerprint verify-hybrid \
    --hybrid hybrid.csv \
    --parent1 parent1.csv \
    --parent2 parent2.csv

# Seed lot purity
fingerprint seed-purity --n-seeds 100 --n-offtypes 2
```

### Query Format

CSV with columns:
```csv
snp_id,genotype
SNP_001,0
SNP_002,2
SNP_003,1
SNP_004,-1
```

Genotype coding:
- `0` = homozygous reference
- `1` = heterozygous
- `2` = homozygous alternate
- `-1` = missing

### Python API

```python
from fingerprint.matcher import load_reference, classify_query, Thresholds
from fingerprint.io import load_reference_panel

# Load reference
ref_gt, ref_var_ids, ref_sample_ids = load_reference_panel("results/reference_panel.parquet")

# Load thresholds
with open("results/tables/calibrated_thresholds.json") as f:
    thresh = Thresholds(**json.load(f))

# Prepare query
query_gt = np.full(len(ref_var_ids), -1, dtype=np.int8)
for snp_id, gt in query_dict.items():
    idx = np.where(ref_var_ids == snp_id)[0]
    if len(idx) > 0:
        query_gt[idx[0]] = gt

# Classify
result = classify_query(query_gt, ref_gt, ref_sample_ids, thresh)
print(result.decision, result.best_match, result.best_distance)
```

---

## Web Interface

Launch the Streamlit app:

```bash
streamlit run app/streamlit_app.py
```

### Features

| Tab | Function |
|-----|----------|
| **Match Sample** | Upload query CSV, get variety match |
| **Verify Hybrid** | Check F1 against declared parents |
| **Seed Purity** | Calculate off-type rate & CI |
| **Design Assays** | Design KASP primers for panel |
| **Reference Panel** | Browse varieties and SNPs |

### Screenshots

*(Add screenshots here after running)*

---

## Interpreting Results

### Match Decision Classes

| Decision | Meaning | Action |
|----------|---------|--------|
| **MATCH** | Sample matches reference variety | Confirm identity |
| **AMBIGUOUS/RELATED** | Close to reference but not exact | Check for relatives, mixtures |
| **NO MATCH** | No close match in database | New variety or foreign seed |
| **INSUFFICIENT DATA** | <90% loci called | Re-genotype sample |

### Key Metrics

- **Distance**: Mismatch rate (0 = identical, 1 = completely different)
- **Loci Used**: Number of SNPs with calls in both query and reference
- **Prob. Identity**: Probability two random varieties match by chance
- **Top Hits**: Best 3 matches with distances

### Hybrid Verification

| Metric | Good F1 | Selfed/Contaminated |
|--------|---------|---------------------|
| Concordance | > 95% | < 90% |
| Selfed-like | < 5% | > 10% |

### Seed Purity

- **Observed rate**: `n_offtypes / n_seeds`
- **95% CI**: Clopper-Pearson exact interval
- **Detection power**: Seeds needed to detect p% off-types at 95% confidence

---

## Troubleshooting

### Common Issues

#### "Reference panel not found"
```bash
# Build reference first
fingerprint build-ref --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz
```

#### "Thresholds not calibrated"
```bash
# Run calibration
fingerprint calibrate --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz
```

#### "Low call rate in query"
- Check DNA quality
- Ensure query uses same SNP set as panel
- Minimum 90% loci required

#### "No match found"
- Sample may be new variety not in reference
- Check for sample swap/contamination
- Try extended panel for close relatives

#### "Assay design failed"
- Check reference FASTA path and indexing
- Verify SNP coordinates match reference build
- Some SNPs may not have unique flanks

### Performance Tips

- Use `--cores N` with Snakemake for parallel execution
- For large datasets, increase memory: `snakemake --resources mem_mb=16000`
- Cache intermediate results: don't delete `data/processed/`

### Getting Help

```bash
# Command help
fingerprint --help
fingerprint match --help

# Check logs
snakemake --cores 4 2>&1 | tee log.txt
```

---

## File Formats

### Panel SNP List
```
SNP_001
SNP_002
SNP_003
```

### Scores Table (TSV)
```tsv
variant_id	chrom	pos	ref	alt	maf	pic	he	pi_locus	flank_uniqueness	flank_gc_content	nearby_variants_count	repeat_overlap	composite_score	call_rate	het_rate
```

### Thresholds (JSON)
```json
{
  "match_threshold": 0.015,
  "related_threshold": 0.045,
  "margin_to_second": 0.05,
  "min_loci_fraction": 0.9,
  "replicate_mean": 0.002,
  "replicate_std": 0.003,
  "nearest_mean": 0.12,
  "nearest_std": 0.04,
  "auc": 0.999,
  "eer": 0.001
}
```

### KASP Design Sheet (Excel)
Columns: Variant_ID, Chrom, Pos, Ref, Alt, Flank_Sequence, Common_Reverse, Allele1_Forward_FAM, Allele2_Forward_HEX, Common_Tm, Allele1_Tm, Allele2_Tm, Amplicon_Size, BLAST_Unique, No_Dimers, GC_OK, Tm_OK, Pass_All