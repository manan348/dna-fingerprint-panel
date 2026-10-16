# DNA Fingerprint Panel for Varietal Identification

A low-compute SNP panel design, validation, and deployment toolkit for seed purity checks and variety certification. Primary crop: maize (extension to wheat described).

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.28+-red.svg)](https://streamlit.io/)
[![Snakemake](https://img.shields.io/badge/Snakemake-7.32+-orange.svg)](https://snakemake.github.io/)
[![Tests](https://img.shields.io/badge/Tests-43%20passing-brightgreen.svg)](tests/)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [System Architecture](#system-architecture)
- [Tech Stack](#tech-stack)
- [Project Structure](#project-structure)
- [Local Development](#local-development)
- [Environment Variables](#environment-variables)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Outputs](#outputs)
- [License](#license)

---

## Overview

Seed is the delivery vehicle for genetic gain. A farmer who buys a bag labelled with a high-yielding hybrid is trusting that the seed inside really is that hybrid, genetically pure and correctly labelled. In practice this trust is hard to verify. Morphology-based identification and grow-out tests take a full season, depend on the environment, require expert judgement and often cannot separate closely related lines. Mislabelled, mixed or counterfeit seed lots reduce farm yields, damage breeders' reputations and undermine certification systems.

DNA fingerprinting replaces subjective, slow identification with objective, season-independent evidence. A small panel of well-chosen SNP markers gives every variety a unique genetic barcode. Public genotype resources are now large enough that the panel can be designed entirely in silico on a laptop, because the computation is mostly filtering, counting and set optimisation.

This project builds a complete, reproducible workflow that:

1. Selects a minimal SNP panel able to separate every variety in a target collection with a defined safety margin against genotyping error
2. Converts that panel into ready-to-order KASP assays
3. Provides a validated decision tool for match, no-match and mixture calls

Designed around a computationally light core (laptop with 16 GB RAM, no GPU needed), with a clear path to wet-lab validation and institutional use.

---

## Features

| Feature | Description |
|---------|-------------|
| **SNP Scoring** | MAF, PIC, expected heterozygosity, probability of identity, flank uniqueness, GC content, nearby variant density |
| **Panel Selection** | Lazy-greedy set multicover with ILP (OR-Tools CP-SAT) refinement; physical spacing and per-chromosome quotas |
| **Robustness Testing** | Simulations for genotyping error, missing data, seed mixtures, and assay dropout |
| **Matching Engine** | Mismatch-rate distance with calibrated thresholds (MATCH / AMBIGUOUS / NO MATCH / INSUFFICIENT DATA) |
| **Hybrid Verification** | F1 concordance and selfed-like fraction against declared parents |
| **Seed Purity** | Off-type rate with Clopper-Pearson confidence intervals |
| **KASP Assay Design** | Primer3-based common reverse + allele-specific forward primers with FAM/HEX tails; BLAST uniqueness and dimer checks |
| **Web Interface** | Streamlit app for interactive matching, hybrid verification, purity testing, and assay design |
| **Reproducible Pipeline** | Snakemake workflow with pinned conda environment and Docker container |
| **Unit Tests** | 43 passing tests covering scoring, selection, matching, and robustness |

---

## System Architecture

```text
Browser
   │  HTTP
   ▼
Streamlit App  (app/streamlit_app.py)
   │  imports
   ▼
fingerprint package (src/fingerprint/)
   ├── io.py          VCF/PLINK/Parquet I/O, sample harmonisation
   ├── qc.py          Sample/SNP QC, clone detection
   ├── scoring.py     MAF, PIC, HE, PI, flank scores
   ├── selection.py   Greedy + ILP panel selection
   ├── robustness.py  Error/missing/mixture/dropout simulations
   ├── matcher.py     Distance, thresholds, hybrid verification, purity
   ├── assay.py       KASP primer design (Primer3)
   └── cli.py         Command-line interface
   │
   ▼
Reference Database (Parquet / SQLite)
   │  varieties × SNPs versioned with panel
   ▼
Raw & Processed Data (data/raw/, data/processed/)
   │  VCF → PLINK2 → filtered genotypes
   ▼
Snakemake Workflow (Snakefile)
   │  reproducible, versioned, scripted
   ▼
Conda Environment + Docker
   │  pinned dependencies, exact reproduction
```

The system is a five-layer pipeline. Each layer reads versioned files from the previous layer and writes versioned files for the next, so any stage can be re-run, audited or replaced independently.

---

## Tech Stack

### Frontend

| Technology | Version | Purpose |
|------------|---------|---------|
| Streamlit | 1.28+ | Web UI for matching, hybrid verification, purity testing, assay design |

### Backend

| Technology | Version | Purpose |
|------------|---------|---------|
| Python | 3.10+ | Core logic and CLI |
| NumPy | 1.24+ | Array operations |
| pandas | 2.0+ | Data frames and I/O |
| SciPy | 1.10+ | Statistics and distributions |
| scikit-allel | 1.3+ | Genotype array handling |
| PyVCF3 | 0.10+ | VCF reading |
| Primer3-Py | 2.0+ | Primer design |
| Biopython | 1.81+ | Sequence utilities |
| OR-Tools | 9.7+ | CP-SAT solver for ILP panel selection |
| Plotly / Matplotlib / Seaborn | 5.17+ / 3.7+ / 0.12+ | Figures and plots |
| PyYAML | 6.0+ | Configuration |
| Click | 8.1+ | CLI framework |
| pytest | 7.4+ | Unit testing |

### Bioinformatics Tools (conda / Docker)

| Tool | Version | Purpose |
|------|---------|---------|
| bcftools | 1.18 | VCF manipulation and indexing |
| samtools | 1.18 | FASTA indexing and region extraction |
| bedtools | 2.31 | Region intersections |
| plink2 | 2.00 | Genotype filtering and conversion |
| CrossMap | 0.6 | Coordinate liftover |
| Snakemake | 7.32 | Reproducible workflow |

### Infrastructure

| Component | Details |
|-----------|---------|
| Packaging | `pyproject.toml`, `pip install -e .` |
| Environment | `environment.yml` (conda), `requirements.txt` (pip) |
| Container | `Dockerfile` for exact reproduction |
| Workflow | `Snakefile` (Snakemake) |
| Linting | `ruff`, `black`, `mypy` |
| License | MIT |

---

## Project Structure

```text
fingerprint-panel/
├── README.md                     # This file
├── LICENSE                       # MIT license
├── CHANGELOG.md                  # Version history
├── Makefile                      # Common tasks (test, lint, build, run)
├── pyproject.toml                # Package config, pytest, ruff, mypy
├── environment.yml               # Conda environment (pinned)
├── requirements.txt              # pip dependencies
├── Dockerfile                    # Reproducible container
├── Snakefile                     # Snakemake workflow
├── config.yaml                   # All pipeline parameters (frozen)
├── .gitignore                    # Git ignore rules
├── .pre-commit-config.yaml       # Pre-commit hooks
│
├── app/                          # Streamlit web interface
│   └── streamlit_app.py          # Match, verify hybrid, purity, design assays
│
├── docs/                         # Documentation
│   ├── protocol.md               # Frozen analysis plan (hypotheses, experiments, gates)
│   ├── user_guide.md             # End-user guide
│   └── lab_sop_kasp.md           # Wet-lab KASP genotyping SOP
│
├── scripts/                      # Utility scripts
│   └── generate_demo_data.py     # Generates demo data for the web app
│
├── src/fingerprint/              # Core Python package
│   ├── __init__.py
│   ├── cli.py                    # Command-line interface
│   ├── io.py                     # VCF/PLINK/Parquet I/O, harmonisation
│   ├── qc.py                     # Sample/SNP QC, clone detection
│   ├── scoring.py                # MAF, PIC, HE, PI, flank scores
│   ├── selection.py              # Greedy + ILP panel selection
│   ├── robustness.py             # Error/missing/mixture/dropout simulations
│   ├── matcher.py                # Distance, thresholds, hybrid, purity
│   └── assay.py                  # KASP primer design (Primer3)
│
├── tests/                        # Unit tests (43 passing)
│   ├── conftest.py               # Pytest fixtures
│   ├── test_scoring.py           # SNP scoring tests
│   ├── test_selection.py         # Panel selection tests
│   ├── test_matcher.py           # Matching engine tests
│   └── test_robustness.py        # Robustness simulation tests
│
├── workflows/                    # Snakemake rule scripts
│   └── scripts/
│       ├── harmonise_samples.py
│       ├── clone_detection.py
│       ├── score_snps.py
│       ├── select_panel.py
│       ├── robustness_test.py
│       ├── calibrate_thresholds.py
│       ├── design_kasp.py
│       ├── plot_resolution.py
│       └── plot_panel_map.py
│
├── data/                         # Data directories (not committed)
│   ├── raw/                      # Downloads (not committed)
│   ├── processed/                # Filtered genotypes
│   ├── demo/                     # Demo data for web app
│   ├── MANIFEST.tsv              # Source, licence, checksum
│   └── metadata.tsv              # Sample metadata
│
├── results/                      # Pipeline outputs (versioned)
│   ├── panels/                   # Core, backup, extended SNP lists
│   ├── tables/                   # Scores, resolution curve, robustness, thresholds
│   ├── figures/                  # Resolution curve, panel map, distances, heatmap
│   └── kasp_design_sheet.xlsx    # KASP primer designs
│
└── notebooks/                    # Exploratory analysis and figures
```

---

## Local Development

### Prerequisites

- **OS:** Linux, macOS, or Windows (WSL2 recommended for full pipeline)
- **Python:** 3.10+
- **Conda:** Miniconda or Anaconda
- **Docker:** Optional, for containerised execution
- **Git LFS:** Required for large binary files (BAM, etc.)

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/dna-fingerprint-panel.git
cd dna-fingerprint-panel
```

### 2. Create the Conda Environment

```bash
conda env create -f environment.yml
conda activate fingerprint-panel
```

> **Note:** On Windows, the environment path contains a space (`C:\Users\Abdul Manan\anaconda3\envs\fingerprint-panel`). Activate the environment before running any executables.

### 3. Install the Package

```bash
pip install -e .
```

> `pysam` fails to build on Windows without a C compiler. The code handles missing `pysam` gracefully; full VCF/Tabix functionality requires conda/Docker/WSL2.

### 4. Run the Test Suite

```bash
pytest tests/ -v
```

**Expected:** 43 passed, 1 error (samtools not available on Windows — expected).

### 5. Run the Pipeline

```bash
# Dry run
snakemake --cores 4 --use-conda -n

# Full pipeline
snakemake --cores 4 --use-conda
```

### 6. Launch the Web App

```bash
streamlit run app/streamlit_app.py
```

Open **http://localhost:8501** in your browser.

### 7. Build the Docker Container (optional)

```bash
docker build -t fingerprint-panel:1.0.0 .
docker run -it --rm -v $(pwd):/work fingerprint-panel:1.0.0
```

---

## Environment Variables

No environment variables are required. All pipeline parameters are in `config.yaml`.

Streamlit reads the following from the UI sidebar:
- **Reference Database:** Path to Parquet reference panel (default: `results/reference_panel.parquet`)
- **Thresholds JSON:** Path to calibrated thresholds (default: `results/tables/calibrated_thresholds.json`)

---

## Usage

### Command-Line Interface

```bash
# Show all commands
fingerprint --help

# Step 4: Score SNPs for informativeness and assay-friendliness
fingerprint score \
    --vcf data/processed/merged.vcf.gz \
    --fasta data/reference/B73_v5.fa

# Step 5: Select fingerprint panel (greedy + ILP)
fingerprint select \
    --scores results/tables/candidate_snp_scores.tsv \
    --vcf data/processed/filtered.vcf.gz

# Step 7: Robustness simulations
fingerprint robustness \
    --core-panel results/panels/core_panel_snps.txt \
    --vcf data/processed/filtered.vcf.gz

# Step 8: Calibrate matching thresholds
fingerprint calibrate \
    --core-panel results/panels/core_panel_snps.txt \
    --vcf data/processed/filtered.vcf.gz

# Step 8: Build reference database for matching
fingerprint build-ref \
    --core-panel results/panels/core_panel_snps.txt \
    --vcf data/processed/filtered.vcf.gz

# Step 8: Match an unknown sample
fingerprint match \
    --reference results/reference_panel.parquet \
    --thresholds results/tables/calibrated_thresholds.json \
    --query query_sample.csv

# Step 8: Verify an F1 hybrid
fingerprint verify-hybrid \
    --hybrid hybrid.csv \
    --parent1 parent1.csv \
    --parent2 parent2.csv

# Step 8: Seed-lot purity statistics
fingerprint seed-purity \
    --n-seeds 100 \
    --n-offtypes 2

# Step 9: Design KASP assays
fingerprint design-assays \
    --panel results/panels/core_panel_details.tsv \
    --fasta data/reference/B73_v5.fa
```

### Web Interface

```bash
streamlit run app/streamlit_app.py
```

**Tabs:**
1. **Match Sample** — Upload a genotype CSV, get a variety match decision
2. **Verify Hybrid** — Upload hybrid and parent genotypes, verify F1
3. **Seed Purity** — Calculate off-type rate and confidence intervals
4. **Design Assays** — Design KASP primers for a SNP panel
5. **Reference Panel** — Browse reference varieties and SNPs

### Demo Data

Try the app immediately with built-in demo data:

```bash
# Generate demo files
python scripts/generate_demo_data.py

# Launch app
streamlit run app/streamlit_app.py

# In the sidebar, click "Load Demo Data"
```

Demo data includes 10 varieties × 20 SNPs, a query sample matching V01, hybrid verification files, and calibrated thresholds.

---

## Pipeline Steps

| Step | Description | Command |
|------|-------------|---------|
| 1 | Define scope, target varieties, success criteria | `docs/protocol.md` (manual) |
| 2 | Acquire, organise, harmonise data | Snakemake / manual download |
| 3 | Quality control and filtering | Snakemake `qc_filter` |
| 4 | Score SNP informativeness | `fingerprint score` |
| 5 | Select panel (greedy + ILP) | `fingerprint select` |
| 6 | Experimental design | Documented in protocol |
| 7 | Robustness testing | `fingerprint robustness` |
| 8 | Calibrate thresholds | `fingerprint calibrate` |
| 9 | Design KASP assays | `fingerprint design-assays` |
| 10 | Package and release | Manual |

---

## Outputs

| Path | Description |
|------|-------------|
| `results/panels/core_panel_snps.txt` | Core fingerprint panel SNP IDs (48–96 SNPs) |
| `results/panels/backup_panel_snps.txt` | Backup panel (~20% extra) |
| `results/panels/extended_panel_snps.txt` | Extended panel for close relatives (192–384 SNPs) |
| `results/tables/candidate_snp_scores.tsv` | All SNP informativeness scores |
| `results/tables/panel_resolution_curve.tsv` | Panel size vs. discrimination |
| `results/tables/robustness_results.tsv` | Simulation results |
| `results/tables/calibrated_thresholds.json` | Match/related/no-match thresholds |
| `results/reference_panel.parquet` | Reference database for matching |
| `results/kasp_design_sheet.xlsx` | KASP primer designs |
| `docs/lab_sop_kasp.md` | Wet-lab genotyping protocol |
| `results/figures/*.png` | Resolution curve, panel map, distance distributions, heatmap |

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

**Prepared by:** Abdul Manan, Founder & CEO, NeuroCrop  
**Affiliation:** Plant Breeding & Genetics, University of the Punjab, Lahore, Pakistan  
**Version:** 1.0 | October 2026  
**Primary crop:** Maize (extension to wheat described)