# Analysis Protocol: DNA Fingerprint Panel for Varietal Identification

**Version:** 1.0
**Date:** October 2026
**Status:** FROZEN - Do not modify after Step 1 completion

---

## 1. Hypotheses

### H1: Panel Size vs. Discrimination
A panel of 48-96 SNPs selected by greedy multicover will separate ≥99% of non-clonal variety pairs at k=3, outperforming PIC-ranked and random panels of the same size.

### H2: Robustness to Errors
The selected panel will maintain ≥95% correct identification rate under:
- Genotyping error rate ≤ 2%
- Missing data rate ≤ 10%
- Seed mixture proportion ≤ 10%
- Assay dropout ≤ 10%

### H3: Threshold Calibration
Match thresholds calibrated from replicate data will achieve:
- False match rate < 0.1%
- False non-match rate < 5%
- AUC > 0.99

### H4: Assay Design Success
≥80% of selected SNPs will yield working KASP assays with:
- Call rate ≥ 95%
- Replicate concordance ≥ 99%

---

## 2. Experiments

### E1: Reference Collection Assembly
**Objective:** Build harmonised genotype matrix from public datasets.
**Data:** Maize HapMap3, NAM founders, 600K array, CIMMYT lines.
**Metrics:** N varieties, N SNPs, call rate, clone groups detected.
**Acceptance:** ≥500 unique varieties, ≥10,000 candidate SNPs after QC.

### E2: SNP Scoring and Filtering
**Objective:** Rank all candidate SNPs by composite score.
**Method:** MAF (≥0.20), PIC, HE, flank uniqueness, GC content, nearby variants.
**Metrics:** Score distributions, filter pass rates.
**Acceptance:** ≥5,000 SNPs pass all filters.

### E3: Panel Selection
**Objective:** Find minimal panel separating all pairs at k=3.
**Method:** Lazy-greedy multicover → ILP refinement.
**Constraints:** 2 Mb spacing, ≤15% per chromosome.
**Metrics:** Panel size, unresolved pairs, resolution curve vs. baselines.
**Acceptance:** ≤96 SNPs, 0 unresolved non-clonal pairs.

### E4: Robustness Simulation
**Objective:** Quantify panel performance under perturbations.
**Design:** Full factorial grid:
- Error rates: 0%, 1%, 2%, 5%
- Missing rates: 0%, 5%, 10%, 15%
- Mixture proportions: 0%, 5%, 10%, 20%
- Dropout rates: 0%, 5%, 10%, 15%
**Replicates:** 200 per condition.
**Metrics:** Accuracy, mean/min separation, backup recovery.
**Acceptance:** ≥95% accuracy at 2% error, 10% missing, 10% mixture, 10% dropout.

### E5: Threshold Calibration
**Objective:** Set match/related thresholds from replicate data.
**Method:** Compute replicate distances and nearest-neighbor distances.
**Thresholds:** Match at 99.9th percentile of replicate distance; related at 3× match threshold.
**Metrics:** ROC AUC, EER, FMR, FNMR.
**Acceptance:** AUC > 0.99, EER < 0.01, clear gap between replicate and nearest distributions.

### E6: Matching Engine Validation
**Objective:** Test matcher on hold-out varieties.
**Design:** 10-fold cross-validation by variety.
**Metrics:** Per-fold accuracy, confusion matrix.
**Acceptance:** Mean accuracy ≥ 98%.

### E7: Assay Design
**Objective:** Design KASP primers for core + backup panels.
**Method:** Primer3 for common reverse; allele-specific forwards with FAM/HEX tails.
**Checks:** BLAST uniqueness, dimer/hairpin, GC, Tm.
**Acceptance:** ≥80% pass all in silico checks.

### E8: Wet-Lab Pilot
**Objective:** Validate assays on real DNA.
**Design:** 48 lines (including 3 relative pairs, 3 replicates, 1 NTC) on 96-well plates.
**Metrics:** Per-assay call rate, replicate concordance, cluster quality.
**Acceptance:** Mean call rate ≥ 95%, mean concordance ≥ 99%, 0 NTC amplification.

---

## 3. Statistical Analysis Plan

### Primary Comparisons
- Designed panel vs. PIC-ranked panel (paired, same data splits)
- Designed panel vs. random panels (1,000 random draws per size)

### Uncertainty Quantification
- Bootstrap over varieties (1,000 resamples) for all metrics
- Report 95% percentile intervals

### Hypothesis Tests
- Paired t-test / Wilcoxon for method comparisons
- Effect sizes (Cohen's d) reported with p-values
- Significance threshold: α = 0.05 (no multiple testing correction for exploratory analyses)

### Proportions
- Clopper-Pearson exact CIs for call rates, concordance, off-type rates

---

## 4. Acceptance Criteria (Go/No-Go Gates)

| Gate | Criterion | Threshold |
|------|-----------|-----------|
| G1 (Week 4) | Unique varieties after QC | ≥ 500 |
| G1 | Candidate SNPs after filters | ≥ 5,000 |
| G2 (Week 8) | Core panel size | ≤ 96 |
| G2 | Unresolved non-clonal pairs | 0 |
| G2 | Resolution at k=3 (96 SNPs) | ≥ 99% |
| G2 | Robustness accuracy (2% err, 10% miss) | ≥ 95% |
| G3 (Week 11) | Matcher AUC | > 0.99 |
| G3 | Assays passing in silico checks | ≥ 80% |
| G4 (Week 14) | Wet-lab call rate | ≥ 95% |
| G4 | Wet-lab replicate concordance | ≥ 99% |

---

## 5. Data Versioning

All data files carry version tags:
- `data/MANIFEST.tsv` - Source datasets, URLs, checksums, licences
- `data/metadata.tsv` - Sample metadata with schema validation
- `results/panels/*_v1.txt` - Panel files with git commit hash
- `config.yaml` - All parameters frozen at protocol freeze

---

## 6. Deviations Log

| Date | Step | Deviation | Reason | Author |
|------|------|-----------|--------|--------|
| | | | | |

*Any deviation from this protocol must be logged here with justification.*