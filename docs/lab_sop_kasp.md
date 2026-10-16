# Lab SOP: KASP Genotyping for Varietal Identification Panel

**Document Version:** 1.0
**Effective Date:** October 2026
**Panel:** Core Fingerprint Panel (96 SNPs)
**Prepared by:** NeuroCrop

---

## 1. Purpose

This SOP describes the laboratory procedure for genotyping samples using the DNA Fingerprint Panel KASP assays for varietal identification and seed purity assessment.

---

## 2. Scope

Applies to all genotyping of maize varieties using the 96-SNP core fingerprint panel (plus backup assays). Covers reaction setup, thermal cycling, data analysis, and quality control.

---

## 3. Materials and Equipment

### 3.1 Reagents

| Reagent | Supplier | Catalogue No. | Storage |
|---------|----------|---------------|---------|
| KASP Master Mix (2×) | LGC Biosearch | KBS-1001 | -20°C |
| Assay primers (custom) | IDT / LGC | See design sheet | -20°C |
| Nuclease-free water | - | - | RT |
| DNA samples (20-50 ng/µL) | - | - | -20°C or 4°C |

### 3.2 Equipment

| Equipment | Specification |
|-----------|---------------|
| qPCR instrument | FAM/HEX channels, endpoint read (e.g., Bio-Rad CFX, QuantStudio, LightCycler) |
| Pipettes | 0.5-10 µL, 2-20 µL, 20-200 µL, calibrated |
| Plate centrifuge | 96/384-well, ≥1000 × g |
| Vortex mixer | - |
| Cold block / ice | - |

### 3.3 Consumables

| Item | Specification |
|------|---------------|
| PCR plates | 96-well or 384-well, skirted, compatible with instrument |
| Optical seals | Pressure-sensitive or heat seal |
| Tips | Filtered, low-retention |

---

## 4. Assay Information

### 4.1 Panel Composition

- **Core panel:** 96 SNPs
- **Backup panel:** 19 SNPs (20%)
- **Total assays:** 115

### 4.2 Primer Concentrations

Prepare assay mix for each SNP:
- Allele-specific forward (FAM): 12 µM
- Allele-specific forward (HEX): 12 µM
- Common reverse: 30 µM
- In TE buffer (10 mM Tris, 0.1 mM EDTA, pH 8.0)

Store at -20°C. Avoid freeze-thaw cycles.

### 4.3 Plate Layout (96-well)

| Well | Content |
|------|---------|
| A1-H1 | Sample 1 (96 assays) |
| A2-H2 | Sample 2 (96 assays) |
| ... | ... |
| A11-H11 | Sample 11 |
| A12 | NTC (no template control) |
| B12 | Positive control (variety A) |
| C12 | Positive control (variety B) |
| D12 | Replicate of sample 1 |
| E12 | Replicate of sample 2 |
| F12-H12 | Empty / additional controls |

*Adjust for 384-well format accordingly.*

---

## 5. Procedure

### 5.1 DNA Preparation

1. Quantify DNA (Qubit dsDNA HS or NanoDrop).
2. Normalize to **20-50 ng/µL** in TE or nuclease-free water.
3. Minimum 5 µL per sample needed.
4. Keep on ice during setup.

### 5.2 Reaction Setup

**Per reaction (5 µL total):**

| Component | Volume | Final Conc. |
|-----------|--------|-------------|
| 2× KASP Master Mix | 2.5 µL | 1× |
| Assay mix | 0.07 µL | See below |
| DNA (20-50 ng/µL) | 2.0 µL | 8-20 ng |
| Nuclease-free water | 0.43 µL | - |
| **Total** | **5.0 µL** | - |

**Assay mix final concentrations in reaction:**
- Each allele-specific primer: ~0.17 µM
- Common reverse: ~0.42 µM

#### Setup Steps

1. Thaw master mix and assay mixes on ice. Vortex briefly, centrifuge.
2. Prepare master mix for N reactions + 10% overage:
   ```
   Master mix = N × 1.1 × (2.5 µL KASP MM + 0.07 µL assay + 0.43 µL water)
   ```
3. Aliquot **3.0 µL** master mix per well.
4. Add **2.0 µL** DNA sample per well.
5. For NTC: add 2.0 µL water instead of DNA.
6. Seal plate securely.
7. Centrifuge at 1000 × g for 30 seconds.

### 5.3 Thermal Cycling Protocol

| Step | Temperature | Time | Cycles |
|------|-------------|------|--------|
| Enzyme activation | 94°C | 15 min | 1 |
| Touchdown | 94°C | 20 sec | 10 |
| | 61°C (-0.6°C/cycle) | 60 sec | |
| Amplification | 94°C | 20 sec | 26 |
| | 55°C | 60 sec | |
| Hold | 4°C | ∞ | 1 |

**Total run time:** ~1 hour 45 minutes

#### Notes:
- Use "fast" ramp rates if available
- Endpoint read only (no real-time monitoring needed)
- Read FAM and HEX channels at end

---

## 6. Data Analysis

### 6.1 Cluster Calling

1. Open run file in instrument software (e.g., CFX Maestro, QuantStudio Design & Analysis).
2. Select endpoint analysis / genotype calling mode.
3. Auto-call clusters for each assay.
4. **Manual review required** for:
   - Assays with >3 clusters
   - Compressed clusters
   - Samples between clusters
   - NTC showing amplification

### 6.2 Export Genotypes

Export as CSV with columns:
```csv
Sample,Assay,Allele_Call
Sample_01,SNP_001,AA
Sample_01,SNP_002,BB
Sample_01,SNP_003,AB
```

Allele calls: AA (hom FAM), BB (hom HEX), AB (het), NN (no call)

### 6.3 Convert to Analysis Format

Convert to fingerprint toolkit format:
```python
# AA -> 0, AB -> 1, BB -> 2, NN -> -1
mapping = {'AA': 0, 'AB': 1, 'BB': 2, 'NN': -1}
```

Save as CSV: `snp_id,genotype`

---

## 7. Quality Control

### 7.1 Per-Assay QC

| Metric | Threshold | Action if Failed |
|--------|-----------|------------------|
| Call rate | ≥ 95% | Repeat assay; replace with backup |
| Cluster separation | Clear, distinct | Review manually; repeat if needed |
| NTC amplification | None | Investigate contamination; repeat plate |

### 7.2 Per-Sample QC

| Metric | Threshold | Action if Failed |
|--------|-----------|------------------|
| Call rate | ≥ 90% | Re-genotype sample |
| Replicate concordance | ≥ 99% | Investigate; re-genotype |

### 7.3 Plate-Level QC

- All positive controls match expected genotypes
- NTC shows no amplification in any assay
- Replicate samples cluster together
- No systematic cluster shifts across plate

### 7.4 QC Report

Document for each plate:
- Plate ID, date, operator
- Number of samples, assays
- Per-assay call rates
- Per-sample call rates
- Replicate concordance
- Any deviations from protocol
- Pass/Fail decision

---

## 8. Troubleshooting

| Problem | Possible Cause | Solution |
|---------|----------------|----------|
| Poor clustering (all assays) | Bad DNA / inhibitors | Re-extract DNA; dilute 1:10 |
| Poor clustering (specific assay) | Primer issue | Replace with backup assay |
| Low call rate | Low DNA concentration | Increase DNA input; re-quantify |
| NTC amplification | Contamination | Clean work area; new reagents |
| Compressed clusters | Suboptimal cycling | Verify thermal profile; check instrument calibration |
| Unexpected genotype | Sample swap | Verify sample tracking; re-run |

---

## 9. Safety

- Follow institutional biosafety guidelines
- KASP master mix: non-hazardous, but wear gloves
- Dispose of PCR waste per local regulations
- UV-decontaminate work surfaces between runs

---

## 10. Documentation and Records

Retain for minimum 5 years:
- Plate setup sheets
- Raw instrument files
- Called genotype files
- QC reports
- Deviation logs

---

## 11. References

1. Semagn et al. (2014). KASP genotyping for crop improvement. *Mol Breeding* 33:1-14.
2. LGC Biosearch KASP Genotyping Chemistry User Guide.
3. Primer3 manual (Untergasser et al. 2012).

---

## 12. Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | Oct 2026 | NeuroCrop | Initial release |

---

## Appendix: Assay Design Sheet

See accompanying `results/kasp_design_sheet.xlsx` for complete primer sequences, Tm, amplicon sizes, and QC status for all 115 assays.

### Example Entry

| Variant_ID | Chrom | Pos | Ref | Alt | Common_Reverse | Allele1_Forward_FAM | Allele2_Forward_HEX | Pass_All |
|------------|-------|-----|-----|-----|----------------|---------------------|---------------------|----------|
| SNP_001 | chr1 | 1234567 | A | G | TGCATGCATGCATGCATGCA | GAAGGTGACCAAGTTCATGCTACGTACGTACG | GAAGGTCGGAGTCAACGGATTTGCATGCATGC | ✅ |

---

**End of SOP**