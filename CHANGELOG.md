# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-15

### Added
- Complete pipeline for DNA fingerprint panel design and validation
- SNP scoring module (MAF, PIC, HE, PI, flank uniqueness, GC content)
- Lazy-greedy set multicover panel selection with ILP refinement
- Robustness simulation (genotyping error, missing data, mixtures, assay dropout)
- Matching engine with calibrated decision thresholds (MATCH/RELATED/NO MATCH)
- Hybrid verification module
- Seed-lot purity statistics (Clopper-Pearson CI)
- KASP assay design with Primer3 and in silico validation
- Command-line interface with all pipeline steps
- Streamlit web interface for interactive use
- Snakemake workflow for reproducible pipeline execution
- Docker container for reproducible deployment
- Comprehensive test suite
- Documentation (protocol, user guide, lab SOP)

### Changed
- N/A (initial release)

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

---

## [Unreleased]

### Planned for v1.1.0
- Wheat (polyploid) support with PolyMarker integration
- Targeted amplicon panel design (DArTag)
- Haplotype-based panel selection
- REST API for matching service
- Integration with breeding databases
- Multi-panel comparison dashboard

### Planned for v1.2.0
- GPU-accelerated ILP solving
- Active learning for panel augmentation
- Federated learning across institutions
- Mobile app for field-based verification