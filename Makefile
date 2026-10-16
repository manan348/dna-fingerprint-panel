# DNA Fingerprint Panel - Makefile for common tasks

.PHONY: help install test lint format clean build docker-build docker-run run-pipeline

# Default target
help:
	@echo "DNA Fingerprint Panel - Available targets:"
	@echo ""
	@echo "  install        - Install package in development mode"
	@echo "  test           - Run unit tests"
	@echo "  test-cov       - Run tests with coverage"
	@echo "  lint           - Run linting (ruff)"
	@echo "  format         - Format code (ruff format)"
	@echo "  typecheck      - Run type checking (mypy)"
	@echo "  clean          - Clean build artifacts"
	@echo "  build          - Build package"
	@echo "  docker-build   - Build Docker image"
	@echo "  docker-run     - Run Docker container"
	@echo "  run-pipeline   - Run full Snakemake pipeline"
	@echo "  run-web        - Launch Streamlit web app"
	@echo "  pre-commit     - Install pre-commit hooks"

# Installation
install:
	pip install -e .[dev]

install-minimal:
	pip install -r requirements.txt

# Testing
test:
	pytest tests/ -v

test-cov:
	pytest tests/ -v --cov=fingerprint --cov-report=html --cov-report=term

test-slow:
	pytest tests/ -v -m slow

# Code quality
lint:
	ruff check src/ tests/

format:
	ruff format src/ tests/

typecheck:
	mypy src/

pre-commit:
	pre-commit install
	pre-commit run --all-files

# Cleaning
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete

# Building
build:
	python -m build

# Docker
docker-build:
	docker build -t fingerprint-panel:1.0.0 .

docker-run:
	docker run -it --rm -v $(PWD):/work fingerprint-panel:1.0.0

docker-shell:
	docker run -it --rm -v $(PWD):/work fingerprint-panel:1.0.0 bash

# Pipeline
run-pipeline:
	snakemake --cores 4 --use-conda

run-pipeline-dry:
	snakemake --cores 4 --use-conda -n

run-pipeline-force:
	snakemake --cores 4 --use-conda -F

# Web app
run-web:
	streamlit run app/streamlit_app.py

# CLI shortcuts
cli-score:
	fingerprint score --vcf data/processed/merged.vcf.gz --fasta data/reference/B73_v5.fa

cli-select:
	fingerprint select --scores results/tables/candidate_snp_scores.tsv --vcf data/processed/filtered.vcf.gz

cli-robustness:
	fingerprint robustness --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz

cli-calibrate:
	fingerprint calibrate --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz

cli-build-ref:
	fingerprint build-ref --core-panel results/panels/core_panel_snps.txt --vcf data/processed/filtered.vcf.gz

cli-design:
	fingerprint design-assays --panel results/panels/core_panel_details.tsv --fasta data/reference/B73_v5.fa

# Data download (manual - placeholder)
download-data:
	@echo "Download datasets manually to data/raw/ and record in data/MANIFEST.tsv"
	@echo "See docs/protocol.md for dataset list"

# Documentation
docs-build:
	cd docs && sphinx-build -b html . _build

docs-serve:
	cd docs && python -m http.server -d _build 8000

# Release
release-check:
	@echo "Checklist for release:"
	@echo "  1. Update version in pyproject.toml and __init__.py"
	@echo "  2. Update CHANGELOG.md"
	@echo "  3. Run: make test-cov"
	@echo "  4. Run: make lint && make typecheck"
	@echo "  5. Run: make build"
	@echo "  6. Tag: git tag v\$(VERSION)"
	@echo "  7. Push: git push origin --tags"

# Development setup
dev-setup: install pre-commit
	@echo "Development environment ready!"

# CI simulation
ci: lint typecheck test-cov
	@echo "CI checks passed!"