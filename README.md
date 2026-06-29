# Fragmentomics NF

![Example prediction analysis report](docs/assets/example_prediction_analysis_report.png)

Production-style Nextflow DSL2 pipeline for cfDNA WGS fragmentomics research.

The primary pipeline lives in [`fragmentomics_nf/`](fragmentomics_nf/), with full usage notes in [`fragmentomics_nf/README.md`](fragmentomics_nf/README.md).

## Architecture

![Architecture workflow](fragmentomics_nf/docs/assets/architecture_workflow_visual.png)

- [Download the architecture workflow PDF](fragmentomics_nf/docs/assets/architecture_workflow_visual.pdf)
- [View the architecture workflow documentation](fragmentomics_nf/docs/architecture_workflow.md)
- [View the full traceability matrix representation](fragmentomics_nf/docs/traceability_matrix.md)

## What This Repository Is

This repository contains a research-oriented Nextflow DSL2 pipeline for cfDNA WGS fragmentomics. It turns sequencing-derived inputs or public fragment files into normalized fragmentomic features, trains transparent baseline models, generates prediction scores, and produces reviewable HTML reports.

The workflow is designed around reproducibility and inspectability:

- accepts FASTQ, BAM, CRAM, or public fragment-file inputs
- validates sample metadata before analysis
- filters and normalizes fragment records
- builds genomic bins and corrects bin-level fragment counts for GC bias
- assembles model-ready feature matrices
- trains and applies baseline scikit-learn models
- summarizes model behavior with metrics, plots, tables, and a rendered report

Outputs are for research use only. Model scores are fragmentome score estimates, not clinical diagnostic calls.

## Example Report

The repository includes a small rendered example report so the output format is visible without committing raw data or heavyweight pipeline artifacts:

- [`fragmentomics_nf/results/finaledb_lung_demo_300/analysis/prediction_analysis_report.html`](fragmentomics_nf/results/finaledb_lung_demo_300/analysis/prediction_analysis_report.html)
- [`fragmentomics_nf/results/finaledb_lung_demo_300/analysis/prediction_analysis_summary.json`](fragmentomics_nf/results/finaledb_lung_demo_300/analysis/prediction_analysis_summary.json)

Large generated files remain ignored, including Nextflow `work/`, bulk `results/` outputs, raw fragments, reference genomes, feature matrices, model binaries, and compressed genomics artifacts.
