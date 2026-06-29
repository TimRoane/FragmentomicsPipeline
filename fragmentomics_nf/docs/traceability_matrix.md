# Traceability Matrix

## Document Purpose

This document provides a representative traceability matrix for the `fragmentomics_nf` pipeline. It maps user needs and system requirements to design elements, implementation artifacts, expected outputs, and verification evidence.

The matrix is intended for research software documentation, technical review, onboarding, and presentation. It is not a finalized regulatory design history file, but it is structured so it can be expanded into one.

## Scope

In scope:

- Sample sheet and manifest validation.
- FASTQ/BAM/CRAM preprocessing.
- FinaleDB-style fragment ingest.
- Fragment-level QC and filtering.
- Genomic bin generation and GC annotation.
- Fragment counting, GC correction, and feature matrix assembly.
- Baseline model training, prediction, and post-analysis reporting.
- Research-only language and output interpretation.

Out of scope:

- Clinical diagnostic claims.
- Locked clinical thresholds.
- Laboratory wet-lab sample handling.
- Production deployment qualification.
- External cohort validation.

## Traceability ID Scheme

| Prefix | Meaning |
| --- | --- |
| `UN` | User need |
| `SR` | System requirement |
| `DR` | Design requirement |
| `IR` | Implementation reference |
| `VER` | Verification evidence |
| `OUT` | Output artifact |

## User Needs

| ID | User Need | Rationale | Related System Requirements |
| --- | --- | --- | --- |
| `UN-001` | The pipeline shall process cfDNA WGS-derived inputs into fragmentomics features. | Enables reproducible fragmentomics analysis from sequencing-derived data. | `SR-001`, `SR-002`, `SR-003`, `SR-004`, `SR-005` |
| `UN-002` | The pipeline shall support public fragment files in addition to conventional FASTQ/BAM/CRAM inputs. | Enables research workflows using public datasets and precomputed fragments. | `SR-006`, `SR-007`, `SR-008` |
| `UN-003` | The pipeline shall produce interpretable tabular outputs suitable for downstream modeling and review. | Users need inspectable feature matrices, QC files, and summaries. | `SR-009`, `SR-010`, `SR-011` |
| `UN-004` | The pipeline shall train and apply a transparent baseline research model. | Enables demonstration, benchmarking, and exploratory classification. | `SR-012`, `SR-013`, `SR-014` |
| `UN-005` | The pipeline shall generate human-readable prediction analysis reports. | Users need presentation-ready and reviewable model performance summaries. | `SR-015`, `SR-016`, `SR-017` |
| `UN-006` | The pipeline shall avoid clinical diagnostic language in research outputs. | Reduces the chance that model scores are misinterpreted as diagnoses. | `SR-018` |
| `UN-007` | The pipeline shall keep heavyweight raw data and runtime artifacts separate from source control. | Prevents accidental publication of raw data, reference genomes, and large generated files. | `SR-019` |

## System Requirements Traceability

| Requirement ID | Requirement | Design / Workflow Element | Implementation Reference | Output Artifact(s) | Verification Evidence | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `SR-001` | The pipeline shall validate input sample metadata before analysis. | `VALIDATE_SAMPLESHEET` process normalizes input rows and emits a validated sample sheet plus run manifest. | `main.nf`; `modules/local/validate_samplesheet.nf`; `bin/validate_samplesheet.py` | `OUT-001`: `results/<run>/validated_samplesheet.csv`; `OUT-002`: `results/<run>/run_manifest.json` | `VER-001`: `tests/test_finaledb_ingest.py::test_validate_samplesheet_accepts_finaledb_manifest_columns` | Implemented |
| `SR-002` | The pipeline shall support paired FASTQ input preprocessing. | FASTQ rows are routed to read trimming, alignment, sorting, duplicate marking, filtering, and fragment extraction. | `main.nf`; `modules/local/fastp.nf`; `modules/local/align_bwamem2.nf`; `modules/local/align_bowtie2.nf`; `modules/local/sort_index_bam.nf` | `OUT-003`: trimmed FASTQ files; `OUT-004`: aligned/sorted BAM files | `VER-002`: Nextflow execution with `tests/tiny/samplesheet.csv` recommended | Implemented; integration verification recommended |
| `SR-003` | The pipeline shall support BAM input preprocessing. | BAM rows bypass FASTQ alignment and enter sorting/indexing, duplicate marking, filtering, and fragment extraction. | `main.nf`; `modules/local/sort_index_bam.nf`; `modules/local/mark_duplicates.nf`; `modules/local/filter_bam.nf` | `OUT-005`: filtered BAM files; `OUT-006`: BAM indexes | `VER-003`: Integration test with representative BAM recommended | Implemented; verification recommended |
| `SR-004` | The pipeline shall support CRAM input preprocessing when reference FASTA is supplied. | CRAM rows are routed into sorting/indexing and downstream BAM-derived fragment extraction. | `main.nf`; `nextflow.config`; `modules/local/sort_index_bam.nf` | `OUT-005`: filtered BAM files; `OUT-006`: BAM indexes | `VER-004`: Integration test with representative CRAM plus FASTA recommended | Implemented; verification recommended |
| `SR-005` | The pipeline shall convert filtered alignments to fragment records. | Filtered BAM files are converted into fragment tables and sample-level fragment QC. | `modules/local/bam_to_fragments.nf`; `bin/bam_to_fragments.py`; `modules/local/fragment_qc.nf`; `bin/fragment_qc.py` | `OUT-007`: `results/<run>/fragments/*.fragments.parquet`; `OUT-008`: `results/<run>/qc/*.fragment_qc.tsv` | `VER-005`: Unit/integration test for BAM-to-fragment conversion recommended | Implemented; verification recommended |
| `SR-006` | The pipeline shall support public five-column fragment files. | `input_type=finaledb_fragments` routes validated fragment files directly to ingest. | `main.nf`; `modules/local/ingest_finaledb_fragments.nf`; `bin/ingest_finaledb_fragments.py` | `OUT-009`: normalized fragment Parquet files; `OUT-010`: fragment ingest QC files | `VER-006`: `tests/test_finaledb_ingest.py::test_ingest_filters_mapq_and_computes_length_from_end_start` | Implemented |
| `SR-007` | The fragment ingest shall interpret five-column fragment files as `chrom`, `start`, `end`, `mapq`, `strand`. | Ingest parser treats column 4 as MAPQ and computes fragment length as `end - start`. | `bin/ingest_finaledb_fragments.py`; `docs/usage.md`; `README.md` | `OUT-009`: normalized fragment Parquet files | `VER-006`: `tests/test_finaledb_ingest.py::test_ingest_filters_mapq_and_computes_length_from_end_start` | Implemented |
| `SR-008` | The fragment ingest shall normalize chromosome naming to match generated bins. | Chromosome style is inferred from bins and input names are converted as needed. | `bin/ingest_finaledb_fragments.py`; `nextflow.config` parameter `chrom_style` | `OUT-009`: normalized fragment Parquet files | `VER-007`: `tests/test_finaledb_ingest.py::test_chrom_normalization` | Implemented |
| `SR-009` | The pipeline shall generate genomic bins from chromosome sizes. | `BUILD_BINS` creates autosome-filtered bins at configured resolutions. | `modules/local/build_bins.nf`; `bin/make_bins.py`; `nextflow.config` parameters `bin_small`, `bin_large`, `autosomes_only` | `OUT-011`: `resources/generated/*.autosomes.filtered.bed` | `VER-008`: Integration test comparing expected bin intervals recommended | Implemented; verification recommended |
| `SR-010` | The pipeline shall annotate bins with GC content when reference sequence is available. | `GC_ANNOTATE_BINS` computes GC annotations for generated bins. | `modules/local/gc_annotate_bins.nf`; `bin/compute_gc.py` | `OUT-012`: `resources/generated/*.bin_gc.tsv` | `VER-009`: Unit test with tiny FASTA and expected GC recommended | Implemented; verification recommended |
| `SR-011` | The pipeline shall count fragments overlapping genomic bins and compute length-class summaries. | `BIN_FRAGMENT_COUNTS` counts total, short, and long fragments by bin. | `modules/local/bin_fragment_counts.nf`; `bin/bin_fragment_counts.py`; `nextflow.config` parameters `short_min`, `short_max`, `long_min`, `long_max` | `OUT-013`: `results/<run>/bin_counts/*.counts.tsv.gz` | `VER-010`: `tests/test_bin_fragment_counts.py::test_vectorized_counter_counts_boundary_overlaps` | Implemented |
| `SR-012` | The pipeline shall perform GC correction on binned fragment counts. | `GC_CORRECT_BINS` applies GC correction to bin-level counts. | `modules/local/gc_correct_bins.nf`; `bin/gc_correct.py`; `nextflow.config` parameter `loess_frac` | `OUT-014`: `results/<run>/gc_corrected/*.gc_corrected.tsv.gz` | `VER-011`: Unit test using known GC/count fixture recommended | Implemented; verification recommended |
| `SR-013` | The pipeline shall assemble sample-level features into a model-ready matrix. | `ASSEMBLE_FEATURE_MATRIX` combines validated labels, GC-corrected bin features, and arm-level features. | `modules/local/assemble_feature_matrix.nf`; `bin/assemble_matrix.py` | `OUT-015`: `results/<run>/features/feature_matrix.tsv`; `OUT-016`: `feature_dictionary.tsv`; `OUT-017`: `feature_summary.json` | `VER-012`: Fixture-based matrix shape/content test recommended | Implemented; verification recommended |
| `SR-014` | The pipeline shall train a baseline model from a labeled feature matrix. | `TRAIN_MODEL` trains a scikit-learn model and stores metrics, metadata, selected features, and predictions. | `modules/local/train_model.nf`; `bin/train_model.py`; `nextflow.config` parameter `train_max_features` | `OUT-018`: `results/<run>/models/model.pkl`; `OUT-019`: `model_metadata.json`; `OUT-020`: `model_metrics.json`; `OUT-021`: `model_cv_predictions.tsv`; `OUT-022`: `model_test_predictions.tsv` | `VER-013`: Reproducibility test with fixed random state and tiny matrix recommended | Implemented; verification recommended |
| `SR-015` | The pipeline shall apply a saved model to a feature matrix and emit prediction scores. | `PREDICT_MODEL` loads a model artifact and writes prediction scores. | `modules/local/predict_model.nf`; `bin/predict.py` | `OUT-023`: `results/<run>/predictions/predictions.tsv` | `VER-014`: Fixture model prediction test recommended | Implemented; verification recommended |
| `SR-016` | The pipeline shall summarize prediction performance and score distributions. | `POST_ANALYZE_PREDICTIONS` combines predictions, labels, model metrics, metadata, CV predictions, and test predictions. | `modules/local/post_analyze_predictions.nf`; `bin/post_analyze_predictions.py` | `OUT-024`: `prediction_analysis_summary.json`; `OUT-025`: ROC/PR plots; `OUT-026`: confusion matrix and ranked prediction tables | `VER-015`: Report-generation test confirming required sections and metrics recommended | Implemented; verification recommended |
| `SR-017` | The pipeline shall produce a human-readable HTML report. | Post-analysis renders an HTML report with held-out metrics, plots, cohort details, and interpretation notes. | `bin/post_analyze_predictions.py`; `modules/local/post_analyze_predictions.nf` | `OUT-027`: `prediction_analysis_report.html` | `VER-016`: Representative report committed at `results/finaledb_lung_demo_300/analysis/prediction_analysis_report.html` | Implemented |
| `SR-018` | The pipeline shall use research-oriented prediction language by default. | Prediction labels are expressed as fragmentome score categories, and README caveats identify scores as research outputs. | `nextflow.config` parameter `research_language`; `bin/predict.py`; `bin/post_analyze_predictions.py`; `README.md`; `docs/usage.md` | `OUT-023`: prediction scores; `OUT-027`: report text | `VER-017`: Text assertion test for report language recommended | Implemented; verification recommended |
| `SR-019` | The repository shall exclude raw data and heavyweight runtime artifacts from source control while allowing small representative reports. | `.gitignore` excludes `work/`, bulk `results/`, raw fragments, reference genomes, model binaries, compressed genomics artifacts, and logs; selected small HTML/JSON reports are allowed. | `.gitignore`; root `README.md`; `docs/assets/example_prediction_analysis_report.png` | `OUT-028`: source-control-safe repository contents | `VER-018`: `git ls-files` size check and ignored-file review | Implemented |

## Design Requirement Traceability

| Design ID | Design Requirement | Satisfies System Requirement(s) | Implementation Reference | Verification |
| --- | --- | --- | --- | --- |
| `DR-001` | Workflow modes shall be explicit and limited to supported values. | `SR-001` through `SR-017` | `main.nf` mode validation; `docs/usage.md` | Run with invalid `--mode` should fail with a clear error. |
| `DR-002` | Required parameters shall be checked before dependent processes run. | `SR-001`, `SR-009`, `SR-014`, `SR-015`, `SR-016` | `main.nf` helper `requireParam` | Missing-parameter tests recommended. |
| `DR-003` | Public fragment ingest shall preserve sample traceability from manifest to normalized sample sheet. | `SR-001`, `SR-006`, `SR-007` | `bin/validate_samplesheet.py`; `bin/finaledb_metadata.py`; `bin/ingest_finaledb_fragments.py` | `VER-001`, `VER-006`. |
| `DR-004` | Feature extraction shall use generated bins and consistent chromosome naming. | `SR-008`, `SR-009`, `SR-011` | `bin/make_bins.py`; `bin/ingest_finaledb_fragments.py`; `bin/bin_fragment_counts.py` | `VER-007`, `VER-010`. |
| `DR-005` | Model outputs shall be stored with enough metadata to reconstruct evaluation context. | `SR-014`, `SR-016`, `SR-017` | `bin/train_model.py`; `bin/post_analyze_predictions.py` | `VER-013`, `VER-015`. |
| `DR-006` | Optional feature families shall be controlled by explicit parameters. | `SR-013` | `main.nf`; `nextflow.config` parameters `run_end_motifs`, `run_wps` | Parameterized workflow test recommended. |
| `DR-007` | The repository shall include a representative report but not bulk result directories. | `SR-019` | `.gitignore`; root `README.md` | `VER-018`. |

## Output Traceability

| Output ID | Artifact | Producing Process | Primary Consumer | Related Requirement(s) |
| --- | --- | --- | --- | --- |
| `OUT-001` | `validated_samplesheet.csv` | `VALIDATE_SAMPLESHEET` | All downstream workflows | `SR-001` |
| `OUT-002` | `run_manifest.json` | `VALIDATE_SAMPLESHEET` | Audit/review | `SR-001` |
| `OUT-003` | `*.trimmed.fastq.gz` | `FASTP` | Alignment | `SR-002` |
| `OUT-004` | `*.aligned.bam`, `*.sorted.bam` | `ALIGN_*`, `SORT_INDEX_BAM` | Duplicate marking/filtering | `SR-002`, `SR-003`, `SR-004` |
| `OUT-005` | `*.filtered.bam` | `FILTER_BAM` | Fragment extraction and QC | `SR-003`, `SR-004`, `SR-005` |
| `OUT-007` | `*.fragments.parquet` | `BAM_TO_FRAGMENTS` or `INGEST_FINALEDB_FRAGMENTS` | Bin counting and optional motif/WPS features | `SR-005`, `SR-006` |
| `OUT-010` | `*.fragment_ingest_qc.tsv` | `INGEST_FINALEDB_FRAGMENTS` | QC review | `SR-006`, `SR-007` |
| `OUT-011` | `*.autosomes.filtered.bed` | `BUILD_BINS` | GC annotation and fragment counting | `SR-009`, `SR-011` |
| `OUT-012` | `*.bin_gc.tsv` | `GC_ANNOTATE_BINS` | Fragment counting and GC correction | `SR-010`, `SR-012` |
| `OUT-013` | `*.counts.tsv.gz` | `BIN_FRAGMENT_COUNTS` | GC correction | `SR-011` |
| `OUT-014` | `*.gc_corrected.tsv.gz` | `GC_CORRECT_BINS` | Feature matrix assembly | `SR-012`, `SR-013` |
| `OUT-015` | `feature_matrix.tsv` | `ASSEMBLE_FEATURE_MATRIX` | Training and prediction | `SR-013`, `SR-014`, `SR-015` |
| `OUT-018` | `model.pkl` | `TRAIN_MODEL` | Prediction | `SR-014`, `SR-015` |
| `OUT-019` | `model_metadata.json` | `TRAIN_MODEL` | Post-analysis/reporting | `SR-014`, `SR-016` |
| `OUT-020` | `model_metrics.json` | `TRAIN_MODEL` | Post-analysis/reporting | `SR-014`, `SR-016` |
| `OUT-023` | `predictions.tsv` | `PREDICT_MODEL` | Post-analysis/reporting | `SR-015`, `SR-016` |
| `OUT-024` | `prediction_analysis_summary.json` | `POST_ANALYZE_PREDICTIONS` | Review and presentation | `SR-016`, `SR-017` |
| `OUT-027` | `prediction_analysis_report.html` | `POST_ANALYZE_PREDICTIONS` | Human review and presentation | `SR-017`, `SR-018` |

## Verification Matrix

| Verification ID | Verification Type | Current Evidence | Covers | Gap / Recommended Expansion |
| --- | --- | --- | --- | --- |
| `VER-001` | Unit test | `tests/test_finaledb_ingest.py::test_validate_samplesheet_accepts_finaledb_manifest_columns` | Manifest normalization and validated sample sheet generation | Add negative tests for missing required columns and invalid input paths. |
| `VER-006` | Unit/integration test | `tests/test_finaledb_ingest.py::test_ingest_filters_mapq_and_computes_length_from_end_start` | Fragment ingest schema, MAPQ filtering, length computation | Add fixture for length upper/lower bound filtering. |
| `VER-007` | Unit test | `tests/test_finaledb_ingest.py::test_chrom_normalization` | Chromosome style normalization | Add bin-style auto-detection regression test. |
| `VER-010` | Unit test | `tests/test_bin_fragment_counts.py::test_vectorized_counter_counts_boundary_overlaps` | Bin overlap counting and fragment length summaries | Add multi-sample and empty-bin cases. |
| `VER-016` | Representative artifact | `results/finaledb_lung_demo_300/analysis/prediction_analysis_report.html` | HTML report format and presentation | Add automated smoke test that validates expected report sections. |
| `VER-018` | Source-control safety check | `.gitignore`; no tracked files over 50 MB during repository review | Raw data and heavyweight artifact exclusion | Add pre-commit or CI check for large files and forbidden paths. |

## Risk Traceability

| Risk ID | Risk | Mitigation | Related Requirement(s) | Residual Risk |
| --- | --- | --- | --- | --- |
| `RISK-001` | Public fragment files may be misparsed if column 4 is treated as fragment length instead of MAPQ. | Parser explicitly treats column 4 as MAPQ and computes length from coordinates; unit test covers behavior. | `SR-006`, `SR-007` | Low for covered input format; new input formats need review. |
| `RISK-002` | Chromosome naming mismatch may silently drop fragments or bins. | Chromosome normalization and bin-style matching are implemented and tested. | `SR-008`, `SR-011` | Medium; add counts of dropped chromosomes to QC review. |
| `RISK-003` | Model metrics may be overinterpreted as clinical performance. | Research language, caveats, and score category wording are used in documentation and reports. | `SR-018` | Medium; reports should continue to identify validation design and cohort limitations. |
| `RISK-004` | Same-cohort training and prediction can inflate performance estimates. | Reports distinguish held-out, cross-validated, and full-cohort views; documentation states caveats. | `SR-014`, `SR-016`, `SR-017` | Medium; independent validation set support is recommended for stronger claims. |
| `RISK-005` | Large raw data or model artifacts may be accidentally committed. | `.gitignore` excludes raw fragments, `work/`, bulk `results/`, references, compressed genomics files, and model binaries. | `SR-019` | Low; add automated large-file checks before push. |

## Change Impact Checklist

Use this checklist when modifying the pipeline:

| Change Area | Required Traceability Update |
| --- | --- |
| New input type | Add or update `UN`, `SR`, parser design, outputs, and verification rows. |
| New feature family | Add design and output rows; map to `ASSEMBLE_FEATURE_MATRIX`; add fixture tests. |
| New model type | Update `SR-014`, model metadata expectations, and prediction verification. |
| Report layout or language changes | Update `SR-016`, `SR-017`, `SR-018`, and representative report evidence. |
| Source-control policy changes | Update `SR-019`, `.gitignore`, and large-file verification evidence. |

