# fragmentomics_nf

Production-style Nextflow DSL2 pipeline for public-literature cfDNA WGS fragmentomics.

This project implements transparent research fragmentomics features: filtered cfDNA fragments, short/long coverage, GC-corrected bins, chromosome-arm summaries, feature matrices, baseline machine-learning models, and research prediction reports.

## What This Pipeline Does

Pipeline modes:

- `preprocess`: FASTQ/BAM/CRAM to filtered fragments.
- `feature_extract`: fragment files to binned fragmentomic feature matrix.
- `train`: feature matrix plus labels to a baseline research model.
- `predict`: feature matrix plus model to prediction scores.
- `full`: intended combined workflow for conventional sequencing inputs.

For FinaleDB fragment files, the usual workflow is:

```text
FinaleDB .frag.tsv.gz/.bgz
  -> validate manifest
  -> ingest/filter fragments
  -> count fragments in genomic bins
  -> GC correction
  -> assemble feature matrix
  -> train model
  -> predict / analyze scores
```

## Current Research Caveats

- Outputs are for research use only.
- `cancer_probability` is a model score, not a clinical diagnosis.
- If you train and predict on the same 50-sample demo cohort, predictions are a sanity check, not an independent validation.
- For real model evaluation, train on one cohort and predict on held-out samples.
- FinaleDB fragments use a five-column schema where column 4 is MAPQ, not length. This pipeline computes `length = end - start`.

## Repository Layout

```text
fragmentomics_nf/
  main.nf
  nextflow.config
  modules/local/
  bin/
  docs/
  assets/
  containers/
  tests/
  finaledb_lung_demo_50/
  results/
```

Important outputs:

```text
results/<run>/validated_samplesheet.csv
results/<run>/run_manifest.json
results/<run>/qc/*.fragment_ingest_qc.tsv
results/<run>/fragments/*.fragments.parquet
results/<run>/bin_counts/*.counts.tsv.gz
results/<run>/gc_corrected/*.gc_corrected.tsv.gz
results/<run>/features/feature_matrix.tsv
results/<run>/features/feature_dictionary.tsv
results/<run>/features/feature_summary.json
results/<run>/models/model.pkl
results/<run>/models/model_metadata.json
results/<run>/models/model_metrics.json
results/<run>/predictions/predictions.tsv
```

## Requirements

Minimum local tools:

- Nextflow
- Python with `pandas`, `numpy`, `pyarrow`, `scikit-learn`, `joblib`, `statsmodels`
- For FASTQ/BAM/CRAM paths: `samtools`, `bwa-mem2` or `bowtie2`, `fastp`
- Optional: Docker/Singularity if running through containers

Check Nextflow:

```bash
nextflow -version
```

Run commands from the pipeline directory:

```bash
cd /home/tim/RoanePortfolio/FragmentomicsPipeline/fragmentomics_nf
```

## Reference Files

For the FinaleDB lung demo, you need hg19 chromosome sizes. FASTA is recommended for GC annotation.

Directory:

```bash
mkdir -p resources/reference/hg19
```

Download hg19 chromosome sizes:

```bash
wget -O resources/reference/hg19/hg19.chrom.sizes \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
```

Download hg19 FASTA:

```bash
wget -O resources/reference/hg19/hg19.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

gunzip resources/reference/hg19/hg19.fa.gz
samtools faidx resources/reference/hg19/hg19.fa
```

If you already have a FASTA index:

```bash
cut -f1,2 hg19.fa.fai > resources/reference/hg19/hg19.chrom.sizes
```

Use genome-build-consistent resources. FinaleDB demo files in this repo are hg19, so use `--assembly hg19 --genome hg19`.

## FinaleDB Input Format

FinaleDB fragment files are five-column, headerless, BED-like files:

```text
chrom    start    end    mapq    strand
```

Example rows:

```text
1    10000    10168    0     +
1    10033    10215    35    -
2    20000    20130    60    +
```

The pipeline:

- treats column 4 as `mapq`
- computes `length = end - start`
- filters `mapq >= --min_mapq`
- filters `--min_fragment_length <= length <= --max_fragment_length`
- normalizes chromosome names to match bins

Accepted extensions:

```text
.frag.tsv.bgz
.frag.tsv.gz
.tsv.bgz
.tsv.gz
.bed.bgz
.bed.gz
.bgz
.gz
```

The included demo manifest is:

```text
finaledb_lung_demo_50/manifest.tsv
```

It contains `finaledb_id`, `cohort`, `local_fragment`, `fragment_url`, and `metadata_json`. The validator normalizes it into the pipeline sample sheet format.

## Step 1: Feature Extraction From FinaleDB Demo Data

Run:

```bash
nextflow run main.nf \
  --mode feature_extract \
  --input_type finaledb_fragments \
  --input finaledb_lung_demo_50/manifest.tsv \
  --assembly hg19 \
  --genome hg19 \
  --chrom_sizes resources/reference/hg19/hg19.chrom.sizes \
  --fasta resources/reference/hg19/hg19.fa \
  --outdir results/finaledb_lung_demo_50 \
  --write_fragment_bed false \
  --bin_count_batch_size 1000000 \
  -process.maxForks 8 \
  -resume
```

Notes:

- `--write_fragment_bed false` is the default and is faster. Feature extraction uses Parquet fragments and does not need BED files.
- `-process.maxForks 8` is a good starting point on a 32-core machine. Raise or lower it based on memory and disk I/O.
- `--chrom_size` is also accepted as an alias, but `--chrom_sizes` is canonical.

Expected key outputs:

```text
results/finaledb_lung_demo_50/qc/*.fragment_ingest_qc.tsv
results/finaledb_lung_demo_50/fragments/*.fragments.parquet
results/finaledb_lung_demo_50/bin_counts/*.counts.tsv.gz
results/finaledb_lung_demo_50/gc_corrected/*.gc_corrected.tsv.gz
results/finaledb_lung_demo_50/features/feature_matrix.tsv
results/finaledb_lung_demo_50/features/feature_dictionary.tsv
results/finaledb_lung_demo_50/features/feature_summary.json
```

Sanity checks:

```bash
wc -l results/finaledb_lung_demo_50/features/feature_matrix.tsv

cat results/finaledb_lung_demo_50/features/feature_summary.json

head -1 results/finaledb_lung_demo_50/features/feature_matrix.tsv | tr '\t' '\n' | wc -l
```

For the 50-sample demo, expect 51 lines including the header.

Inspect fragment ingest QC:

```bash
head -5 results/finaledb_lung_demo_50/qc/*.fragment_ingest_qc.tsv
```

Important QC columns:

```text
input_fragments
fragments_passing_mapq
fragments_passing_mapq_and_length
mean_length
median_length
fraction_short_100_150
fraction_long_151_220
short_long_ratio
```

## Step 2: Train A Research Model

Run:

```bash
nextflow run main.nf \
  --mode train \
  --feature_matrix results/finaledb_lung_demo_50/features/feature_matrix.tsv \
  --outdir results/finaledb_lung_demo_50 \
  --train_max_features 5000 \
  -resume
```

Training behavior:

- `control`, `healthy`, `normal`, `negative`, `0`, `false` are treated as negative.
- Other disease/cohort labels such as `lung_cancer` are treated as positive in binary mode.
- Constant and all-missing features are dropped.
- The trainer performs a stratified train/test split before model fitting.
- Feature selection is performed inside the model pipeline on training folds to reduce leakage.
- If the matrix is very wide, the trainer keeps up to `--train_max_features` features inside the model pipeline.
- Baseline models currently include logistic L2, random forest, and gradient boosting.
- The model with the best training-partition cross-validated ROC-AUC is saved.
- The saved model is fit on the training partition only, and held-out test predictions are written separately.

Expected outputs:

```text
results/finaledb_lung_demo_50/models/model.pkl
results/finaledb_lung_demo_50/models/model_metadata.json
results/finaledb_lung_demo_50/models/model_metrics.json
results/finaledb_lung_demo_50/models/model_cv_predictions.tsv
results/finaledb_lung_demo_50/models/model_test_predictions.tsv
```

Inspect model metadata:

```bash
cat results/finaledb_lung_demo_50/models/model_metadata.json
```

Inspect model metrics:

```bash
cat results/finaledb_lung_demo_50/models/model_metrics.json
```

Quick summary:

```bash
python - <<'PY'
import json
m = json.load(open("results/finaledb_lung_demo_50/models/model_metrics.json"))
for model, vals in m.items():
    print(model, "ROC-AUC=", vals.get("roc_auc"), "PR-AUC=", vals.get("pr_auc"))
PY
```

## Step 3: Predict Scores

Run:

```bash
nextflow run main.nf \
  --mode predict \
  --feature_matrix results/finaledb_lung_demo_50/features/feature_matrix.tsv \
  --model results/finaledb_lung_demo_50/models/model.pkl \
  --outdir results/finaledb_lung_demo_50 \
  -resume
```

Expected output:

```text
results/finaledb_lung_demo_50/predictions/predictions.tsv
```

Prediction columns:

```text
sample_id
cancer_probability
classification
threshold_used
model_version
qc_pass
qc_flags
top_contributing_features
```

Interpretation:

- `cancer_probability` is a model probability/fragmentome score.
- `classification` is research language: `elevated_fragmentome_score` or `lower_fragmentome_score`.
- `threshold_used` is taken from the trained model metadata.
- These are not clinical calls.

## Step 4: Post-Analyze Predictions

Run the built-in post-analysis mode:

```bash
nextflow run main.nf \
  --mode post_analyze \
  --feature_matrix results/finaledb_lung_demo_50/features/feature_matrix.tsv \
  --predictions results/finaledb_lung_demo_50/predictions/predictions.tsv \
  --model_metrics results/finaledb_lung_demo_50/models/model_metrics.json \
  --model_metadata results/finaledb_lung_demo_50/models/model_metadata.json \
  --outdir results/finaledb_lung_demo_50 \
  -resume
```

This generates:

```text
results/finaledb_lung_demo_50/predictions/predictions_with_labels.tsv
results/finaledb_lung_demo_50/predictions/ranked_predictions.tsv
results/finaledb_lung_demo_50/predictions/prediction_analysis_summary.json
results/finaledb_lung_demo_50/predictions/score_summary_by_label.tsv
results/finaledb_lung_demo_50/predictions/confusion_matrix.tsv
results/finaledb_lung_demo_50/predictions/score_histogram.png
results/finaledb_lung_demo_50/predictions/scores_by_label.png
results/finaledb_lung_demo_50/predictions/ranked_scores.png
results/finaledb_lung_demo_50/predictions/roc_curve.png
results/finaledb_lung_demo_50/predictions/precision_recall_curve.png
results/finaledb_lung_demo_50/predictions/cv_roc_curve.png
results/finaledb_lung_demo_50/predictions/cv_precision_recall_curve.png
results/finaledb_lung_demo_50/predictions/test_score_histogram.png
results/finaledb_lung_demo_50/predictions/test_roc_curve.png
results/finaledb_lung_demo_50/predictions/test_precision_recall_curve.png
results/finaledb_lung_demo_50/predictions/prediction_analysis_report.html
```

Open the HTML report:

```bash
xdg-open results/finaledb_lung_demo_50/predictions/prediction_analysis_report.html
```

Or inspect the summary JSON:

```bash
cat results/finaledb_lung_demo_50/predictions/prediction_analysis_summary.json
```

The post-analysis step:

- joins predictions back to sample labels and metadata
- writes ranked predictions
- computes score summaries by label
- computes confusion matrix, sensitivity, specificity, ROC-AUC, and PR-AUC when labels are available
- creates histogram, boxplot, ranked-score, ROC, and precision-recall plots
- includes separate cross-validated and held-out test ROC/PR curves when training outputs are present
- writes a single HTML report

If all samples receive the same classification, the summary JSON includes a `classification_warning`. This can happen when the threshold is poorly calibrated, when many scores tie at the threshold, or when predicting on the same cohort used for training.

## What To Look For After Prediction

For a labeled demo cohort:

- Do lung cancer samples rank above controls by `cancer_probability`?
- Are controls mostly below the model threshold?
- Are there outlier controls with high scores?
- Are there cancer samples with low scores?
- Does performance agree with `model_metrics.json`?
- Prefer held-out test and cross-validated prediction outputs over full-matrix predictions for evaluation.

For a new unlabeled cohort:

- Focus on score distribution and outliers.
- Do not compute ROC-AUC without labels.
- Treat `elevated_fragmentome_score` as a research flag for follow-up analysis, not a diagnosis.
- Check QC before trusting extreme scores.

Useful files to inspect together:

```text
results/<run>/qc/*.fragment_ingest_qc.tsv
results/<run>/features/feature_summary.json
results/<run>/models/model_metadata.json
results/<run>/models/model_metrics.json
results/<run>/predictions/predictions.tsv
```

## Conventional FASTQ/BAM/CRAM Inputs

Canonical sample sheet columns:

```text
sample_id,fastq_1,fastq_2,bam,cram,label,cancer_type,batch,sex,age,smoking_status
```

Rules:

- Provide exactly one input type per sample: paired FASTQ, BAM, or CRAM.
- FASTQ requires both `fastq_1` and `fastq_2`.
- CRAM requires `--fasta`.
- Labels are required for training.

Preprocess FASTQ/BAM/CRAM:

```bash
nextflow run main.nf \
  --mode preprocess \
  --samplesheet samplesheet.csv \
  --genome hg38 \
  --fasta resources/reference/hg38/hg38.fa \
  --chrom_sizes resources/reference/hg38/hg38.chrom.sizes \
  --aligner bwamem2 \
  --outdir results/my_run \
  -resume
```

Feature extraction from preprocessed outputs:

```bash
nextflow run main.nf \
  --mode feature_extract \
  --samplesheet samplesheet.csv \
  --genome hg38 \
  --fasta resources/reference/hg38/hg38.fa \
  --chrom_sizes resources/reference/hg38/hg38.chrom.sizes \
  --outdir results/my_run \
  -resume
```

Full mode for conventional inputs:

```bash
nextflow run main.nf \
  --mode full \
  --samplesheet samplesheet.csv \
  --genome hg38 \
  --fasta resources/reference/hg38/hg38.fa \
  --chrom_sizes resources/reference/hg38/hg38.chrom.sizes \
  --aligner bwamem2 \
  --outdir results/my_run \
  -process.maxForks 8 \
  -resume
```

## Important Parameters

Common parameters:

```text
--mode                         preprocess | feature_extract | train | predict | post_analyze | full
--outdir                       output directory
--genome                       hg19 or hg38
--assembly                     assembly for FinaleDB metadata, usually hg19
--chrom_sizes                  chromosome sizes file
--fasta                        reference FASTA
--input_type                   auto or finaledb_fragments
--input                        FinaleDB manifest/sample sheet
--samplesheet                  conventional sample sheet
--feature_matrix               feature matrix for train, predict, and post_analyze
--model                        model.pkl for prediction
--predictions                  predictions.tsv for post_analyze
--model_metrics                model_metrics.json for post_analyze
--model_metadata               model_metadata.json for post_analyze
```

Fragment filters:

```text
--min_mapq                     default 30
--min_fragment_length          default 50
--max_fragment_length          default 500
--short_min                    default 100
--short_max                    default 150
--long_min                     default 151
--long_max                     default 220
```

Performance:

```text
--write_fragment_bed           default false
--bin_count_batch_size         default 1000000
--bin_count_hist_max_length    default 1000
--train_max_features           default 5000
-process.maxForks              max concurrent Nextflow tasks
```

For a 32-core local machine, start with:

```text
-process.maxForks 8
```

Increase only if memory and disk I/O remain comfortable.

## Troubleshooting

Missing chromosome sizes:

```text
FileNotFoundError: hg19.chrom.sizes
```

Fix by passing the real path:

```bash
--chrom_sizes resources/reference/hg19/hg19.chrom.sizes
```

Relative FinaleDB metadata path not found:

```text
could not parse FinaleDB metadata_json
```

The validator resolves paths relative to the Nextflow launch directory. Run from the project directory and use `-resume`.

Training cannot find `feature_matrix.tsv`:

```text
FileNotFoundError: feature_matrix.tsv
```

Feature extraction did not finish or matrix assembly did not publish. Check:

```bash
ls -lh results/<run>/features/feature_matrix.tsv
cat results/<run>/features/feature_summary.json
```

Prediction says many samples are elevated:

- Check whether you predicted on the training set.
- Inspect `threshold_used`.
- Compare against `model_metrics.json`.
- Use an independent held-out cohort for meaningful model validation.

## More Documentation

See:

- [docs/usage.md](docs/usage.md)
- [docs/features.md](docs/features.md)
- [docs/data_sources.md](docs/data_sources.md)
- [docs/architecture_workflow.md](docs/architecture_workflow.md)
- [docs/traceability_matrix.md](docs/traceability_matrix.md)
