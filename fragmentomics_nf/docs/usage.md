# Usage

This pipeline has five modes:

- `preprocess`: FASTQ, BAM, or CRAM to filtered BAM and fragment files.
- `feature_extract`: existing fragment/BAM outputs to feature vectors.
- `train`: feature matrix plus labels to a transparent baseline model.
- `predict`: feature matrix plus saved model to research prediction output.
- `post_analyze`: prediction scores plus feature matrix/model metadata to reports, plots, and summary metrics.
- `full`: preprocess, feature extraction, QC, reports, and optional prediction.

Example:

```bash
nextflow run main.nf \
  -profile docker \
  --mode full \
  --samplesheet tests/tiny/samplesheet.csv \
  --genome hg38 \
  --fasta /path/to/hg38.fa \
  --chrom_sizes /path/to/hg38.chrom.sizes \
  --chromosome_arms_bed /path/to/chromosome_arms.bed \
  --outdir results
```

The sample sheet columns are:

```text
sample_id,fastq_1,fastq_2,bam,cram,label,cancer_type,batch,sex,age,smoking_status
```

Provide exactly one input type per row: paired FASTQ, BAM, or CRAM. CRAM requires `--fasta`.

For FinaleDB fragment ingest, use `--input_type finaledb_fragments`. The sample sheet may provide `finaledb_fragment`, `fragment_file`, `fragment`, `fragments`, or `metadata_json`. Metadata JSON is parsed structurally from `analysis.{assembly}[]` and only `desc == "fragment"` / `type == "tsv"` entries are selected.

The `finaledb_lung_demo_50/manifest.tsv` format is also accepted directly. It is normalized from `finaledb_id`, `cohort`, `local_fragment`, `fragment_url`, and `metadata_json` into the pipeline's canonical validated sample sheet.

```bash
nextflow run main.nf \
  --mode feature_extract \
  --input_type finaledb_fragments \
  --input finaledb_lung_demo_50/manifest.tsv \
  --assembly hg19 \
  --chrom_sizes /path/to/hg19.chrom.sizes \
  --outdir results/finaledb_lung_demo_50
```

Performance knobs for local FinaleDB runs:

```bash
nextflow run main.nf \
  --mode feature_extract \
  --input_type finaledb_fragments \
  --input finaledb_lung_demo_50/manifest.tsv \
  --assembly hg19 \
  --genome hg19 \
  --chrom_sizes /path/to/hg19.chrom.sizes \
  --outdir results/finaledb_lung_demo_50 \
  --bin_count_batch_size 1000000 \
  --write_fragment_bed false \
  -process.maxForks 8 \
  -resume
```

`--write_fragment_bed false` is the default for FinaleDB inputs; feature extraction uses Parquet fragments and does not need the extra BED file. Increase `-process.maxForks` cautiously on large runs because disk I/O can become the limiting factor before CPU.

Prediction output uses research language by default. It reports `cancer_probability` as a model probability and class labels as fragmentome score categories, not diagnostic calls.
