# Architecture Workflow

This diagram summarizes the `fragmentomics_nf` pipeline architecture from inputs through preprocessing, feature extraction, modeling, and reporting.

## Rendered Visual

- [Download the architecture workflow PDF](assets/architecture_workflow_visual.pdf)
- [Open the architecture workflow PNG](assets/architecture_workflow_visual.png)
- [Open the printable HTML source](assets/architecture_workflow_visual.html)

![Architecture workflow visual](assets/architecture_workflow_visual.png)

## Mermaid Source

```mermaid
flowchart TD
  subgraph INPUTS["Input Layer"]
    SS["Sample sheet or manifest"]
    FASTQ["Paired FASTQ"]
    BAM["BAM"]
    CRAM["CRAM"]
    FRAG["Public fragment files<br/>.frag.tsv.gz / .bgz"]
    REF["Reference resources<br/>FASTA, chrom sizes, optional masks"]
  end

  subgraph VALIDATION["Validation And Routing"]
    VS["VALIDATE_SAMPLESHEET<br/>validated_samplesheet.csv<br/>run_manifest.json"]
    ROUTE{"Input type / mode"}
  end

  subgraph PREPROCESS["Preprocess Workflow"]
    FASTP["FASTP<br/>read trimming and QC"]
    ALIGN["ALIGN_BWAMEM2 or ALIGN_BOWTIE2<br/>alignment"]
    SORT["SORT_INDEX_BAM"]
    DEDUP["MARK_DUPLICATES"]
    FILTER["FILTER_BAM<br/>MAPQ / flags / length-ready alignments"]
    BAMFRAG["BAM_TO_FRAGMENTS<br/>fragment Parquet"]
    FRAGQC["FRAGMENT_QC<br/>sample QC"]
  end

  subgraph PUBLIC_FRAGMENTS["Public Fragment Ingest"]
    INGEST["INGEST_FINALEDB_FRAGMENTS<br/>schema normalization<br/>MAPQ and length filtering<br/>chromosome style matching"]
    INGESTQC["fragment_ingest_qc.tsv"]
  end

  subgraph FEATURE_ENGINE["Feature Engineering"]
    BINS["BUILD_BINS<br/>autosomal genomic bins"]
    GCANN["GC_ANNOTATE_BINS<br/>bin GC content"]
    COUNTS["BIN_FRAGMENT_COUNTS<br/>total / short / long counts"]
    GCCORR["GC_CORRECT_BINS<br/>GC-adjusted bin counts"]
    ARM["CHROMOSOME_ARM_FEATURES"]
    MITO["MITOCHONDRIAL_FEATURES<br/>BAM-backed workflows"]
    MOTIF["END_MOTIF_FEATURES<br/>optional"]
    WPS["WPS_FEATURES<br/>optional"]
    MATRIX["ASSEMBLE_FEATURE_MATRIX<br/>feature_matrix.tsv<br/>feature_dictionary.tsv<br/>feature_summary.json"]
  end

  subgraph MODELING["Modeling And Analysis"]
    TRAIN["TRAIN_MODEL<br/>model.pkl<br/>model_metrics.json<br/>model_metadata.json"]
    PREDICT["PREDICT_MODEL<br/>predictions.tsv"]
    ANALYZE["POST_ANALYZE_PREDICTIONS<br/>summary JSON<br/>plots<br/>ranked tables"]
    REPORT["prediction_analysis_report.html"]
  end

  subgraph RUNTIME["Execution And Runtime Controls"]
    NF["Nextflow DSL2 orchestration<br/>main.nf"]
    CONFIG["nextflow.config<br/>params and process resources"]
    PROFILES["Profiles<br/>docker / singularity / slurm / aws / gcp"]
    OUTDIR["results/&lt;run&gt;/"]
    WORK["work/<br/>ignored runtime cache"]
  end

  SS --> VS
  VS --> ROUTE
  ROUTE -->|"FASTQ"| FASTP
  ROUTE -->|"BAM / CRAM"| SORT
  ROUTE -->|"finaledb_fragments"| INGEST

  FASTQ --> FASTP
  FASTP --> ALIGN
  ALIGN --> SORT
  BAM --> SORT
  CRAM --> SORT
  SORT --> DEDUP
  DEDUP --> FILTER
  FILTER --> BAMFRAG
  FILTER --> FRAGQC

  FRAG --> INGEST
  INGEST --> INGESTQC
  INGEST --> COUNTS
  BAMFRAG --> COUNTS

  REF --> BINS
  REF --> GCANN
  BINS --> GCANN
  GCANN --> COUNTS
  COUNTS --> GCCORR
  GCCORR --> ARM
  BAMFRAG --> MOTIF
  BAMFRAG --> WPS
  FILTER --> MITO
  ARM --> MATRIX
  GCCORR --> MATRIX
  MOTIF -. optional .-> MATRIX
  WPS -. optional .-> MATRIX
  MITO -. when BAM available .-> MATRIX

  MATRIX --> TRAIN
  MATRIX --> PREDICT
  TRAIN --> PREDICT
  TRAIN --> ANALYZE
  PREDICT --> ANALYZE
  ANALYZE --> REPORT

  NF -. orchestrates .-> VALIDATION
  NF -. orchestrates .-> PREPROCESS
  NF -. orchestrates .-> PUBLIC_FRAGMENTS
  NF -. orchestrates .-> FEATURE_ENGINE
  NF -. orchestrates .-> MODELING
  CONFIG -. controls .-> NF
  PROFILES -. configure runtime .-> NF
  MODELING --> OUTDIR
  FEATURE_ENGINE --> OUTDIR
  PREPROCESS --> OUTDIR
  PUBLIC_FRAGMENTS --> OUTDIR
  NF --> WORK
```

## Layer Summary

| Layer | Responsibility | Key Files |
| --- | --- | --- |
| Input layer | Accept user-provided sample metadata, sequencing files, fragment files, and reference resources. | `README.md`, `docs/usage.md`, `assets/schema_input.json`, `assets/schema_params.json` |
| Validation and routing | Normalize input metadata and choose the correct workflow path. | `main.nf`, `modules/local/validate_samplesheet.nf`, `bin/validate_samplesheet.py` |
| Preprocess workflow | Convert FASTQ/BAM/CRAM inputs into filtered alignments and fragment records. | `modules/local/fastp.nf`, `align_*.nf`, `sort_index_bam.nf`, `mark_duplicates.nf`, `filter_bam.nf`, `bam_to_fragments.nf` |
| Public fragment ingest | Normalize public fragment files into the same fragment representation used by downstream feature extraction. | `modules/local/ingest_finaledb_fragments.nf`, `bin/ingest_finaledb_fragments.py`, `bin/finaledb_metadata.py` |
| Feature engineering | Build bins, annotate GC, count fragments, correct counts, and assemble a feature matrix. | `modules/local/build_bins.nf`, `gc_annotate_bins.nf`, `bin_fragment_counts.nf`, `gc_correct_bins.nf`, `assemble_feature_matrix.nf` |
| Modeling and analysis | Train baseline models, generate predictions, summarize performance, and render reports. | `modules/local/train_model.nf`, `predict_model.nf`, `post_analyze_predictions.nf`, `bin/train_model.py`, `bin/predict.py`, `bin/post_analyze_predictions.py` |
| Runtime controls | Configure Nextflow execution, process resources, profiles, output locations, and ignored runtime cache. | `nextflow.config`, `conf/*.config`, `.gitignore` |

## Primary Data Flow

```text
sample metadata
  -> validation
  -> fragments
  -> bin-level counts
  -> GC-corrected features
  -> feature matrix
  -> model training / prediction
  -> prediction analysis report
```

## Architecture Notes

- `main.nf` is the orchestration boundary; most domain behavior lives in Python scripts under `bin/`.
- Nextflow modules under `modules/local/` are thin process wrappers that publish outputs into `results/<run>/`.
- Public fragment ingest and BAM-derived fragment extraction converge on Parquet fragment files before feature extraction.
- Optional feature families are parameter-controlled and can be added to the matrix without changing the core count and GC-correction path.
- Runtime cache and heavyweight generated outputs remain outside source control; only small representative report artifacts are tracked.
