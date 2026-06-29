#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_SAMPLESHEET } from './modules/local/validate_samplesheet'
include { FASTP } from './modules/local/fastp'
include { ALIGN_BWAMEM2 } from './modules/local/align_bwamem2'
include { ALIGN_BOWTIE2 } from './modules/local/align_bowtie2'
include { SORT_INDEX_BAM } from './modules/local/sort_index_bam'
include { MARK_DUPLICATES } from './modules/local/mark_duplicates'
include { FILTER_BAM } from './modules/local/filter_bam'
include { BAM_TO_FRAGMENTS } from './modules/local/bam_to_fragments'
include { INGEST_FINALEDB_FRAGMENTS } from './modules/local/ingest_finaledb_fragments'
include { FRAGMENT_QC } from './modules/local/fragment_qc'
include { BUILD_BINS } from './modules/local/build_bins'
include { GC_ANNOTATE_BINS } from './modules/local/gc_annotate_bins'
include { BIN_FRAGMENT_COUNTS } from './modules/local/bin_fragment_counts'
include { GC_CORRECT_BINS } from './modules/local/gc_correct_bins'
include { CHROMOSOME_ARM_FEATURES } from './modules/local/chromosome_arm_features'
include { MITOCHONDRIAL_FEATURES } from './modules/local/mitochondrial_features'
include { END_MOTIF_FEATURES } from './modules/local/end_motif_features'
include { WPS_FEATURES } from './modules/local/wps_features'
include { ASSEMBLE_FEATURE_MATRIX } from './modules/local/assemble_feature_matrix'
include { TRAIN_MODEL } from './modules/local/train_model'
include { PREDICT_MODEL } from './modules/local/predict_model'
include { POST_ANALYZE_PREDICTIONS } from './modules/local/post_analyze_predictions'
include { REPORT_SAMPLE } from './modules/local/report_sample'
include { MULTIQC } from './modules/local/multiqc'

def requireParam(name, value) {
  if (value == null || value.toString().trim() == '') {
    error "Missing required parameter: --${name}"
  }
}

def chromSizesParam() {
  return params.chrom_sizes ?: params.chrom_size
}

workflow PREPROCESS {
  take:
  samplesheet

  main:
  requireParam('samplesheet', params.samplesheet)
  VALIDATE_SAMPLESHEET(samplesheet)

  samples = VALIDATE_SAMPLESHEET.out.validated
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id as String, row) }

  fastq_samples = samples.filter { sample_id, row -> row.fastq_1 && row.fastq_2 }
  bam_samples = samples.filter { sample_id, row -> row.bam }
  cram_samples = samples.filter { sample_id, row -> row.cram }

  FASTP(fastq_samples)

  if (params.aligner == 'bowtie2') {
    ALIGN_BOWTIE2(FASTP.out.reads)
    aligned = ALIGN_BOWTIE2.out.bam.mix(bam_samples.map { id, row -> tuple(id, file(row.bam)) })
  } else {
    ALIGN_BWAMEM2(FASTP.out.reads)
    aligned = ALIGN_BWAMEM2.out.bam.mix(bam_samples.map { id, row -> tuple(id, file(row.bam)) })
  }

  cram_bams = cram_samples.map { id, row -> tuple(id, file(row.cram)) }
  SORT_INDEX_BAM(aligned.mix(cram_bams))
  MARK_DUPLICATES(SORT_INDEX_BAM.out.bam)
  FILTER_BAM(MARK_DUPLICATES.out.bam)
  BAM_TO_FRAGMENTS(FILTER_BAM.out.bam)
  FRAGMENT_QC(BAM_TO_FRAGMENTS.out.fragments, FILTER_BAM.out.bam)

  emit:
  validated = VALIDATE_SAMPLESHEET.out.validated
  manifest = VALIDATE_SAMPLESHEET.out.manifest
  fragments = BAM_TO_FRAGMENTS.out.fragments
  filtered_bams = FILTER_BAM.out.bam
  qc = FRAGMENT_QC.out.qc
}

workflow FEATURE_EXTRACT {
  take:
  fragments
  bams
  validated

  main:
  requireParam('chrom_sizes', chromSizesParam())
  BUILD_BINS(file(chromSizesParam()))
  GC_ANNOTATE_BINS(BUILD_BINS.out.bins)
  BIN_FRAGMENT_COUNTS(fragments, GC_ANNOTATE_BINS.out.annotated)
  GC_CORRECT_BINS(BIN_FRAGMENT_COUNTS.out.counts)
  CHROMOSOME_ARM_FEATURES(GC_CORRECT_BINS.out.corrected)
  MITOCHONDRIAL_FEATURES(bams)

  if (params.run_end_motifs) {
    END_MOTIF_FEATURES(fragments)
  }

  if (params.run_wps) {
    WPS_FEATURES(fragments)
  }

  ASSEMBLE_FEATURE_MATRIX(
    validated,
    GC_CORRECT_BINS.out.corrected.map { sample_id, corrected -> corrected }.collect(),
    CHROMOSOME_ARM_FEATURES.out.features.map { sample_id, features -> features }.collect()
  )

  emit:
  matrix = ASSEMBLE_FEATURE_MATRIX.out.matrix
  features = ASSEMBLE_FEATURE_MATRIX.out.sample_features
}

workflow {
  mode = params.mode
  if (!['preprocess', 'feature_extract', 'train', 'predict', 'post_analyze', 'full'].contains(mode)) {
    error "Unsupported --mode '${mode}'. Use preprocess, feature_extract, train, predict, post_analyze, or full."
  }

  if (mode in ['preprocess', 'full']) {
    PREPROCESS(file(params.samplesheet ?: params.input))
  }

  if (mode == 'preprocess') {
    MULTIQC(Channel.fromPath("${params.outdir}/**/*", checkIfExists: false))
  } else if (mode == 'feature_extract') {
    requireParam('samplesheet', params.samplesheet ?: params.input)
    VALIDATE_SAMPLESHEET(file(params.samplesheet ?: params.input))
    if (params.input_type == 'finaledb_fragments') {
      finaledb_ch = VALIDATE_SAMPLESHEET.out.validated
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id as String, file(row.finaledb_fragment ?: row.local_fragment ?: row.fragment ?: row.fragments ?: row.fragment_file)) }
      bins_for_style = Channel.value(file(chromSizesParam()))
      INGEST_FINALEDB_FRAGMENTS(finaledb_ch, bins_for_style)
      fragment_ch = INGEST_FINALEDB_FRAGMENTS.out.fragments
      bam_ch = Channel.empty()
    } else {
      fragment_ch = Channel.fromPath("${params.outdir}/fragments/*.fragments.parquet").map { f -> tuple(f.baseName.replace('.fragments', ''), f) }
      bam_ch = Channel.fromPath("${params.outdir}/bam/*.filtered.bam").map { f -> tuple(f.baseName.replace('.filtered', ''), f) }
    }
    FEATURE_EXTRACT(fragment_ch, bam_ch, VALIDATE_SAMPLESHEET.out.validated)
  } else if (mode == 'full') {
    FEATURE_EXTRACT(PREPROCESS.out.fragments, PREPROCESS.out.filtered_bams, PREPROCESS.out.validated)
    if (params.model) {
      PREDICT_MODEL(FEATURE_EXTRACT.out.matrix, file(params.model))
    }
    REPORT_SAMPLE(FEATURE_EXTRACT.out.features)
    MULTIQC(Channel.fromPath("${params.outdir}/**/*", checkIfExists: false))
  } else if (mode == 'train') {
    requireParam('feature_matrix', params.feature_matrix)
    TRAIN_MODEL(file(params.feature_matrix))
  } else if (mode == 'predict') {
    requireParam('feature_matrix', params.feature_matrix)
    requireParam('model', params.model)
    PREDICT_MODEL(file(params.feature_matrix), file(params.model))
  } else if (mode == 'post_analyze') {
    requireParam('feature_matrix', params.feature_matrix)
    requireParam('predictions', params.predictions)
    requireParam('model_metrics', params.model_metrics)
    requireParam('model_metadata', params.model_metadata)
    cvPredictions = params.cv_predictions ?: "${params.outdir}/models/model_cv_predictions.tsv"
    testPredictions = params.test_predictions ?: "${params.outdir}/models/model_test_predictions.tsv"
    POST_ANALYZE_PREDICTIONS(
      file(params.predictions),
      file(params.feature_matrix),
      file(params.model_metrics),
      file(params.model_metadata),
      file(cvPredictions),
      file(testPredictions)
    )
  }
}
