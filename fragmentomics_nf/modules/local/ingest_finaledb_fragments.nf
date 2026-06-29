process INGEST_FINALEDB_FRAGMENTS {
  tag "$sample_id"
  publishDir "${params.outdir}/fragments", mode: 'copy', pattern: '*.fragments.*'
  publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.fragment_ingest_qc.tsv'

  input:
  tuple val(sample_id), path(fragment_file)
  path bins_for_style

  output:
  tuple val(sample_id), path("${sample_id}.fragments.parquet"), emit: fragments
  path "${sample_id}.fragments.bed.gz", optional: true, emit: bed
  tuple val(sample_id), path("${sample_id}.fragment_ingest_qc.tsv"), emit: qc

  script:
  def minLen = params.min_fragment_length ?: params.fragment_min_len
  def maxLen = params.max_fragment_length ?: params.fragment_max_len
  def bedOpt = params.write_fragment_bed ? "--out-bed ${sample_id}.fragments.bed.gz" : ""
  """
  ingest_finaledb_fragments.py \
    --sample-id ${sample_id} \
    --fragments ${fragment_file} \
    --out-parquet ${sample_id}.fragments.parquet \
    ${bedOpt} \
    --qc-out ${sample_id}.fragment_ingest_qc.tsv \
    --min-mapq ${params.min_mapq} \
    --min-fragment-length ${minLen} \
    --max-fragment-length ${maxLen} \
    --short-min ${params.short_min} \
    --short-max ${params.short_max} \
    --long-min ${params.long_min} \
    --long-max ${params.long_max} \
    --chrom-style ${params.chrom_style} \
    ${bins_for_style ? "--bins ${bins_for_style}" : ""}
  """
}
