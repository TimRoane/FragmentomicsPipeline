process BAM_TO_FRAGMENTS {
  tag "$sample_id"
  publishDir "${params.outdir}/fragments", mode: 'copy', pattern: '*'

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.fragments.parquet"), emit: fragments
  path "${sample_id}.fragments.bed.gz", emit: bed

  script:
  """
  bam_to_fragments.py \
    --bam ${bam} \
    --sample-id ${sample_id} \
    --out-parquet ${sample_id}.fragments.parquet \
    --out-bed ${sample_id}.fragments.bed.gz \
    --min-len ${params.fragment_min_len} \
    --max-len ${params.fragment_max_len} \
    --min-mapq ${params.min_mapq} \
    ${params.autosomes_only ? "--autosomes-only" : ""}
  """
}

