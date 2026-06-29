process FRAGMENT_QC {
  tag "$sample_id"
  publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*'

  input:
  tuple val(sample_id), path(fragments)
  tuple val(sample_id2), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.fragment_qc.tsv"), emit: qc

  script:
  """
  fragment_qc.py \
    --sample-id ${sample_id} \
    --fragments ${fragments} \
    --bam ${bam} \
    --out ${sample_id}.fragment_qc.tsv \
    --short-min ${params.short_min} --short-max ${params.short_max} \
    --long-min ${params.long_min} --long-max ${params.long_max}
  """
}

