process MITOCHONDRIAL_FEATURES {
  tag "$sample_id"
  publishDir "${params.outdir}/features", mode: 'copy', pattern: '*.mito_features.tsv'

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.mito_features.tsv"), emit: features

  script:
  """
  mito_features.py --sample-id ${sample_id} --bam ${bam} --out ${sample_id}.mito_features.tsv
  """
}

