process WPS_FEATURES {
  tag "$sample_id"
  publishDir "${params.outdir}/wps", mode: 'copy', pattern: '*.wps.tsv.gz'

  input:
  tuple val(sample_id), path(fragments)

  output:
  tuple val(sample_id), path("${sample_id}.wps.tsv.gz"), emit: features

  script:
  """
  wps.py --sample-id ${sample_id} --fragments ${fragments} --regions ${params.regulatory_bed ?: ''} --out ${sample_id}.wps.tsv.gz
  """
}

