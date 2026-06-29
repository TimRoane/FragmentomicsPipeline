process REPORT_SAMPLE {
  tag "$sample_id"
  publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.html'

  input:
  tuple val(sample_id), path(features)

  output:
  path "${sample_id}.html", emit: html

  script:
  """
  report_sample.py --sample-id ${sample_id} --features ${features} --out ${sample_id}.html
  """
}

