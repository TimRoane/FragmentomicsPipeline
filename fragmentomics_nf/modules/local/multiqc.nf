process MULTIQC {
  tag 'multiqc'
  publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'multiqc_report.html'

  input:
  path inputs

  output:
  path 'multiqc_report.html', emit: report

  script:
  """
  multiqc . --filename multiqc_report.html
  """
}

