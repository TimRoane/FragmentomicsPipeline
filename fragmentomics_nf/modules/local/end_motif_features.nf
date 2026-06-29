process END_MOTIF_FEATURES {
  tag "$sample_id"
  publishDir "${params.outdir}/features", mode: 'copy', pattern: '*.end_motifs.tsv'

  input:
  tuple val(sample_id), path(fragments)

  output:
  tuple val(sample_id), path("${sample_id}.end_motifs.tsv"), emit: features

  script:
  """
  end_motifs.py --sample-id ${sample_id} --fragments ${fragments} --fasta ${params.fasta} --out ${sample_id}.end_motifs.tsv
  """
}

