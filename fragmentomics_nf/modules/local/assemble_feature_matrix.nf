process ASSEMBLE_FEATURE_MATRIX {
  tag 'cohort'
  publishDir "${params.outdir}/features", mode: 'copy', pattern: '*'

  input:
  path validated
  path corrected
  path arms

  output:
  path 'feature_matrix.tsv', emit: matrix
  path 'feature_dictionary.tsv', emit: dictionary
  path 'feature_summary.json', emit: summary
  tuple val('cohort'), path('feature_matrix.tsv'), emit: sample_features

  script:
  """
  assemble_matrix.py \
    --samplesheet ${validated} \
    --corrected ${corrected.join(' ')} \
    --arms ${arms.join(' ')} \
    --out-matrix feature_matrix.tsv \
    --out-dictionary feature_dictionary.tsv \
    --out-summary feature_summary.json
  """
}
