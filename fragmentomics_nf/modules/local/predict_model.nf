process PREDICT_MODEL {
  tag 'predict'
  publishDir "${params.outdir}/predictions", mode: 'copy', pattern: '*'

  input:
  path matrix
  path model

  output:
  path 'predictions.tsv', emit: predictions

  script:
  """
  predict.py --matrix ${matrix} --model ${model} --out predictions.tsv --research-language ${params.research_language}
  """
}

