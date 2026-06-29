process TRAIN_MODEL {
  tag 'train'
  publishDir "${params.outdir}/models", mode: 'copy', pattern: '*'

  input:
  path matrix

  output:
  path 'model.pkl', emit: model
  path 'model_metadata.json', emit: metadata
  path 'model_metrics.json', emit: metrics
  path 'model_cv_predictions.tsv', emit: cv_predictions
  path 'model_test_predictions.tsv', emit: test_predictions

  script:
  """
  train_model.py \
    --matrix ${matrix} \
    --out-model model.pkl \
    --out-metadata model_metadata.json \
    --out-metrics model_metrics.json \
    --out-cv-predictions model_cv_predictions.tsv \
    --out-test-predictions model_test_predictions.tsv \
    --max-features ${params.train_max_features}
  """
}
