process POST_ANALYZE_PREDICTIONS {
  tag 'post_analyze'
  publishDir "${params.outdir}/predictions", mode: 'copy', pattern: '*'

  input:
  path predictions
  path matrix
  path model_metrics
  path model_metadata
  path cv_predictions
  path test_predictions

  output:
  path 'predictions_with_labels.tsv', emit: joined
  path 'ranked_predictions.tsv', emit: ranked
  path 'prediction_analysis_summary.json', emit: summary
  path 'score_summary_by_label.tsv', emit: score_summary
  path 'confusion_matrix.tsv', optional: true, emit: confusion
  path 'score_histogram.png', emit: score_histogram
  path 'scores_by_label.png', emit: scores_by_label
  path 'ranked_scores.png', emit: ranked_scores
  path 'roc_curve.png', emit: roc_curve
  path 'precision_recall_curve.png', emit: precision_recall_curve
  path 'cv_roc_curve.png', emit: cv_roc_curve
  path 'cv_precision_recall_curve.png', emit: cv_precision_recall_curve
  path 'test_roc_curve.png', emit: test_roc_curve
  path 'test_precision_recall_curve.png', emit: test_precision_recall_curve
  path 'test_score_histogram.png', emit: test_score_histogram
  path 'prediction_analysis_report.html', emit: html

  script:
  """
  export MPLCONFIGDIR=\$PWD/.matplotlib
  post_analyze_predictions.py \
    --predictions ${predictions} \
    --matrix ${matrix} \
    --model-metrics ${model_metrics} \
    --model-metadata ${model_metadata} \
    --cv-predictions ${cv_predictions} \
    --test-predictions ${test_predictions} \
    --outdir .
  """
}
