process CHROMOSOME_ARM_FEATURES {
  tag "$sample_id"
  publishDir "${params.outdir}/features", mode: 'copy', pattern: '*.arm_features.tsv'

  input:
  tuple val(sample_id), path(corrected)

  output:
  tuple val(sample_id), path("${sample_id}.arm_features.tsv"), emit: features

  script:
  def arms = params.chromosome_arms_bed ? "--arms ${params.chromosome_arms_bed}" : ""
  def ref = params.healthy_reference_panel ? "--healthy-reference ${params.healthy_reference_panel}" : ""
  """
  arm_features.py \
    --sample-id ${sample_id} \
    --corrected ${corrected} \
    ${arms} \
    ${ref} \
    --out ${sample_id}.arm_features.tsv
  """
}
