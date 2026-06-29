process VALIDATE_SAMPLESHEET {
  tag 'samplesheet'
  publishDir "${params.outdir}", mode: 'copy', pattern: 'validated_samplesheet.csv'
  publishDir "${params.outdir}", mode: 'copy', pattern: 'run_manifest.json'

  input:
  path samplesheet

  output:
  path 'validated_samplesheet.csv', emit: validated
  path 'run_manifest.json', emit: manifest

  script:
  """
  validate_samplesheet.py \
    --samplesheet ${samplesheet} \
    --out validated_samplesheet.csv \
    --manifest run_manifest.json \
    --mode ${params.mode} \
    --genome ${params.genome} \
    --assembly ${params.assembly ?: params.genome} \
    --input-type ${params.input_type} \
    --path-base ${launchDir} \
    ${params.fasta ? "--fasta ${params.fasta}" : ""} \
    ${params.tissue_of_origin ? "--tissue-of-origin" : ""}
  """
}
