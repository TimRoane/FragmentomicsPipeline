process FASTP {
  tag "$sample_id"
  publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.fastp.*'
  publishDir "${params.outdir}/fastq", mode: 'copy', pattern: '*.trimmed.fastq.gz'

  input:
  tuple val(sample_id), val(row)

  output:
  tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: reads
  path "${sample_id}.fastp.html", emit: html
  path "${sample_id}.fastp.json", emit: json

  script:
  def detect = params.fastp_detect_adapter_for_pe ? '--detect_adapter_for_pe' : ''
  """
  fastp \
    --in1 ${row.fastq_1} \
    --in2 ${row.fastq_2} \
    --out1 ${sample_id}_R1.trimmed.fastq.gz \
    --out2 ${sample_id}_R2.trimmed.fastq.gz \
    --html ${sample_id}.fastp.html \
    --json ${sample_id}.fastp.json \
    --length_required ${params.fastp_min_length} \
    --qualified_quality_phred ${params.fastp_qualified_quality_phred} \
    ${detect}
  """
}

