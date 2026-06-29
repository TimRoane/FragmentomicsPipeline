process SORT_INDEX_BAM {
  tag "$sample_id"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.sorted.bam*'

  input:
  tuple val(sample_id), path(input_bam)

  output:
  tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam

  script:
  """
  samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${input_bam}
  samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
  """
}

