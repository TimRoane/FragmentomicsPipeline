process MARK_DUPLICATES {
  tag "$sample_id"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.dedup.bam*'
  publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.duplication_metrics.txt'

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: bam
  path "${sample_id}.duplication_metrics.txt", emit: metrics

  script:
  """
  samtools fixmate -@ ${task.cpus} -m ${bam} ${sample_id}.fixmate.bam
  samtools sort -@ ${task.cpus} -o ${sample_id}.positionsort.bam ${sample_id}.fixmate.bam
  samtools markdup -@ ${task.cpus} -s ${sample_id}.positionsort.bam ${sample_id}.dedup.bam 2> ${sample_id}.duplication_metrics.txt
  samtools index -@ ${task.cpus} ${sample_id}.dedup.bam
  """
}

