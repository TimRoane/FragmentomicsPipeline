process ALIGN_BWAMEM2 {
  tag "$sample_id"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.aligned.bam*'

  input:
  tuple val(sample_id), path(read1), path(read2)

  output:
  tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: bam

  script:
  """
  bwa-mem2 mem -t ${task.cpus} ${params.bwa_index ?: params.fasta} ${read1} ${read2} \
    | samtools sort -@ ${task.cpus} -o ${sample_id}.aligned.bam
  samtools index -@ ${task.cpus} ${sample_id}.aligned.bam
  """
}

