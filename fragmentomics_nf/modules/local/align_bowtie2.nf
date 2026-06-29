process ALIGN_BOWTIE2 {
  tag "$sample_id"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.aligned.bam*'

  input:
  tuple val(sample_id), path(read1), path(read2)

  output:
  tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: bam

  script:
  """
  bowtie2 -x ${params.bowtie2_index} -1 ${read1} -2 ${read2} -p ${task.cpus} \
    | samtools sort -@ ${task.cpus} -o ${sample_id}.aligned.bam
  samtools index -@ ${task.cpus} ${sample_id}.aligned.bam
  """
}

