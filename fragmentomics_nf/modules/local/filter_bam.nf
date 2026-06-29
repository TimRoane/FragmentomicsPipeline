process FILTER_BAM {
  tag "$sample_id"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.filtered.bam*'

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.filtered.bam"), emit: bam

  script:
  def regions = [params.blacklist_bed, params.low_mappability_bed, params.centromere_telomere_bed].findAll { it }
  def exclude = regions ? regions.collect { "-L ^${it}" }.join(' ') : ''
  """
  samtools view -@ ${task.cpus} -b -f 2 -F ${params.samtools_exclude_flags} -q ${params.min_mapq} ${exclude} ${bam} \
    | samtools sort -@ ${task.cpus} -o ${sample_id}.filtered.bam
  samtools index -@ ${task.cpus} ${sample_id}.filtered.bam
  """
}

