process GC_CORRECT_BINS {
  tag "$sample_id"
  publishDir "${params.outdir}/gc_corrected", mode: 'copy', pattern: '*.tsv.gz'

  input:
  tuple val(sample_id), path(counts)

  output:
  tuple val(sample_id), path("${sample_id}.gc_corrected.tsv.gz"), emit: corrected

  script:
  """
  gc_correct.py --counts ${counts} --out ${sample_id}.gc_corrected.tsv.gz --loess-frac ${params.loess_frac}
  """
}

