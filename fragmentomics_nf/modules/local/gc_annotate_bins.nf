process GC_ANNOTATE_BINS {
  tag 'gc'
  publishDir "resources/generated", mode: 'copy', pattern: '*.bin_gc.tsv'

  input:
  path bins

  output:
  path '*.bin_gc.tsv', emit: annotated

  script:
  """
  compute_gc.py --bins ${bins} --fasta ${params.fasta} --out ${bins.simpleName}.bin_gc.tsv
  """
}

