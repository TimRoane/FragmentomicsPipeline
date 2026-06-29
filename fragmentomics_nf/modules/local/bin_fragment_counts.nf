process BIN_FRAGMENT_COUNTS {
  tag "$sample_id"
  publishDir "${params.outdir}/bin_counts", mode: 'copy', pattern: '*.tsv.gz'

  input:
  tuple val(sample_id), path(fragments)
  path bins_gc

  output:
  tuple val(sample_id), path("${sample_id}.${bins_gc.simpleName}.counts.tsv.gz"), emit: counts

  script:
  """
  bin_fragment_counts.py \
    --sample-id ${sample_id} \
    --fragments ${fragments} \
    --bins-gc ${bins_gc} \
    --out ${sample_id}.${bins_gc.simpleName}.counts.tsv.gz \
    --short-min ${params.short_min} --short-max ${params.short_max} \
    --long-min ${params.long_min} --long-max ${params.long_max} \
    --batch-size ${params.bin_count_batch_size} \
    --hist-max-length ${params.bin_count_hist_max_length}
  """
}
