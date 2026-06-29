process BUILD_BINS {
  tag "${params.genome}"
  publishDir "resources/generated", mode: 'copy', pattern: '*.bed'

  input:
  path chrom_sizes

  output:
  path "${params.genome}.${params.bin_small}.autosomes.filtered.bed", emit: bins
  path "${params.genome}.${params.bin_large}.autosomes.filtered.bed", emit: bins_large

  script:
  """
  make_bins.py --chrom-sizes ${chrom_sizes} --bin-size ${params.bin_small} --out ${params.genome}.${params.bin_small}.autosomes.filtered.bed --autosomes-only ${params.autosomes_only}
  make_bins.py --chrom-sizes ${chrom_sizes} --bin-size ${params.bin_large} --out ${params.genome}.${params.bin_large}.autosomes.filtered.bed --autosomes-only ${params.autosomes_only}
  """
}

