# Data Sources

Required reference inputs:

- Reference FASTA and `.fai` index for alignment, CRAM decoding, GC annotation, and end motifs.
- Chromosome sizes file with `chrom` and `size` columns.

Recommended BED inputs:

- ENCODE blacklist regions for the selected genome build.
- Low-mappability regions matched to read length and genome build.
- Centromere and telomere regions.
- Chromosome arm definitions.
- Optional regulatory regions such as TSS, DHS/ATAC, TFBS, or CpG islands.

Use genome-build-consistent files. Do not mix hg19 and hg38 resources.

## FinaleDB Fragment Ingest

FinaleDB fragment files are five-column, headerless, BED-like files:

```text
chrom    start    end    mapq    strand
```

The fourth column is MAPQ, not fragment length. The pipeline derives fragment length as `end - start` and applies MAPQ and length filters before binning. Real FinaleDB files may contain MAPQ 0 rows, so `--min_mapq 30` remains the default.

FinaleDB metadata JSON points to fragment files under `analysis.{assembly}[]`. The parser selects the item where `desc == "fragment"` and `type == "tsv"`, then uses `key` to build:

```text
http://finaledb.research.cchmc.org/data/{key}
```

Valid fragment inputs include `.frag.tsv.bgz`, `.frag.tsv.gz`, `.tsv.bgz`, `.tsv.gz`, `.bed.bgz`, `.bed.gz`, `.bgz`, and `.gz`.

Example:

```bash
nextflow run main.nf \
  -profile docker \
  --mode feature_extract \
  --input finaledb_lung_demo_50/nextflow_samplesheet.csv \
  --input_type finaledb_fragments \
  --assembly hg19 \
  --min_mapq 30 \
  --outdir results/finaledb_lung_demo_50
```
