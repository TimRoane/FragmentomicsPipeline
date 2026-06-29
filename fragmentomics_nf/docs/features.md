# Features

The pipeline implements transparent public-literature fragmentomic features.

Core features:

- 100 kb and 5 Mb genomic bin counts.
- Short fragments: 100-150 bp by default.
- Long fragments: 151-220 bp by default.
- Short/long ratios.
- Per-sample LOESS GC correction of total, short, and long coverage.
- Per-sample centering and scaling for selected binned feature groups.
- Chromosome-arm GC-corrected depth summaries and placeholder Z-score columns when a reference panel is absent.
- Mitochondrial read proportion and log10 mitochondrial fraction.

Optional features:

- 4-mer end motifs from both fragment ends.
- Snyder-style WPS over supplied regulatory regions.

This implementation is intended for reproducible research, QC, and transparent baseline modeling.
