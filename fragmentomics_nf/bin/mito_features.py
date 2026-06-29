#!/usr/bin/env python
import argparse
import math
import pandas as pd


def is_mito(name):
    return name in {"chrM", "MT", "M", "chrMT"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    try:
        import pysam
        bam = pysam.AlignmentFile(args.bam)
        refs = list(bam.references)
        mito = sum(bam.count(r) for r in refs if is_mito(r))
        nuclear = sum(bam.count(r) for r in refs if not is_mito(r))
    except Exception:
        mito = 0
        nuclear = 0
    total = mito + nuclear
    frac = mito / total if total else 0.0
    pd.DataFrame([{
        "sample_id": args.sample_id,
        "mito_reads": mito,
        "nuclear_reads": nuclear,
        "mito_fraction": frac,
        "log10_mito_fraction": math.log10(max(frac, 1e-12)),
    }]).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

