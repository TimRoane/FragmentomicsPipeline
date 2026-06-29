#!/usr/bin/env python
import argparse
from collections import Counter
import math
import pandas as pd


def entropy(counts):
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * math.log2(c / total) for c in counts.values() if c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    try:
        import pysam
        fasta = pysam.FastaFile(args.fasta)
    except Exception:
        fasta = None
    frags = pd.read_parquet(args.fragments)
    counts = Counter()
    if fasta is not None:
        for row in frags.itertuples(index=False):
            left = fasta.fetch(row.chrom, int(row.start), int(row.start) + 4).upper()
            right = fasta.fetch(row.chrom, max(int(row.end) - 4, int(row.start)), int(row.end)).upper()
            if len(left) == 4:
                counts[f"left_4mer_{left}"] += 1
            if len(right) == 4:
                counts[f"right_4mer_{right}"] += 1
            if len(left) == 4 and len(right) == 4:
                counts[f"combined_4mer_{left}_{right}"] += 1
    total = max(sum(counts.values()), 1)
    out = {"sample_id": args.sample_id, "motif_diversity_score": entropy(counts)}
    out.update({k: v / total for k, v in counts.items()})
    pd.DataFrame([out]).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

