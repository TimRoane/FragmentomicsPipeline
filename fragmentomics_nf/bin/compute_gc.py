#!/usr/bin/env python
import argparse
import pandas as pd


def gc_for_seq(seq):
    seq = seq.upper()
    denom = sum(1 for b in seq if b in "ACGT")
    return 0.0 if denom == 0 else (seq.count("G") + seq.count("C")) / denom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins", required=True)
    ap.add_argument("--fasta")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    bins = pd.read_csv(args.bins, sep="\t", names=["chrom", "start", "end", "bin_id"])
    if args.fasta:
        try:
            import pysam
            fasta = pysam.FastaFile(args.fasta)
            bins["gc_content"] = [gc_for_seq(fasta.fetch(r.chrom, int(r.start), int(r.end))) for r in bins.itertuples()]
        except Exception:
            bins["gc_content"] = 0.5
    else:
        bins["gc_content"] = 0.5
    bins["effective_bin_size"] = bins["end"] - bins["start"]
    bins.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

