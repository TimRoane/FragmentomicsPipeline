#!/usr/bin/env python
import argparse


def is_autosome(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return c.isdigit() and 1 <= int(c) <= 22


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--bin-size", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--autosomes-only", default="true")
    args = ap.parse_args()
    autosomes_only = str(args.autosomes_only).lower() in {"1", "true", "yes"}
    with open(args.chrom_sizes) as inp, open(args.out, "w") as out:
        for line in inp:
            if not line.strip():
                continue
            chrom, size = line.rstrip().split("\t")[:2]
            size = int(size)
            if autosomes_only and not is_autosome(chrom):
                continue
            for start in range(0, size, args.bin_size):
                end = min(start + args.bin_size, size)
                if end - start >= args.bin_size * 0.5:
                    out.write(f"{chrom}\t{start}\t{end}\t{chrom}:{start}-{end}\n")


if __name__ == "__main__":
    main()

