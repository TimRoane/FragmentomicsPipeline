#!/usr/bin/env python
import argparse
import gzip
import hashlib
import sys

import pandas as pd


def is_autosome(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return c.isdigit() and 1 <= int(c) <= 22


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--out-parquet", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--min-len", type=int, default=50)
    ap.add_argument("--max-len", type=int, default=500)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--autosomes-only", action="store_true")
    args = ap.parse_args()

    try:
        import pysam
    except ImportError:
        raise SystemExit("pysam is required for bam_to_fragments.py")

    rows = []
    bam = pysam.AlignmentFile(args.bam)
    for read in bam.fetch(until_eof=True):
        if not read.is_read1:
            continue
        if read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair:
            continue
        if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
            continue
        if read.mapping_quality < args.min_mapq:
            continue
        chrom = bam.get_reference_name(read.reference_id)
        if args.autosomes_only and not is_autosome(chrom):
            continue
        mate_start = read.next_reference_start
        tlen = abs(read.template_length)
        if tlen <= 0:
            continue
        start = min(read.reference_start, mate_start)
        end = start + tlen
        length = end - start
        if length < args.min_len or length > args.max_len:
            continue
        rows.append({
            "sample_id": args.sample_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "mapq": int(read.mapping_quality),
            "strand": "-" if read.is_reverse else "+",
            "length": length,
            "is_short": 100 <= length <= 150,
            "is_long": 151 <= length <= 220,
            "read_name_hash": hashlib.sha1(read.query_name.encode()).hexdigest()[:16],
        })

    df = pd.DataFrame(rows, columns=[
        "sample_id", "chrom", "start", "end", "mapq", "strand", "length",
        "is_short", "is_long", "read_name_hash",
    ])
    df.to_parquet(args.out_parquet, index=False)
    with gzip.open(args.out_bed, "wt") as bed:
        for row in df.itertuples(index=False):
            bed.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.sample_id}:{row.read_name_hash}\t{row.mapq}\t{row.strand}\n")


if __name__ == "__main__":
    main()
