#!/usr/bin/env python
import argparse
import gzip
from collections import Counter

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ACCEPTED_EXTENSIONS = (
    ".frag.tsv.bgz", ".frag.tsv.gz", ".tsv.bgz", ".tsv.gz",
    ".bed.bgz", ".bed.gz", ".bgz", ".gz",
)


def accepted_fragment_path(path):
    return str(path).endswith(ACCEPTED_EXTENSIONS)


def open_text(path):
    path = str(path)
    if path.endswith((".bgz", ".gz")):
        return gzip.open(path, "rt")
    return open(path, "rt")


def normalize_chrom(chrom, style):
    c = str(chrom)
    raw = c[3:] if c.startswith("chr") else c
    if raw in {"M", "MT"}:
        raw = "MT"
    if style == "chr":
        return "chrM" if raw == "MT" else f"chr{raw}"
    if style == "no_chr":
        return raw
    return c


def infer_style_from_bins(path):
    if not path:
        return "no_chr"
    with open_text(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                chrom = line.split("\t", 1)[0]
                return "chr" if chrom.startswith("chr") else "no_chr"
    return "no_chr"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--out-parquet", required=True)
    ap.add_argument("--out-bed")
    ap.add_argument("--qc-out", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-fragment-length", type=int, default=50)
    ap.add_argument("--max-fragment-length", type=int, default=500)
    ap.add_argument("--short-min", type=int, default=100)
    ap.add_argument("--short-max", type=int, default=150)
    ap.add_argument("--long-min", type=int, default=151)
    ap.add_argument("--long-max", type=int, default=220)
    ap.add_argument("--chrom-style", choices=["match_bins", "chr", "no_chr"], default="match_bins")
    ap.add_argument("--bins", default="")
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    args = ap.parse_args()

    if not accepted_fragment_path(args.fragments):
        raise SystemExit(f"Unsupported FinaleDB fragment extension: {args.fragments}")

    chrom_style = infer_style_from_bins(args.bins) if args.chrom_style == "match_bins" else args.chrom_style
    reader = pd.read_csv(
        open_text(args.fragments),
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "mapq", "strand"],
        usecols=[0, 1, 2, 3, 4],
        dtype={"chrom": "string", "start": "int64", "end": "int64", "mapq": "int64", "strand": "string"},
        chunksize=args.chunksize,
    )
    schema = pa.schema([
        ("sample_id", pa.string()),
        ("chrom", pa.string()),
        ("start", pa.int64()),
        ("end", pa.int64()),
        ("mapq", pa.int64()),
        ("strand", pa.string()),
        ("length", pa.int64()),
        ("is_short", pa.bool_()),
        ("is_long", pa.bool_()),
    ])
    writer = None
    input_count = 0
    mapq_count = 0
    pass_count = 0
    short_count = 0
    long_count = 0
    length_sum = 0.0
    length_hist = Counter()
    bed = gzip.open(args.out_bed, "wt") if args.out_bed else None
    try:
        for df in reader:
            input_count += int(len(df))
            df["length"] = df["end"] - df["start"]
            mapq_df = df[df["mapq"] >= args.min_mapq].copy()
            mapq_count += int(len(mapq_df))
            length_df = mapq_df[
                (mapq_df["length"] >= args.min_fragment_length)
                & (mapq_df["length"] <= args.max_fragment_length)
            ].copy()
            if length_df.empty:
                continue
            length_df["chrom"] = length_df["chrom"].map(lambda c: normalize_chrom(c, chrom_style))
            length_df.insert(0, "sample_id", args.sample_id)
            length_df["is_short"] = length_df["length"].between(args.short_min, args.short_max)
            length_df["is_long"] = length_df["length"].between(args.long_min, args.long_max)
            length_df = length_df[["sample_id", "chrom", "start", "end", "mapq", "strand", "length", "is_short", "is_long"]]
            pass_count += int(len(length_df))
            short_count += int(length_df["is_short"].sum())
            long_count += int(length_df["is_long"].sum())
            length_sum += float(length_df["length"].sum())
            length_hist.update(length_df["length"].astype(int).tolist())
            table = pa.Table.from_pandas(length_df, schema=schema, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(args.out_parquet, schema)
            writer.write_table(table)
            if bed is not None:
                for row in length_df.itertuples(index=False):
                    bed.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.sample_id}\t{row.mapq}\t{row.strand}\n")
    finally:
        if bed is not None:
            bed.close()
    if writer is None:
        pq.write_table(pa.Table.from_pydict({name: [] for name in schema.names}, schema=schema), args.out_parquet)
    else:
        writer.close()

    n = pass_count
    median_length = 0.0
    if n:
        midpoint = (n - 1) / 2
        cumulative = 0
        lower = upper = None
        for length in sorted(length_hist):
            count = length_hist[length]
            prev = cumulative
            cumulative += count
            if lower is None and prev <= midpoint < cumulative:
                lower = length
            if cumulative > n / 2:
                upper = length
                break
        median_length = float((lower + upper) / 2) if lower is not None and upper is not None else 0.0
    qc = pd.DataFrame([{
        "sample_id": args.sample_id,
        "input_fragments": input_count,
        "fragments_passing_mapq": mapq_count,
        "fragments_passing_mapq_and_length": pass_count,
        "mean_length": length_sum / n if n else 0.0,
        "median_length": median_length,
        "fraction_short_100_150": short_count / n if n else 0.0,
        "fraction_long_151_220": long_count / n if n else 0.0,
        "short_long_ratio": short_count / max(long_count, 1),
        "min_mapq_used": args.min_mapq,
        "min_fragment_length_used": args.min_fragment_length,
        "max_fragment_length_used": args.max_fragment_length,
    }])
    qc.to_csv(args.qc_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
