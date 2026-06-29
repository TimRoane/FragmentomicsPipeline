#!/usr/bin/env python
import argparse
import math
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--bam")
    ap.add_argument("--out", required=True)
    ap.add_argument("--short-min", type=int, default=100)
    ap.add_argument("--short-max", type=int, default=150)
    ap.add_argument("--long-min", type=int, default=151)
    ap.add_argument("--long-max", type=int, default=220)
    args = ap.parse_args()

    df = pd.read_parquet(args.fragments)
    lengths = df["length"] if len(df) else pd.Series(dtype=float)
    n = int(len(df))
    short = int(lengths.between(args.short_min, args.short_max).sum()) if n else 0
    long = int(lengths.between(args.long_min, args.long_max).sum()) if n else 0
    mode = int(lengths.mode().iloc[0]) if n and not lengths.mode().empty else 0
    metrics = {
        "sample_id": args.sample_id,
        "filtered_fragments": n,
        "median_fragment_length": float(lengths.median()) if n else 0,
        "mode_fragment_length": mode,
        "fraction_short_100_150": short / n if n else 0,
        "fraction_long_151_220": long / n if n else 0,
        "short_long_ratio_global": short / max(long, 1),
        "fraction_lt_100": int((lengths < 100).sum()) / n if n else 0,
        "fraction_100_150": short / n if n else 0,
        "fraction_151_220": long / n if n else 0,
        "fraction_221_320": int(lengths.between(221, 320).sum()) / n if n else 0,
        "fraction_gt_320": int((lengths > 320).sum()) / n if n else 0,
    }
    pd.DataFrame([metrics]).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

