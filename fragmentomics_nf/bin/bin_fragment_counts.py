#!/usr/bin/env python
import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def add_overlaps(frags, bins, offsets, counts, short_counts, long_counts, length_sums, length_hist, args):
    if frags.empty or bins.empty:
        return

    starts = bins["start"].to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64)
    global_idx = bins["_global_idx"].to_numpy(dtype=np.int64)

    f_start = frags["start"].to_numpy(dtype=np.int64)
    f_end = frags["end"].to_numpy(dtype=np.int64)
    lengths = frags["length"].to_numpy(dtype=np.int64)

    first = np.searchsorted(ends, f_start, side="right")
    last = np.searchsorted(starts, f_end, side="left")
    span = last - first
    keep = span > 0
    if not keep.any():
        return

    first = first[keep]
    span = span[keep]
    lengths = lengths[keep]

    # Most cfDNA fragments hit one bin; a tiny boundary fraction hits two.
    # Handle those vectorized cases directly and keep a general fallback.
    if np.all(span == 1):
        local_hits = first
        hit_lengths = lengths
    elif np.all(span <= 2):
        second_mask = span == 2
        local_hits = np.concatenate([first, first[second_mask] + 1])
        hit_lengths = np.concatenate([lengths, lengths[second_mask]])
    else:
        local_hits = np.concatenate([np.arange(s, s + n, dtype=np.int64) for s, n in zip(first, span)])
        hit_lengths = np.repeat(lengths, span)
    hit_bins = global_idx[local_hits]

    np.add.at(counts, hit_bins, 1)
    np.add.at(short_counts, hit_bins, ((hit_lengths >= args.short_min) & (hit_lengths <= args.short_max)).astype(np.int64))
    np.add.at(long_counts, hit_bins, ((hit_lengths >= args.long_min) & (hit_lengths <= args.long_max)).astype(np.int64))
    np.add.at(length_sums, hit_bins, hit_lengths)

    hist_lengths = np.clip(hit_lengths, 0, args.hist_max_length)
    np.add.at(length_hist, (hit_bins, hist_lengths), 1)


def medians_from_hist(length_hist, counts):
    med = np.zeros(len(counts), dtype=float)
    nonzero = np.flatnonzero(counts)
    for i in nonzero:
        total = int(counts[i])
        midpoint = (total - 1) / 2
        cumsum = np.cumsum(length_hist[i])
        lower = int(np.searchsorted(cumsum, midpoint, side="right"))
        upper = int(np.searchsorted(cumsum, total / 2, side="right"))
        med[i] = (lower + upper) / 2
    return med


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--bins-gc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--short-min", type=int, default=100)
    ap.add_argument("--short-max", type=int, default=150)
    ap.add_argument("--long-min", type=int, default=151)
    ap.add_argument("--long-max", type=int, default=220)
    ap.add_argument("--batch-size", type=int, default=1_000_000)
    ap.add_argument("--hist-max-length", type=int, default=1000)
    args = ap.parse_args()

    bins = pd.read_csv(args.bins_gc, sep="\t").reset_index(drop=True)
    bins["_global_idx"] = np.arange(len(bins), dtype=np.int64)
    bins_by_chrom = {chrom: sub.sort_values("start") for chrom, sub in bins.groupby("chrom", sort=False)}

    n_bins = len(bins)
    counts = np.zeros(n_bins, dtype=np.int64)
    short_counts = np.zeros(n_bins, dtype=np.int64)
    long_counts = np.zeros(n_bins, dtype=np.int64)
    length_sums = np.zeros(n_bins, dtype=np.float64)
    length_hist = np.zeros((n_bins, args.hist_max_length + 1), dtype=np.int32)

    columns = ["chrom", "start", "end", "length"]
    parquet = pq.ParquetFile(args.fragments)
    for batch in parquet.iter_batches(batch_size=args.batch_size, columns=columns):
        frags = batch.to_pandas()
        for chrom, chrom_frags in frags.groupby("chrom", sort=False):
            chrom_bins = bins_by_chrom.get(chrom)
            if chrom_bins is not None:
                add_overlaps(chrom_frags, chrom_bins, None, counts, short_counts, long_counts, length_sums, length_hist, args)

    medians = medians_from_hist(length_hist, counts)
    eff = bins["effective_bin_size"].clip(lower=1).to_numpy(dtype=float)
    mean_lengths = np.divide(length_sums, counts, out=np.zeros_like(length_sums), where=counts > 0)

    out = bins.drop(columns=["_global_idx"]).copy()
    out.insert(0, "sample_id", args.sample_id)
    out["n_total_fragments"] = counts
    out["n_short_100_150"] = short_counts
    out["n_long_151_220"] = long_counts
    out["short_long_ratio"] = short_counts / np.maximum(long_counts, 1)
    out["mean_fragment_length"] = mean_lengths
    out["median_fragment_length"] = medians
    out["coverage_total"] = counts / eff
    out["coverage_short"] = short_counts / eff
    out["coverage_long"] = long_counts / eff

    ordered = [
        "sample_id", "chrom", "start", "end", "bin_id", "n_total_fragments",
        "n_short_100_150", "n_long_151_220", "short_long_ratio",
        "mean_fragment_length", "median_fragment_length", "coverage_total",
        "coverage_short", "coverage_long", "gc_content", "effective_bin_size",
    ]
    out[ordered].to_csv(args.out, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
