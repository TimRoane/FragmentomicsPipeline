#!/usr/bin/env python
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--regions", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=120)
    args = ap.parse_args()
    frags = pd.read_parquet(args.fragments)
    if args.regions:
        regions = pd.read_csv(args.regions, sep="\t", names=["chrom", "start", "end", "name"])
    else:
        regions = frags[["chrom", "start", "end"]].head(100).copy()
        regions["name"] = [f"region_{i}" for i in range(len(regions))]
    rows = []
    half = args.window // 2
    for r in regions.itertuples(index=False):
        center = (int(r.start) + int(r.end)) // 2
        wstart, wend = center - half, center + half
        sub = frags[frags.chrom == r.chrom]
        spanning = int(((sub.start <= wstart) & (sub.end >= wend)).sum())
        ending = int(((sub.start.between(wstart, wend)) | (sub.end.between(wstart, wend))).sum())
        rows.append({"sample_id": args.sample_id, "chrom": r.chrom, "start": wstart, "end": wend, "name": r.name, "wps": spanning - ending})
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()

