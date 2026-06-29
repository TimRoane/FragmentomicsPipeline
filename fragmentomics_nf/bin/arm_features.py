#!/usr/bin/env python
import argparse
import math
import pandas as pd


def read_arms(path):
    if not path:
        return None
    return pd.read_csv(path, sep="\t", names=["chrom", "start", "end", "arm"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--corrected", required=True)
    ap.add_argument("--arms", default="")
    ap.add_argument("--healthy-reference", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    df = pd.read_csv(args.corrected, sep="\t")
    arms = read_arms(args.arms)
    features = {"sample_id": args.sample_id}
    if arms is not None:
        for arm in arms.itertuples(index=False):
            sub = df[(df.chrom == arm.chrom) & (df.start < arm.end) & (df.end > arm.start)]
            depth = float(sub["coverage_total_gc_corrected"].mean()) if len(sub) else 0.0
            name = str(arm.arm).replace("chr", "arm_chr")
            features[f"{name}_gc_corrected_depth"] = depth
            features[f"{name}_z"] = 0.0
    else:
        for chrom, sub in df.groupby("chrom"):
            depth = float(sub["coverage_total_gc_corrected"].mean()) if len(sub) else 0.0
            features[f"arm_{chrom}_z"] = 0.0
            features[f"arm_{chrom}_gc_corrected_depth"] = depth
    z_cols = [v for k, v in features.items() if k.endswith("_z") and isinstance(v, (int, float))]
    features["plasma_aneuploidy_score"] = float(sum(abs(v) for v in z_cols) / max(len(z_cols), 1))
    pd.DataFrame([features]).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

