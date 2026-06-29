#!/usr/bin/env python
import argparse
import glob
import json
from pathlib import Path
import pandas as pd

META = ["sample_id", "label", "cancer_type", "batch", "sex", "age", "smoking_status"]


def wide_from_corrected(path):
    df = pd.read_csv(path, sep="\t")
    sid = str(df["sample_id"].iloc[0])
    out = {"sample_id": sid}
    for row in df.itertuples(index=False):
        bid = str(row.bin_id).replace(":", "_").replace("-", "_")
        out[f"bin_{bid}_total_gc_corrected"] = row.coverage_total_gc_corrected
        out[f"bin_{bid}_short_gc_corrected"] = row.coverage_short_gc_corrected
        out[f"bin_{bid}_long_gc_corrected"] = row.coverage_long_gc_corrected
        out[f"bin_{bid}_short_long_ratio"] = row.short_long_ratio
        out[f"bin_{bid}_mean_fragment_length"] = row.mean_fragment_length
    # Per-sample center/scale selected binned groups.
    for prefix in ["total_gc_corrected", "short_gc_corrected", "short_long_ratio"]:
        cols = [c for c in out if c.endswith(prefix)]
        vals = pd.Series([out[c] for c in cols], dtype=float)
        sd = vals.std() or 1.0
        for c, v in zip(cols, vals):
            out[f"{c}_sample_z"] = (v - vals.mean()) / sd
    return out


def read_one(path):
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t")
    return df.iloc[0].to_dict() if len(df) else {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--corrected", nargs="+", required=True)
    ap.add_argument("--arms", nargs="*", default=[])
    ap.add_argument("--mito", nargs="*", default=[])
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-dictionary", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()
    samples = pd.read_csv(args.samplesheet)
    by_sample = {r["sample_id"]: {k: r.get(k, "") for k in META if k in r} for _, r in samples.iterrows()}
    feature_rows = {}
    for p in args.corrected:
        row = wide_from_corrected(p)
        feature_rows[row["sample_id"]] = row
    for group_paths in [args.arms, args.mito]:
        for p in group_paths:
            row = read_one(p)
            sid = row.pop("sample_id", None)
            if sid:
                feature_rows.setdefault(sid, {"sample_id": sid}).update(row)
    rows = []
    seen = set()
    for sid in samples["sample_id"].astype(str):
        feats = feature_rows.get(sid)
        if feats is None:
            continue
        row = by_sample.get(sid, {"sample_id": sid})
        row.update(feats)
        rows.append(row)
        seen.add(sid)
    for sid in sorted(set(feature_rows) - seen):
        row = by_sample.get(sid, {"sample_id": sid})
        row.update(feature_rows[sid])
        rows.append(row)
    matrix = pd.DataFrame(rows)
    meta_cols = [c for c in META if c in matrix.columns]
    feature_cols = sorted([c for c in matrix.columns if c not in meta_cols])
    matrix = matrix[meta_cols + feature_cols]
    matrix.to_csv(args.out_matrix, sep="\t", index=False)
    dictionary = []
    for c in feature_cols:
        group = c.split("_")[0]
        dictionary.append({
            "feature_name": c,
            "feature_group": group,
            "description": "Transparent public-literature fragmentomics feature",
            "genomic_region": "bin_or_summary",
            "normalization": "GC corrected and/or sample standardized where indicated",
            "requires_reference_panel": c.endswith("_z"),
        })
    pd.DataFrame(dictionary).to_csv(args.out_dictionary, sep="\t", index=False)
    Path(args.out_summary).write_text(json.dumps({
        "n_samples": int(len(matrix)),
        "n_features": int(len(feature_cols)),
        "feature_groups": sorted(set(d["feature_group"] for d in dictionary)),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
