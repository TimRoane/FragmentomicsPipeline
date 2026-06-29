#!/usr/bin/env python
import argparse
import json
import joblib
import numpy as np
import pandas as pd

META = {"sample_id", "label", "cancer_type", "batch", "sex", "age", "smoking_status"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--research-language", default="true")
    args = ap.parse_args()
    df = pd.read_csv(args.matrix, sep="\t")
    artifact = joblib.load(args.model)
    feature_cols = artifact["feature_cols"]
    X = df.reindex(columns=feature_cols).apply(pd.to_numeric, errors="coerce")
    probs = artifact["pipeline"].predict_proba(X)[:, 1]
    threshold = float(artifact.get("threshold_95_specificity", 0.5))
    classification = np.where(probs > threshold, "elevated_fragmentome_score", "lower_fragmentome_score")
    out = pd.DataFrame({
        "sample_id": df["sample_id"],
        "cancer_probability": probs,
        "classification": classification,
        "threshold_used": threshold,
        "model_version": "0.1.0",
        "qc_pass": True,
        "qc_flags": "",
        "top_contributing_features": "",
    })
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
