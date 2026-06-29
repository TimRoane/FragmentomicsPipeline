#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd


def correct(y, gc, frac):
    y = np.asarray(y, dtype=float)
    gc = np.asarray(gc, dtype=float)
    if len(y) < 5 or np.nanstd(gc) == 0:
        pred = np.repeat(np.nanmedian(y), len(y))
    else:
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            fit = lowess(y, gc, frac=frac, return_sorted=True)
            pred = np.interp(gc, fit[:, 0], fit[:, 1])
        except Exception:
            coeff = np.polyfit(gc, y, deg=1)
            pred = np.polyval(coeff, gc)
    return y - pred + np.nanmedian(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--loess-frac", type=float, default=0.30)
    args = ap.parse_args()
    df = pd.read_csv(args.counts, sep="\t")
    for col in ["coverage_total", "coverage_short", "coverage_long"]:
        df[f"{col}_gc_corrected"] = correct(df[col], df["gc_content"], args.loess_frac)
    df.to_csv(args.out, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()

