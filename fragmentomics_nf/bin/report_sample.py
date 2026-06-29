#!/usr/bin/env python
import argparse
import html
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    df = pd.read_csv(args.features, sep="\t")
    row = df[df["sample_id"] == args.sample_id].iloc[0] if "sample_id" in df else df.iloc[0]
    items = "\n".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>" for k, v in row.items())
    Path(args.out).write_text(f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(args.sample_id)} fragmentomics report</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #1f2933; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
    th, td {{ border-bottom: 1px solid #d9e2ec; padding: 0.45rem; text-align: left; }}
    th {{ width: 22rem; }}
  </style>
</head>
<body>
  <h1>{html.escape(args.sample_id)} fragmentomics report</h1>
  <p>Research-use public-literature fragmentomics summary. This report is not a clinical diagnostic interpretation.</p>
  <h2>Sample summary</h2>
  <table>{items}</table>
</body>
</html>
""")


if __name__ == "__main__":
    main()
