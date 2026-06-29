#!/usr/bin/env python
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"


class BinFragmentCountsTests(unittest.TestCase):
    def test_vectorized_counter_counts_boundary_overlaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pd.DataFrame({
                "sample_id": ["s1", "s1", "s1", "s1"],
                "chrom": ["chr1", "chr1", "chr1", "chr2"],
                "start": [10, 90, 100, 5],
                "end": [50, 110, 170, 140],
                "mapq": [60, 60, 60, 60],
                "strand": ["+", "+", "-", "+"],
                "length": [40, 20, 70, 135],
                "is_short": [False, False, False, True],
                "is_long": [False, False, False, False],
            }).to_parquet(tmp / "frags.parquet", index=False)
            pd.DataFrame({
                "chrom": ["chr1", "chr1", "chr2"],
                "start": [0, 100, 0],
                "end": [100, 200, 100],
                "bin_id": ["chr1:0-100", "chr1:100-200", "chr2:0-100"],
                "gc_content": [0.5, 0.5, 0.5],
                "effective_bin_size": [100, 100, 100],
            }).to_csv(tmp / "bins_gc.tsv", sep="\t", index=False)
            subprocess.run([
                sys.executable, str(BIN / "bin_fragment_counts.py"),
                "--sample-id", "s1",
                "--fragments", str(tmp / "frags.parquet"),
                "--bins-gc", str(tmp / "bins_gc.tsv"),
                "--out", str(tmp / "counts.tsv.gz"),
                "--short-min", "100",
                "--short-max", "150",
                "--long-min", "151",
                "--long-max", "220",
            ], check=True)
            out = pd.read_csv(tmp / "counts.tsv.gz", sep="\t")
            self.assertEqual(out["n_total_fragments"].tolist(), [2, 2, 1])
            self.assertEqual(out["n_short_100_150"].tolist(), [0, 0, 1])
            self.assertEqual(out["median_fragment_length"].tolist(), [30.0, 45.0, 135.0])


if __name__ == "__main__":
    unittest.main()

