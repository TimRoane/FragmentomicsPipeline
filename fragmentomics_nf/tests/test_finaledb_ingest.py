#!/usr/bin/env python
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from finaledb_metadata import build_download_url, find_fragment_key, is_accepted_fragment_key
from ingest_finaledb_fragments import normalize_chrom


class FinaleDbIngestTests(unittest.TestCase):
    def test_metadata_parser_selects_fragment_not_profile(self):
        record = {
            "analysis": {
                "hg19": [
                    {"desc": "fragment profile", "type": "tsv", "key": "entries/x/hg19/profile.tsv.bgz"},
                    {"desc": "coverage", "type": "bigWig", "key": "entries/x/hg19/cov.bw"},
                    {"desc": "fragment", "type": "tsv", "key": "entries/EE85756/hg19/EE85756.hg19.frag.tsv.bgz"},
                ]
            }
        }
        key = find_fragment_key(record, "hg19")
        self.assertEqual(key, "entries/EE85756/hg19/EE85756.hg19.frag.tsv.bgz")
        self.assertEqual(
            build_download_url(key),
            "http://finaledb.research.cchmc.org/data/entries/EE85756/hg19/EE85756.hg19.frag.tsv.bgz",
        )

    def test_bgz_extension_is_accepted(self):
        self.assertTrue(is_accepted_fragment_key("sample.frag.tsv.bgz"))
        self.assertTrue(is_accepted_fragment_key("sample.bgz"))
        self.assertFalse(is_accepted_fragment_key("sample.bigWig"))

    def test_chrom_normalization(self):
        self.assertEqual(normalize_chrom("1", "chr"), "chr1")
        self.assertEqual(normalize_chrom("MT", "chr"), "chrM")
        self.assertEqual(normalize_chrom("chr1", "no_chr"), "1")
        self.assertEqual(normalize_chrom("chrM", "no_chr"), "MT")

    def test_ingest_filters_mapq_and_computes_length_from_end_start(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            frag = tmp / "tiny.frag.tsv.bgz"
            with gzip.open(frag, "wt") as fh:
                fh.write("1\t10000\t10168\t0\t+\n")
                fh.write("1\t10032\t10175\t23\t-\n")
                fh.write("1\t10033\t10215\t35\t-\n")
                fh.write("2\t20000\t20130\t60\t+\n")
            bins = tmp / "bins.bed"
            bins.write_text("chr1\t0\t100000\tbin1\nchr2\t0\t100000\tbin2\n")
            out = tmp / "out.parquet"
            bed = tmp / "out.bed.gz"
            qc = tmp / "qc.tsv"
            subprocess.run([
                sys.executable, str(BIN / "ingest_finaledb_fragments.py"),
                "--sample-id", "s1",
                "--fragments", str(frag),
                "--out-parquet", str(out),
                "--out-bed", str(bed),
                "--qc-out", str(qc),
                "--min-mapq", "30",
                "--bins", str(bins),
            ], check=True)
            df = pd.read_parquet(out)
            self.assertEqual(df["mapq"].tolist(), [35, 60])
            self.assertEqual(df["length"].tolist(), [182, 130])
            self.assertEqual(df["is_short"].tolist(), [False, True])
            self.assertEqual(df["chrom"].tolist(), ["chr1", "chr2"])
            q = pd.read_csv(qc, sep="\t").iloc[0]
            self.assertEqual(int(q["input_fragments"]), 4)
            self.assertEqual(int(q["fragments_passing_mapq"]), 2)
            self.assertEqual(int(q["fragments_passing_mapq_and_length"]), 2)

    def test_validate_samplesheet_accepts_finaledb_manifest_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            frag = tmp / "sample.frag.tsv.bgz"
            frag.write_bytes(b"")
            meta = tmp / "sample.json"
            meta.write_text(json.dumps({
                "analysis": {
                    "hg19": [
                        {"desc": "fragment profile", "type": "tsv", "key": "entries/x/profile.tsv.bgz"},
                        {"desc": "fragment", "type": "tsv", "key": "entries/EE1/hg19/EE1.hg19.frag.tsv.bgz"},
                    ]
                }
            }))
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "cohort\tfinaledb_id\tsample_name\tdisease\ttissue\tassembly\tfragment_url\tlocal_fragment\tmetadata_json\n"
                f"control\tEE1\tC1\tHealthy\tblood plasma\thg19\thttp://example.invalid/x\t{frag}\t{meta}\n"
            )
            out = tmp / "validated.csv"
            run_manifest = tmp / "run_manifest.json"
            subprocess.run([
                sys.executable, str(BIN / "validate_samplesheet.py"),
                "--samplesheet", str(manifest),
                "--out", str(out),
                "--manifest", str(run_manifest),
                "--mode", "feature_extract",
                "--genome", "hg19",
                "--assembly", "hg19",
                "--input-type", "finaledb_fragments",
            ], check=True)
            df = pd.read_csv(out)
            self.assertEqual(df.loc[0, "sample_id"], "EE1")
            self.assertEqual(df.loc[0, "label"], "control")
            self.assertEqual(df.loc[0, "finaledb_fragment"], str(frag))
            self.assertTrue(df.loc[0, "finaledb_url"].endswith("/entries/EE1/hg19/EE1.hg19.frag.tsv.bgz"))


if __name__ == "__main__":
    unittest.main()
