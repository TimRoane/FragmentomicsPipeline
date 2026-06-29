#!/usr/bin/env python
import argparse
import json


ACCEPTED_EXTENSIONS = (
    ".frag.tsv.bgz", ".frag.tsv.gz", ".tsv.bgz", ".tsv.gz",
    ".bed.bgz", ".bed.gz", ".bgz", ".gz",
)


def is_accepted_fragment_key(key):
    return str(key).endswith(ACCEPTED_EXTENSIONS)


def find_fragment_key(record, assembly):
    analyses = record.get("analysis", {}).get(assembly, [])
    for item in analyses:
        desc = str(item.get("desc", "")).strip().lower()
        typ = str(item.get("type", "")).strip().lower()
        key = item.get("key", "")
        if desc == "fragment" and typ == "tsv" and is_accepted_fragment_key(key):
            return key
    raise ValueError(f"No FinaleDB fragment TSV found under analysis.{assembly}")


def build_download_url(key):
    return f"http://finaledb.research.cchmc.org/data/{key}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata-json", required=True)
    ap.add_argument("--assembly", default="hg19")
    ap.add_argument("--url", action="store_true")
    args = ap.parse_args()
    with open(args.metadata_json) as fh:
        record = json.load(fh)
    key = find_fragment_key(record, args.assembly)
    print(build_download_url(key) if args.url else key)


if __name__ == "__main__":
    main()

