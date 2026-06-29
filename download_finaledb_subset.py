#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import urlopen, Request

BASE = "http://finaledb.research.cchmc.org"
API = f"{BASE}/api/v1/seqrun"


def fetch_json(url: str) -> dict:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def query_finaledb(disease: str, tissue: str, limit: int, offset: int = 0, frag_min=None, frag_max=None) -> list[dict]:
    params = {
        "disease": disease,
        "tissue": tissue,
        "limit": limit,
        "offset": offset,
    }

    # Optional: avoid tiny samples. Good starting range for cfDNA WGS.
    if frag_min is not None and frag_max is not None:
        params["frag_num"] = f"{frag_min},{frag_max}"

    url = f"{API}?{urlencode(params)}"
    print(f"[query] {url}", file=sys.stderr)
    payload = fetch_json(url)
    return payload.get("results", [])


def walk_strings(obj):
    """Yield every string value from a nested dict/list."""
    if isinstance(obj, dict):
        for v in obj.values():
            yield from walk_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk_strings(v)
    elif isinstance(obj, str):
        yield obj


def looks_like_fragment_file(s: str, assembly: str) -> bool:
    lower = s.lower()

    # FinaleDB fragment files are gzipped/indexed TSV-like fragment files.
    has_frag_word = any(x in lower for x in ["frag", "fragment"])
    has_good_ext = any(lower.endswith(x) for x in [".gz", ".tsv.gz", ".bed.gz", ".frag.gz"])
    has_assembly = assembly.lower() in lower

    # Avoid QC reports/bigWigs if possible.
    bad = any(x in lower for x in [".bw", ".bigwig", "fastqc", "multiqc", "gc_bias", "complexity", ".html", ".pdf"])

    return has_frag_word and has_good_ext and not bad and (has_assembly or assembly == "any")


def normalize_download_url(path_or_url: str) -> str:
    if path_or_url.startswith("http://") or path_or_url.startswith("https://"):
        return path_or_url

    # Portal exposes S3 public files through /data/*
    if path_or_url.startswith("/data/"):
        return BASE + path_or_url

    if path_or_url.startswith("data/"):
        return BASE + "/" + path_or_url

    # Fallback. This handles returned paths like hg19/sample/file.frag.gz.
    return BASE + "/data/" + path_or_url.lstrip("/")


def find_fragment_urls(record: dict, assembly: str) -> list[str]:
    """
    Prefer the structured FinaleDB analysis block:
      analysis.hg19[] / analysis.hg38[]
      item.desc == "fragment"
      item.key == entries/EE.../hg19/EE....frag.tsv.bgz
    """
    urls = []

    analysis = record.get("analysis") or {}

    if assembly == "any":
        assemblies = list(analysis.keys())
    else:
        assemblies = [assembly]

    for asm in assemblies:
        for item in analysis.get(asm, []):
            desc = str(item.get("desc", "")).strip().lower()
            key = item.get("key")
            typ = str(item.get("type", "")).strip().lower()

            # Important: avoid "fragment profile"; we only want raw fragments.
            if desc == "fragment" and key and typ in {"tsv", "bed"}:
                urls.append(normalize_download_url(key))

    # Fallback string scan, now including .bgz.
    if not urls:
        for s in walk_strings(record):
            lower = s.lower()
            if (
                ".frag." in lower
                and lower.endswith((".bgz", ".gz", ".tsv.bgz", ".tsv.gz", ".bed.bgz", ".bed.gz"))
                and ("bedgraph" not in lower)
                and ("frag_profile" not in lower)
            ):
                urls.append(normalize_download_url(s))

    # Deduplicate.
    seen = set()
    out = []
    for u in urls:
        if u not in seen:
            seen.add(u)
            out.append(u)

    return out

def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")


def download(url: str, out_path: Path, dry_run: bool = False):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        print(f"[dry-run] {url} -> {out_path}")
        return

    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"[skip] {out_path} exists")
        return

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry", "5",
        "--retry-delay", "5",
        "-o", str(out_path),
        url,
    ]
    print("[download]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="finaledb_subset")
    p.add_argument("--assembly", default="hg19", choices=["hg19", "hg38", "any"])
    p.add_argument("--n-controls", type=int, default=25)
    p.add_argument("--n-lung", type=int, default=25)
    p.add_argument("--tissue", default="blood plasma")
    p.add_argument("--frag-min", type=int, default=5_000_000)
    p.add_argument("--frag-max", type=int, default=100_000_000)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    outdir = Path(args.outdir)
    manifest_path = outdir / "manifest.tsv"
    raw_json_dir = outdir / "metadata_json"
    raw_json_dir.mkdir(parents=True, exist_ok=True)

    cohorts = [
        ("control", "Healthy", args.n_controls),
        ("lung_cancer", "Lung cancer", args.n_lung),
    ]

    manifest_rows = []

    for cohort, disease, n in cohorts:
        records = query_finaledb(
            disease=disease,
            tissue=args.tissue,
            limit=n,
            frag_min=args.frag_min,
            frag_max=args.frag_max,
        )

        print(f"[info] {disease}: got {len(records)} records", file=sys.stderr)

        for rec in records:
            entry_id = rec.get("id")
            finaledb_id = f"EE{entry_id}" if entry_id is not None else "unknown"

            sample = rec.get("sample") or {}
            sample_name = sample.get("name") or finaledb_id
            sample_disease = sample.get("disease") or disease
            sample_tissue = sample.get("tissue") or args.tissue

            raw_json_path = raw_json_dir / f"{safe_name(finaledb_id)}.json"
            raw_json_path.write_text(json.dumps(rec, indent=2))

            urls = find_fragment_urls(rec, args.assembly)

            if not urls:
                print(f"[warn] no fragment URL found for {finaledb_id}. Inspect {raw_json_path}", file=sys.stderr)
                continue

            # Use the first best-looking fragment file.
            url = urls[0]

            ext = ".frag.tsv.gz"
            if url.endswith(".bed.gz"):
                ext = ".bed.gz"
            elif url.endswith(".frag.gz"):
                ext = ".frag.gz"
            elif url.endswith(".tsv.gz"):
                ext = ".tsv.gz"

            out_path = outdir / "fragments" / cohort / f"{safe_name(finaledb_id)}_{safe_name(sample_name)}_{args.assembly}{ext}"
            download(url, out_path, dry_run=args.dry_run)

            # Try to download index too if it follows normal tabix naming.
            index_url = url + ".tbi"
            index_path = Path(str(out_path) + ".tbi")
            try:
                download(index_url, index_path, dry_run=args.dry_run)
            except subprocess.CalledProcessError:
                print(f"[warn] no index found at {index_url}", file=sys.stderr)

            manifest_rows.append({
                "cohort": cohort,
                "finaledb_id": finaledb_id,
                "sample_name": sample_name,
                "disease": sample_disease,
                "tissue": sample_tissue,
                "assembly": args.assembly,
                "fragment_url": url,
                "local_fragment": str(out_path),
                "metadata_json": str(raw_json_path),
            })

    outdir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "cohort",
            "finaledb_id",
            "sample_name",
            "disease",
            "tissue",
            "assembly",
            "fragment_url",
            "local_fragment",
            "metadata_json",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"[done] wrote {manifest_path}")
    print(f"[done] downloaded or selected {len(manifest_rows)} fragment files")


if __name__ == "__main__":
    main()
