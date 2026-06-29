#!/usr/bin/env python
import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from finaledb_metadata import build_download_url, find_fragment_key, is_accepted_fragment_key

REQUIRED_COLUMNS = [
    "sample_id", "fastq_1", "fastq_2", "bam", "cram", "label", "cancer_type",
    "batch", "sex", "age", "smoking_status",
]
OPTIONAL_COLUMNS = [
    "fragment", "fragments", "fragment_file", "finaledb_fragment", "metadata_json", "finaledb_url",
    "local_fragment", "fragment_url", "finaledb_id", "sample_name", "cohort", "disease", "tissue", "assembly",
]


def present(value):
    return value is not None and str(value).strip() not in {"", "NA", "na", "null", "None"}


def resolve_local_path(value, base_dir):
    if not present(value):
        return ""
    value = str(value).strip()
    if value.startswith(("http://", "https://", "entries/")):
        return value
    path = Path(value)
    if path.is_absolute():
        return str(path)
    if path.exists():
        return str(path.resolve())
    based = Path(base_dir) / value
    return str(based.resolve()) if based.exists() else value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--mode", default="full")
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--assembly")
    ap.add_argument("--input-type", default="auto")
    ap.add_argument("--path-base", default=".")
    ap.add_argument("--fasta")
    ap.add_argument("--tissue-of-origin", action="store_true")
    args = ap.parse_args()

    with open(args.samplesheet, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        rows = list(csv.DictReader(fh, dialect=dialect))
    if not rows:
        raise SystemExit("samplesheet is empty")
    required = [] if args.input_type == "finaledb_fragments" else REQUIRED_COLUMNS
    missing_cols = [c for c in required if c not in rows[0]]
    if missing_cols:
        raise SystemExit(f"samplesheet missing required columns: {', '.join(missing_cols)}")

    seen = set()
    errors = []
    modes = {"fastq": 0, "bam": 0, "cram": 0, "finaledb_fragments": 0}
    for i, row in enumerate(rows, start=2):
        if args.input_type == "finaledb_fragments":
            row.setdefault("sample_id", row.get("finaledb_id") or row.get("sample_name") or "")
            row.setdefault("label", row.get("cohort") or row.get("disease") or "")
            if not row.get("cancer_type", "").strip() and row.get("cohort", "").strip() not in {"", "control"}:
                row["cancer_type"] = row.get("cohort", "")
        sid = row["sample_id"].strip()
        if not sid:
            errors.append(f"line {i}: sample_id is required")
        if sid in seen:
            errors.append(f"line {i}: duplicate sample_id '{sid}'")
        seen.add(sid)

        if args.input_type == "finaledb_fragments":
            direct = next((row.get(c, "").strip() for c in ["finaledb_fragment", "local_fragment", "fragment", "fragments", "fragment_file"] if present(row.get(c))), "")
            direct = resolve_local_path(direct, args.path_base)
            metadata_json = resolve_local_path(row.get("metadata_json", "").strip(), args.path_base)
            if direct:
                if direct.startswith("entries/"):
                    row["finaledb_url"] = build_download_url(direct)
                if not row.get("finaledb_url", "").strip() and row.get("fragment_url", "").strip():
                    row["finaledb_url"] = row["fragment_url"].strip()
                if not (direct.startswith("http://") or direct.startswith("https://")) and not is_accepted_fragment_key(direct):
                    errors.append(f"line {i}: unsupported FinaleDB fragment extension '{direct}'")
                row["finaledb_fragment"] = direct
                if row.get("local_fragment"):
                    row["local_fragment"] = direct
            if metadata_json and Path(metadata_json).exists():
                try:
                    with open(metadata_json) as fh:
                        key = find_fragment_key(json.load(fh), args.assembly or args.genome)
                    row["finaledb_fragment"] = direct or key
                    row["finaledb_url"] = build_download_url(key)
                    row["metadata_json"] = metadata_json
                except Exception as exc:
                    errors.append(f"line {i}: could not parse FinaleDB metadata_json: {exc}")
            elif metadata_json and not direct:
                errors.append(f"line {i}: could not parse FinaleDB metadata_json: no such file '{metadata_json}'")
            elif not direct:
                errors.append(f"line {i}: finaledb_fragments input requires finaledb_fragment, fragment_file, fragment, fragments, or metadata_json")
            modes["finaledb_fragments"] += 1
        else:
            has_fastq = present(row["fastq_1"]) or present(row["fastq_2"])
            has_fastq_pair = present(row["fastq_1"]) and present(row["fastq_2"])
            has_bam = present(row["bam"])
            has_cram = present(row["cram"])
            source_count = int(has_fastq_pair) + int(has_bam) + int(has_cram)
            if source_count != 1:
                errors.append(f"line {i}: provide exactly one input type: paired FASTQ, BAM, or CRAM")
            if has_fastq and not has_fastq_pair:
                errors.append(f"line {i}: FASTQ input requires fastq_1 and fastq_2")
            if has_cram and not present(args.fasta):
                errors.append(f"line {i}: CRAM input requires --fasta")
            if has_fastq_pair:
                modes["fastq"] += 1
            if has_bam:
                modes["bam"] += 1
            if has_cram:
                modes["cram"] += 1
        if args.mode == "train" and not present(row["label"]):
            errors.append(f"line {i}: label is required in train mode")
        if args.tissue_of_origin and not present(row["cancer_type"]):
            errors.append(f"line {i}: cancer_type is required for tissue-of-origin training")

    if errors:
        raise SystemExit("\n".join(errors))

    with open(args.out, "w", newline="") as out:
        fieldnames = REQUIRED_COLUMNS + [c for c in OPTIONAL_COLUMNS if any(c in r for r in rows)]
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "").strip() for c in fieldnames})

    manifest = {
        "pipeline": "fragmentomics_nf",
        "analysis_type": "public-literature fragmentomics",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode": args.mode,
        "input_type": args.input_type,
        "genome": args.genome,
        "assembly": args.assembly or args.genome,
        "sample_count": len(rows),
        "input_counts": modes,
        "samplesheet": str(Path(args.samplesheet).resolve()),
    }
    Path(args.manifest).write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
