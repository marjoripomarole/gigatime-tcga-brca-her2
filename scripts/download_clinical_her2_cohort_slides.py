#!/usr/bin/env python3
"""Download selected clinical HER2 cohort slides from GDC."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from gdc_query_tcga_brca import download_file


def str_to_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide-table", default="data/tcga_brca/clinical_her2_cohort_slides_files.csv")
    parser.add_argument("--out-dir", default="data/tcga_brca/slides")
    parser.add_argument("--status-out", default="data/tcga_brca/clinical_her2_cohort_slide_download_status.json")
    parser.add_argument("--only-missing", action="store_true", help="Skip slides that already exist locally.")
    parser.add_argument("--max-downloads", type=int, default=0, help="Maximum number of slides to download. Use 0 for all.")
    parser.add_argument("--download-attempts", type=int, default=8)
    return parser.parse_args()


def read_slide_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"case_submitter_id", "slide_file_id", "slide_file_name", "slide_file_size"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def destination_for(row: dict[str, str], out_dir: Path) -> Path:
    return out_dir / row["case_submitter_id"] / row["slide_file_name"]


def write_status(path: Path, records: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "total_records": len(records),
        "downloaded": sum(1 for record in records if record["status"] == "downloaded"),
        "already_present": sum(1 for record in records if record["status"] == "already_present"),
        "failed": sum(1 for record in records if record["status"] == "failed"),
        "records": records,
    }
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    slide_rows = read_slide_rows(Path(args.slide_table))
    out_dir = Path(args.out_dir)
    selected_rows = []
    for row in slide_rows:
        destination = destination_for(row, out_dir)
        if args.only_missing and destination.exists() and destination.stat().st_size > 0:
            continue
        selected_rows.append(row)
    if args.max_downloads:
        selected_rows = selected_rows[: args.max_downloads]

    records: list[dict[str, object]] = []
    for index, row in enumerate(selected_rows, start=1):
        destination = destination_for(row, out_dir)
        expected_size = int(float(row["slide_file_size"])) if row.get("slide_file_size") else None
        record = {
            "case_submitter_id": row["case_submitter_id"],
            "cohort_group": row.get("cohort_group", ""),
            "slide_file_id": row["slide_file_id"],
            "slide_file_name": row["slide_file_name"],
            "destination": str(destination),
            "expected_size": expected_size,
            "status": "",
            "actual_size": 0,
            "error": "",
        }
        if destination.exists() and expected_size and destination.stat().st_size == expected_size:
            record["status"] = "already_present"
            record["actual_size"] = destination.stat().st_size
            records.append(record)
            print(f"[{index}/{len(selected_rows)}] Already present {destination}")
            write_status(Path(args.status_out), records)
            continue
        try:
            print(f"[{index}/{len(selected_rows)}] Downloading {row['case_submitter_id']} {row['slide_file_name']}")
            downloaded = download_file(
                row["slide_file_id"],
                destination,
                expected_size=expected_size,
                max_attempts=args.download_attempts,
            )
            record["status"] = "downloaded"
            record["actual_size"] = downloaded.stat().st_size
        except Exception as exc:  # noqa: BLE001 - keep status for long download batches.
            record["status"] = "failed"
            record["actual_size"] = destination.stat().st_size if destination.exists() else 0
            record["error"] = str(exc)
            records.append(record)
            write_status(Path(args.status_out), records)
            raise
        records.append(record)
        write_status(Path(args.status_out), records)

    existing_rows = [row for row in slide_rows if destination_for(row, out_dir).exists()]
    print(f"Processed download records: {len(records)}")
    print(f"Slides present under {out_dir}: {len(existing_rows)}/{len(slide_rows)}")
    print(f"Wrote status to {args.status_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
