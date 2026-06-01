#!/usr/bin/env python3
"""Download STAR-count RNA-seq files for a selected TCGA case list."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from gdc_query_tcga_brca import download_file, extract_erbb2_from_star_counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-counts", default="data/tcga_brca_full_query/tcga_brca_star_counts_files.csv")
    parser.add_argument("--cases", default="data/tcga_brca/clinical_her2_cohort_cases.csv")
    parser.add_argument("--out-dir", default="data/tcga_brca/expression_files")
    parser.add_argument("--expression-out", default="data/tcga_brca/erbb2_expression_selected.csv")
    parser.add_argument("--status-out", default="data/tcga_brca/selected_star_counts_download_status.json")
    parser.add_argument("--case-column", default="case_submitter_id")
    parser.add_argument("--sample-type", default="Primary Tumor")
    parser.add_argument("--download-attempts", type=int, default=8)
    return parser.parse_args()


def read_case_ids(path: Path, case_column: str) -> list[str]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        if case_column not in reader.fieldnames:
            raise ValueError(f"{path} does not contain case column {case_column!r}.")
        case_ids = []
        for row in reader:
            case_id = (row.get(case_column) or "").strip()
            if case_id and case_id not in case_ids:
                case_ids.append(case_id)
        return case_ids


def read_star_count_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"file_id", "file_name", "file_size", "case_submitter_id", "sample_type"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def choose_one_star_count_per_case(rows: list[dict[str, str]], case_ids: list[str], sample_type: str) -> list[dict[str, str]]:
    rows_by_case: dict[str, list[dict[str, str]]] = {case_id: [] for case_id in case_ids}
    for row in rows:
        case_id = row.get("case_submitter_id", "")
        if case_id in rows_by_case:
            rows_by_case[case_id].append(row)

    selected = []
    for case_id in case_ids:
        case_rows = rows_by_case.get(case_id, [])
        if not case_rows:
            continue
        case_rows = sorted(
            case_rows,
            key=lambda row: (
                0 if row.get("sample_type") == sample_type else 1,
                row.get("sample_type", ""),
                row.get("file_name", ""),
            ),
        )
        selected.append(case_rows[0])
    return selected


def write_expression(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "case_submitter_id",
        "sample_submitter_id",
        "sample_type",
        "file_id",
        "file_name",
        "gene_id",
        "gene_name",
        "unstranded",
        "tpm_unstranded",
        "fpkm_unstranded",
        "fpkm_uq_unstranded",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_status(path: Path, records: list[dict[str, object]], requested_cases: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "requested_cases": len(requested_cases),
        "records": records,
        "downloaded": sum(1 for record in records if record["status"] == "downloaded"),
        "already_present": sum(1 for record in records if record["status"] == "already_present"),
        "failed": sum(1 for record in records if record["status"] == "failed"),
        "missing_star_count_metadata": len(requested_cases)
        - len({str(record["case_submitter_id"]) for record in records}),
    }
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    case_ids = read_case_ids(Path(args.cases), args.case_column)
    selected_rows = choose_one_star_count_per_case(
        read_star_count_rows(Path(args.star_counts)),
        case_ids,
        args.sample_type,
    )
    out_dir = Path(args.out_dir)
    expression_rows: list[dict[str, str]] = []
    records: list[dict[str, object]] = []

    for index, row in enumerate(selected_rows, start=1):
        case_id = row["case_submitter_id"]
        destination = out_dir / case_id / row["file_name"]
        expected_size = int(float(row["file_size"])) if row.get("file_size") else None
        record = {
            "case_submitter_id": case_id,
            "file_id": row["file_id"],
            "file_name": row["file_name"],
            "destination": str(destination),
            "expected_size": expected_size,
            "status": "",
            "actual_size": 0,
            "error": "",
        }
        try:
            if destination.exists() and destination.stat().st_size > 0:
                record["status"] = "already_present"
                downloaded = destination
            else:
                print(f"[{index}/{len(selected_rows)}] Downloading expression {case_id} {row['file_name']}")
                downloaded = download_file(
                    row["file_id"],
                    destination,
                    expected_size=expected_size,
                    max_attempts=args.download_attempts,
                )
                record["status"] = "downloaded"
            record["actual_size"] = downloaded.stat().st_size
            erbb2 = extract_erbb2_from_star_counts(downloaded)
            expression_rows.append(
                {
                    "case_submitter_id": case_id,
                    "sample_submitter_id": row.get("sample_submitter_id", ""),
                    "sample_type": row.get("sample_type", ""),
                    "file_id": row.get("file_id", ""),
                    "file_name": row.get("file_name", ""),
                    **erbb2,
                }
            )
        except Exception as exc:  # noqa: BLE001 - keep long-batch status.
            record["status"] = "failed"
            record["actual_size"] = destination.stat().st_size if destination.exists() else 0
            record["error"] = str(exc)
            records.append(record)
            write_status(Path(args.status_out), records, case_ids)
            raise
        records.append(record)
        write_status(Path(args.status_out), records, case_ids)

    write_expression(Path(args.expression_out), expression_rows)
    print(f"Wrote ERBB2 expression rows for {len(expression_rows)} cases to {args.expression_out}")
    print(f"Wrote download status to {args.status_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
