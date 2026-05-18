#!/usr/bin/env python3
"""Query GDC TCGA-BRCA slide/RNA-seq files and extract ERBB2 expression."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path
from typing import Any

import requests

GDC_API = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-BRCA"
ERBB2_GENE_ID = "ENSG00000141736"
ERBB2_SYMBOL = "ERBB2"


def eq_filter(field: str, value: str | list[str]) -> dict[str, Any]:
    return {"op": "=", "content": {"field": field, "value": value}}


def and_filter(*parts: dict[str, Any]) -> dict[str, Any]:
    return {"op": "and", "content": list(parts)}


def query_files(filters: dict[str, Any], fields: list[str], page_size: int = 2000) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    offset = 0
    while True:
        payload = {
            "filters": filters,
            "fields": ",".join(fields),
            "format": "JSON",
            "size": page_size,
            "from": offset,
            "sort": "cases.submitter_id:asc,file_name:asc",
        }
        response = requests.post(f"{GDC_API}/files", json=payload, timeout=120)
        response.raise_for_status()
        data = response.json()["data"]
        batch = data["hits"]
        hits.extend(batch)
        pagination = data["pagination"]
        if len(hits) >= pagination["total"] or not batch:
            break
        offset += page_size
    return hits


def first_case(hit: dict[str, Any]) -> dict[str, Any]:
    cases = hit.get("cases") or []
    return cases[0] if cases else {}


def first_sample(hit: dict[str, Any]) -> dict[str, Any]:
    case = first_case(hit)
    samples = case.get("samples") or []
    return samples[0] if samples else {}


def case_submitter_id(hit: dict[str, Any]) -> str:
    return first_case(hit).get("submitter_id", "")


def sample_type(hit: dict[str, Any]) -> str:
    return first_sample(hit).get("sample_type", "")


def case_barcode_from_file_name(file_name: str) -> str:
    match = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", file_name)
    return match.group(1) if match else ""


def limit_to_cases(hits: list[dict[str, Any]], case_limit: int | None) -> list[dict[str, Any]]:
    if not case_limit:
        return hits
    keep: set[str] = set()
    limited: list[dict[str, Any]] = []
    for hit in hits:
        case_id = case_submitter_id(hit) or case_barcode_from_file_name(hit.get("file_name", ""))
        if not case_id:
            continue
        if case_id not in keep and len(keep) >= case_limit:
            continue
        keep.add(case_id)
        limited.append(hit)
    return limited


def select_one_file_per_case(hits: list[dict[str, Any]], preferred_sample_type: str | None) -> list[dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    for hit in hits:
        case_id = case_submitter_id(hit) or case_barcode_from_file_name(hit.get("file_name", ""))
        if not case_id:
            continue
        if case_id not in selected:
            selected[case_id] = hit
            continue
        if preferred_sample_type and sample_type(hit) == preferred_sample_type:
            selected[case_id] = hit
    return list(selected.values())


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def write_manifest(path: Path, hits: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "filename", "md5", "size", "state"], delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(
                {
                    "id": hit.get("file_id", ""),
                    "filename": hit.get("file_name", ""),
                    "md5": hit.get("md5sum", ""),
                    "size": hit.get("file_size", ""),
                    "state": hit.get("state", "released"),
                }
            )


def read_case_file(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if "," in sample:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                return []
            case_field = "case_submitter_id" if "case_submitter_id" in reader.fieldnames else reader.fieldnames[0]
            return [row[case_field] for row in reader if row.get(case_field)]
        return [line.strip() for line in handle if line.strip() and not line.startswith("#")]


def write_file_table(path: Path, hits: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "file_id",
            "file_name",
            "file_size",
            "data_type",
            "data_format",
            "experimental_strategy",
            "case_submitter_id",
            "sample_submitter_id",
            "sample_type",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for hit in hits:
            sample = first_sample(hit)
            writer.writerow(
                {
                    "file_id": hit.get("file_id", ""),
                    "file_name": hit.get("file_name", ""),
                    "file_size": hit.get("file_size", ""),
                    "data_type": hit.get("data_type", ""),
                    "data_format": hit.get("data_format", ""),
                    "experimental_strategy": hit.get("experimental_strategy", ""),
                    "case_submitter_id": case_submitter_id(hit),
                    "sample_submitter_id": sample.get("submitter_id", ""),
                    "sample_type": sample.get("sample_type", ""),
                }
            )


def download_file(file_id: str, destination: Path, expected_size: int | None = None, max_attempts: int = 8) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        if expected_size is None or destination.stat().st_size == expected_size:
            return destination
        destination.unlink()
    tmp_path = destination.with_suffix(destination.suffix + ".part")
    if expected_size and tmp_path.exists() and tmp_path.stat().st_size > expected_size:
        tmp_path.unlink()
    for attempt in range(1, max_attempts + 1):
        resume_from = tmp_path.stat().st_size if tmp_path.exists() else 0
        headers = {"Range": f"bytes={resume_from}-"} if resume_from else None
        mode = "ab" if resume_from else "wb"
        try:
            with requests.get(f"{GDC_API}/data/{file_id}", headers=headers, stream=True, timeout=300) as response:
                if resume_from and response.status_code == 200:
                    mode = "wb"
                    resume_from = 0
                response.raise_for_status()
                with tmp_path.open(mode) as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if expected_size and tmp_path.stat().st_size != expected_size:
                if tmp_path.stat().st_size > expected_size:
                    tmp_path.unlink()
                raise requests.RequestException(
                    f"incomplete download: got {tmp_path.stat().st_size if tmp_path.exists() else 0} bytes, expected {expected_size}"
                )
            tmp_path.replace(destination)
            return destination
        except (requests.RequestException, OSError) as exc:
            if attempt == max_attempts:
                raise
            partial = tmp_path.stat().st_size if tmp_path.exists() else 0
            print(
                f"Download failed for {destination.name} on attempt {attempt}/{max_attempts} "
                f"after {partial} bytes: {exc}. Retrying...",
                file=sys.stderr,
            )
    if destination.exists():
        return destination
    raise FileNotFoundError(destination)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def extract_erbb2_from_star_counts(path: Path) -> dict[str, str]:
    header: list[str] | None = None
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            columns = line.split("\t")
            if columns[0] == "gene_id":
                header = columns
                continue
            if header is None:
                continue
            row = dict(zip(header, columns))
            gene_id = row.get("gene_id", "").split(".")[0]
            gene_name = row.get("gene_name", "")
            if gene_id == ERBB2_GENE_ID or gene_name == ERBB2_SYMBOL:
                return {
                    "gene_id": row.get("gene_id", ""),
                    "gene_name": gene_name,
                    "unstranded": row.get("unstranded", ""),
                    "tpm_unstranded": row.get("tpm_unstranded", ""),
                    "fpkm_unstranded": row.get("fpkm_unstranded", ""),
                    "fpkm_uq_unstranded": row.get("fpkm_uq_unstranded", ""),
                }
    raise ValueError(f"Could not find {ERBB2_SYMBOL}/{ERBB2_GENE_ID} in {path}")


def write_erbb2_expression(path: Path, expression_rows: list[dict[str, str]]) -> None:
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(expression_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="data/tcga_brca", help="Directory for manifests, metadata, and ERBB2 expression.")
    parser.add_argument("--case-limit", type=int, default=None, help="Limit to the first N TCGA cases for a pilot.")
    parser.add_argument("--sample-type", default="Primary Tumor", help="Preferred sample type for one RNA-seq file per case.")
    parser.add_argument("--download-expression", action="store_true", help="Download selected STAR-count files and extract ERBB2.")
    parser.add_argument("--download-slides", action="store_true", help="Directly download selected SVS slides. Use only for small pilots.")
    parser.add_argument("--max-slide-downloads", type=int, default=0, help="Maximum slides to download when --download-slides is set.")
    parser.add_argument("--download-attempts", type=int, default=8, help="Maximum attempts per downloaded file.")
    parser.add_argument(
        "--slide-case-id",
        action="append",
        default=[],
        help="Restrict direct slide downloads to one or more TCGA case submitter IDs. Can be repeated.",
    )
    parser.add_argument(
        "--slide-case-file",
        default=None,
        help="CSV or newline-delimited file of TCGA case submitter IDs for direct slide downloads.",
    )
    parser.add_argument(
        "--one-slide-per-case",
        action="store_true",
        help="For direct SVS downloads, keep only one selected slide per case after ordering.",
    )
    parser.add_argument(
        "--slide-download-order",
        default="smallest",
        choices=["smallest", "manifest"],
        help="Order for direct SVS downloads. The manifest itself is always written in GDC query order.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    common_fields = [
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "state",
        "data_type",
        "data_format",
        "experimental_strategy",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
    ]

    slide_filters = and_filter(
        eq_filter("cases.project.project_id", PROJECT_ID),
        eq_filter("data_type", "Slide Image"),
        eq_filter("data_format", "SVS"),
        eq_filter("access", "open"),
    )
    expression_filters = and_filter(
        eq_filter("cases.project.project_id", PROJECT_ID),
        eq_filter("data_category", "Transcriptome Profiling"),
        eq_filter("data_type", "Gene Expression Quantification"),
        eq_filter("analysis.workflow_type", "STAR - Counts"),
        eq_filter("access", "open"),
    )

    print("Querying GDC diagnostic slide images...", file=sys.stderr)
    slide_hits = limit_to_cases(query_files(slide_filters, common_fields), args.case_limit)
    print(f"Found {len(slide_hits)} slide files after case limit.", file=sys.stderr)
    write_json(out_dir / "file_metadata_slides.json", slide_hits)
    write_manifest(out_dir / "tcga_brca_diagnostic_slides_manifest.tsv", slide_hits)
    write_file_table(out_dir / "tcga_brca_diagnostic_slides_files.csv", slide_hits)

    print("Querying GDC STAR-count RNA-seq files...", file=sys.stderr)
    expression_hits_all = limit_to_cases(query_files(expression_filters, common_fields), args.case_limit)
    expression_hits = select_one_file_per_case(expression_hits_all, args.sample_type)
    print(f"Selected {len(expression_hits)} RNA-seq files for ERBB2 extraction.", file=sys.stderr)
    write_json(out_dir / "file_metadata_star_counts.json", expression_hits)
    write_manifest(out_dir / "tcga_brca_star_counts_manifest.tsv", expression_hits)
    write_file_table(out_dir / "tcga_brca_star_counts_files.csv", expression_hits)

    if args.download_expression:
        expression_rows: list[dict[str, str]] = []
        expression_dir = out_dir / "expression_files"
        for index, hit in enumerate(expression_hits, start=1):
            case_id = case_submitter_id(hit)
            file_name = hit["file_name"]
            destination = expression_dir / case_id / file_name
            print(f"[{index}/{len(expression_hits)}] Downloading expression {case_id} {file_name}", file=sys.stderr)
            downloaded = download_file(
                hit["file_id"],
                destination,
                expected_size=int(hit["file_size"]) if hit.get("file_size") else None,
                max_attempts=args.download_attempts,
            )
            erbb2 = extract_erbb2_from_star_counts(downloaded)
            sample = first_sample(hit)
            expression_rows.append(
                {
                    "case_submitter_id": case_id,
                    "sample_submitter_id": sample.get("submitter_id", ""),
                    "sample_type": sample.get("sample_type", ""),
                    "file_id": hit.get("file_id", ""),
                    "file_name": file_name,
                    **erbb2,
                }
            )
        write_erbb2_expression(out_dir / "erbb2_expression.csv", expression_rows)

    if args.download_slides:
        slide_download_hits = slide_hits
        wanted_cases = set(args.slide_case_id)
        if args.slide_case_file:
            wanted_cases.update(read_case_file(Path(args.slide_case_file)))
        if wanted_cases:
            slide_download_hits = [hit for hit in slide_download_hits if case_submitter_id(hit) in wanted_cases]
        if args.slide_download_order == "smallest":
            slide_download_hits = sorted(slide_download_hits, key=lambda hit: int(hit.get("file_size") or 0))
        if args.one_slide_per_case:
            first_by_case: dict[str, dict[str, Any]] = {}
            for hit in slide_download_hits:
                case_id = case_submitter_id(hit)
                if case_id and case_id not in first_by_case:
                    first_by_case[case_id] = hit
            slide_download_hits = list(first_by_case.values())
        selected_slides = slide_download_hits[: args.max_slide_downloads or len(slide_download_hits)]
        slides_dir = out_dir / "slides"
        for index, hit in enumerate(selected_slides, start=1):
            case_id = case_submitter_id(hit) or case_barcode_from_file_name(hit["file_name"])
            destination = slides_dir / case_id / hit["file_name"]
            print(f"[{index}/{len(selected_slides)}] Downloading slide {case_id} {hit['file_name']}", file=sys.stderr)
            download_file(
                hit["file_id"],
                destination,
                expected_size=int(hit["file_size"]) if hit.get("file_size") else None,
                max_attempts=args.download_attempts,
            )

    print(f"Done. Outputs are in {out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
